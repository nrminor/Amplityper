#!/usr/bin/env python3

"""
Script for visualizing a pileup of amplicon-derived haplotypes with an emphasis on their
mutations from one another.

usage: visualize_haplotypes.py [-h] --sample_id SAMPLE_ID --amplicon AMPLICON
                            --long_table LONG_TABLE --short_table SHORT_TABLE

options:
    -h, --help              show this help message and exit
    --sample_id SAMPLE_ID, -i SAMPLE_ID
                            Sample identifier string.
    --amplicon AMPLICON, -a AMPLICON
                            Amplicon identifier string.
    --long_table LONG_TABLE, -l LONG_TABLE
                            Apache Arrow-formatted table where each row is a single mutation.
    --short_table SHORT_TABLE, -s SHORT_TABLE
                            Excel formatted final Amplityper report with depth-of-coverage data.
"""

import argparse
import asyncio
from pathlib import Path
from typing import List, Tuple

import matplotlib.colors
import numpy as np
import patchworklib as pw  # type: ignore
import polars as pl
from matplotlib import cm  # type: ignore
from plotnine import (
    aes,
    coord_flip,
    element_blank,
    geom_bar,
    geom_tile,
    ggplot,
    labs,
    scale_fill_identity,
    theme,
    theme_minimal,
)


def parse_command_line_args() -> Tuple[str, str, Path, Path]:
    """
        Parse command line arguments while passing errors onto main.

    Args:
        `None`

    Returns:
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sample_id",
        "-i",
        type=str,
        required=True,
        help="Sample identifier string.",
    )
    parser.add_argument(
        "--amplicon",
        "-a",
        type=str,
        required=True,
        default="_LEFT",
        help="Amplicon identifier string.",
    )
    parser.add_argument(
        "--long_table",
        "-l",
        type=Path,
        required=True,
        default="_RIGHT",
        help="Apache Arrow-formatted table where each row is a single mutation.",
    )
    parser.add_argument(
        "--short_table",
        "-s",
        type=Path,
        required=True,
        help="Excel formatted final Amplityper report with depth-of-coverage data.",
    )
    args = parser.parse_args()

    return args.sample_id, args.amplicon, args.long_table, args.short_table


async def read_long_table(long_table_path: Path) -> pl.LazyFrame:
    """
    Asynchronously read the long table while allowing other steps to proceed.
    """
    long_lf = pl.scan_ipc(long_table_path, memory_map=False)

    return long_lf


async def collect_sample_data(
    long_table: pl.LazyFrame, sample_id: str, amplicon: str
) -> pl.LazyFrame:
    """
    Collect the per-sample data that is relevant to the haplotype plot.
    """
    sample_lf = (
        long_table.filter(pl.col("Sample ID") == sample_id)
        .filter(pl.col("Amplicon") == amplicon)
        .with_columns(pl.col("NUC_SUB").str.splitn("-", 3).alias("Position"))
        .unnest("Position")
        .rename(
            {
                "field_0": "Ref",
                "field_1": "Position",
                "field_2": "Alt",
            }
        )
        .with_columns(pl.col("Position").cast(pl.Int32).alias("Position"))
        .with_columns(
            (
                (pl.col("Alt").str.len_chars() > 1) & ~(pl.col("Alt").str.contains("-"))
            ).alias("whether insertion")
        )
        .with_columns(pl.lit(0).alias("insertion offset"))
    )

    return sample_lf


async def adjust_with_offset(accum_df: pl.DataFrame) -> pl.DataFrame:
    """
    Adjust insertion positions with an offset based on their length.
    """

    tmp_df = accum_df.clear()

    for _, data in accum_df.group_by("Contig"):
        new_df = data.with_columns(
            pl.col("insertion offset")
            .cum_sum()
            .cast(pl.Int32)
            .alias("insertion offset")
        )
        tmp_df.vstack(new_df, in_place=True)

    adjusted_df = (
        tmp_df.with_columns(
            (pl.col("Position") + pl.col("insertion offset")).cast(pl.Int32).alias("Position")
        )
        .sort("Position")
        .drop("whether insertion")
    )

    return adjusted_df


async def handle_indels(sample_lf: pl.LazyFrame) -> pl.DataFrame:
    """
    "Explode" out each insertion and deletion so that all positions
    are represented in their own rows.
    """

    mnv_lf = sample_lf.filter(pl.col("Alt").str.len_chars() > 1)
    snv_lf = sample_lf.filter(~(pl.col("Alt").str.len_chars() > 1))

    del_lf = mnv_lf.filter(pl.col("Alt").str.contains("-")).with_columns(
        pl.col("Alt").str.replace("=", "").str.replace("-", "").str.split("").alias("Alt")
    )

    in_lf = (
        mnv_lf.filter(pl.col("whether insertion") == True)  # pylint: disable=C0121
        .with_columns(pl.col("Alt").str.replace("-", "").alias("Alt"))
        .with_columns(pl.col("Alt").str.len_chars().alias("insertion offset"))
        .with_columns(pl.col("Alt").str.split("").alias("Alt"))
    )

    accum_df = snv_lf.collect()

    # deletions
    for _, data in del_lf.collect().group_by("Contig"):
        new_df = (
            data.explode("Alt")
            .filter(pl.col("Alt") != "")
            .filter(pl.col("Alt") != "=")
            .with_columns(
                pl.col("Alt").str.replace("(.*)", "-").alias("Alt")
            )
            .with_row_count()
            .with_columns(
                (pl.col("Position") + pl.col("row_nr")).cast(pl.Int32).alias("Position")
            )
            .drop("row_nr")
        )
        accum_df.vstack(new_df, in_place=True)

    # insertions
    for _, data in in_lf.collect().group_by("Contig"):
        new_df = (
            data.explode("Alt")
            .filter(pl.col("Alt") != "")
            .with_row_count()
            .with_columns(
                (pl.col("Position") + pl.col("row_nr")).cast(pl.Int32).alias("Position")
            )
            .with_columns(
                pl.col("insertion offset").cast(pl.Int32).alias("insertion offset")
            )
            .drop("row_nr")
        )
        accum_df.vstack(new_df, in_place=True)

    accum_df = accum_df.unique()

    adjusted_df = await adjust_with_offset(accum_df)

    return adjusted_df


async def define_all_positions(
    sample_lf: pl.LazyFrame, adjusted_df: pl.DataFrame
) -> Tuple[list[int], pl.DataFrame]:
    """
    Make sure that all positions across the amplicons are represented and
    therefore plot-able.
    """

    ref_start = sample_lf.select("Start Position").unique().collect().item()
    ref_stop = sample_lf.select("Stop Position").unique().collect().item()

    final_df = adjusted_df.select("Position", "Alt", "Contig").clear()

    for haplotype, df in adjusted_df.group_by("Contig"):
        new_stop = ref_stop + df.select("insertion offset").max().item()
        amplicon_positions = pl.int_range(ref_start, new_stop, eager=True).to_list()

        final_df.vstack(
            pl.DataFrame({"Position": amplicon_positions})
            .with_columns(
                pl.col("Position").cast(pl.Int32)
            )
            .join(
                df.select("Position", "Alt", "Contig"),
                how="left",
                on="Position",
            )
            .with_columns(pl.lit(haplotype).alias("Contig"))
            .unique(),
            in_place=True,
        )

    return amplicon_positions, final_df


async def visual_stylings() -> pl.LazyFrame:
    """
    Define plot graphical parameters.
    """

    num_colors = 4
    palette = [cm.viridis(i) for i in np.linspace(0, 1, num_colors)]  # type: ignore # pylint: disable=E1101
    hex_palette = [matplotlib.colors.to_hex(color) for color in palette]

    col_lf = pl.LazyFrame(
        {
            "Alt": ["A", "T", "C", "G", "N", "-"],
            "Color": hex_palette + ["#d3d3d3", "#FFFFFF"],
        }
    )

    return col_lf


async def compile_plotting_data(
    final_df: pl.DataFrame, col_lf: pl.LazyFrame, amplicon_positions: List[int]
):
    """
    Compute plotting dataset for haplotype pileup plot and add colors.
    """

    plotting_lf = (
        pl.LazyFrame({"Position": amplicon_positions})
        .with_columns(
            pl.col("Position").cast(pl.Int32).alias("Position")
        )
        .join(
            final_df.select("Position", "Alt", "Contig").lazy(),
            how="left",
            on="Position",
        )
        .unique()
        .join(col_lf, on="Alt", how="left")
        .with_columns(
            pl.when(pl.col("Color").is_null())
            .then(pl.lit("#d3d3d3"))
            .otherwise(pl.col("Color"))
            .alias("Color")
        )
    )

    return plotting_lf


async def compile_frequency_df(
    depths_df: pl.DataFrame, sample_id: str, amplicon: str
) -> pl.DataFrame:
    """
    Compile frequency information for each haplotype.
    """

    freq_df = (
        depths_df.filter(pl.col("Sample ID") == sample_id)
        .filter(pl.col("Amplicon") == amplicon)
        .with_columns(
            pl.col("Haplotype")
            .str.replace(amplicon, "")
            .str.to_lowercase()
            .str.replace_all(" ", "")
            .alias("Haplotype")
        )
        .select("Haplotype", "Depth of Coverage")
        .with_columns(
            (pl.col("Depth of Coverage") / pl.col("Depth of Coverage").sum()).alias(
                "Frequency"
            )
        )
        .sort("Frequency", descending=True)
    )

    return freq_df


async def create_joint_df(
    plotting_lf: pl.LazyFrame, freq_df: pl.DataFrame
) -> pl.LazyFrame:
    """
    Combine the two dataframes to make plotting straightforward.
    """
    joint_lf = (
        plotting_lf.join(
            freq_df.lazy().select("Haplotype", "Frequency"),
            how="left",
            left_on="Contig",
            right_on="Haplotype",
        )
        .unique()
        .drop_nulls(subset="Frequency")
        .sort("Frequency", descending=True)
    )

    return joint_lf


async def export_plotting_data(
    joint_lf: pl.LazyFrame, sample_id: str, amplicon: str
) -> None:
    """
    Export plotting data as a TSV for transparency and reproducibility.
    """
    joint_lf.sink_csv(
        f"{sample_id}_{amplicon}_plotting_data_with_freqs.tsv", separator="\t"
    )


async def render_haplo_stack(joint_df):
    """
    Render the haplotype stack plot.
    """
    haplo_stack_with_freq = (
        ggplot(joint_df.collect(), aes(x="Position", y="reorder(Contig, Frequency)"))
        + geom_tile(aes(fill="Color", width=0.9, height=0.8))
        + theme_minimal(base_family="helvetica neue")
        + scale_fill_identity()
        + labs(x="Position", y="Haplotypes", fill="Alternate\nAllele")
        + theme(
            legend_position="none",
            figure_size=(15, 10),
            dpi=300,
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
        )
    )

    return haplo_stack_with_freq


async def render_depth_bars(freq_df: pl.DataFrame):
    """
    Render the depth barplot.
    """

    depth_plot = (
        ggplot(freq_df, aes(x="reorder(Haplotype, Frequency)", y="Frequency"))
        + geom_bar(stat="identity")
        + coord_flip()
        + labs(y="Frequency")
        + theme_minimal()
        + theme(
            legend_position="none",
            figure_size=(6, 10),
            dpi=300,
            axis_text_y=element_blank(),
            axis_title_y=element_blank(),
        )
    )

    return depth_plot


async def render_multipanel(haplo_stack_with_freq, depth_plot) -> None:
    """
    Render the multipanel and export as a PDF.
    """
    gg1 = pw.load_ggplot(haplo_stack_with_freq, figsize=(8, 4))
    gg2 = pw.load_ggplot(depth_plot, figsize=(2, 4))
    gg12 = gg1 | gg2
    gg12.savefig(fname="joint_plot.pdf")


async def main() -> None:
    """
    Manage the flow of data through the above functions within an
    asynchronous runtime.
    """

    # parse command line arguments
    sample_id, amplicon, long_table_path, short_table_path = parse_command_line_args()

    # read in the long table and set plot style parameters while waiting
    long_lf = await read_long_table(long_table_path)
    col_lf = await visual_stylings()

    # collect sample data
    sample_lf = await collect_sample_data(long_lf, sample_id, amplicon)

    # handle insertions and deletions so that they are represented as empty
    adjusted_df = await handle_indels(sample_lf)

    # flesh out all amplicon positions for the plotting data
    amp_positions, final_df = await define_all_positions(sample_lf, adjusted_df)

    # compile plotting data for haplotype pileup plot
    plotting_lf = await compile_plotting_data(final_df, col_lf, amp_positions)

    # read the table with depth information
    depths_df = pl.read_excel(short_table_path)

    # compute frequencies for each haplotype and join it with the pileup table
    freq_df = await compile_frequency_df(depths_df, sample_id, amplicon)
    joint_lf = await create_joint_df(plotting_lf, freq_df)

    # export the final plotting data for transparency while allowing the
    # remaining steps to proceed
    await export_plotting_data(joint_lf, sample_id, amplicon)

    # render the haplotype pileup "stack"
    haplo_stack_with_freq = await render_haplo_stack(joint_lf)

    # render the bar plot of frequencies
    depth_plot = await render_depth_bars(freq_df)

    # render and export the multipanel plot
    await render_multipanel(haplo_stack_with_freq, depth_plot)


if __name__ == "__main__":
    asyncio.run(main())
