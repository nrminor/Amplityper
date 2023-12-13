#!/usr/bin/env Rscript

renv::activate()
renv::restore(prompt = FALSE)

# load some helpful tidyverse packages
suppressPackageStartupMessages(c(
  library(readr),
  library(dplyr),
  library(stringr),
  library(tidyr)
))
renv::snapshot()

# load the input and output file names from the config file
source("config.R")

# set the columns for the data frame that will be made by splitting up the header
columns <- c(
  "gene number",
  "gene",
  "tag",
  "db_xref",
  "location",
  "gbkey"
)

# read the file, process it, and write the bed file in a pipeline
read_delim(gene_fasta, delim = "\t",
           col_names = FALSE, show_col_types = FALSE) %>%
  mutate(`X1` = str_remove_all(`X1`, "\\[|\\]"),
         `X1` = str_remove_all(`X1`, "gene="),
         `X1` = str_remove_all(`X1`, "location=")) %>%
  filter(grepl("gene", `X1`),
         grepl(">", `X1`, fixed = TRUE)) %>%
  separate(`X1`, " ", into = columns, remove = TRUE) %>%
  mutate(`gene number` = str_remove_all(`gene number`, fixed(">lcl|")),
         accession = str_extract(`gene number`, "^[^_]+_[^_]+")) %>%
  select(c(gene, location, accession)) %>%
  separate(location, into = c("start", "stop"),
           remove = TRUE, convert = TRUE) %>%
  select(c(accession, start, stop, gene)) %>%
  write_delim(output_name, delim = "\t", col_names = FALSE)
