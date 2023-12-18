#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {

	println()
    println("Note: This workflow currently only supports amplicons sequenced on an Illumina paired-end platform.")
    println("Support for long reads (PacBio and Oxford Nanopore) will be added in the future.")
	println("-----------------------------------------------------------------------------------------------")
	println()

	assert (
		( params.long_reads == true && params.illumina_pe == false ) ||
		( params.long_reads == false && params.illumina_pe == true ) ||
		( params.long_reads == true && params.single_short_reads == false ) ||
		( params.long_reads == false && params.single_short_reads == true ) ||
		( params.illumina_pe == false && params.single_short_reads == true ) ||
		( params.illumina_pe == true && params.single_short_reads == false )
	) : "Please make sure that only one out of the three parameters long_reads, illumina_pe,\nand single_short_reads is set to true to avoid confusing the workflow."

	// input channels
	if ( params.long_reads == true ) {

		ch_reads = Channel
			.fromPath( "${params.fastq_dir}/*.fastq.gz" )
			.map { reads -> tuple( file(reads).getSimpleName(), file(reads) ) }

	} else if ( params.single_short_reads == true ) {

		ch_reads = Channel
			.fromPath( "${params.fastq_dir}/*_001.fastq.gz" )
			.map { reads -> tuple( file(reads).getSimpleName(), file(reads) ) }

	} else {

		ch_reads = Channel
			.fromFilePairs( "${params.fastq_dir}/*{_R1,_R2}_001.fastq.gz", flat: true )

	}

    ch_primer_bed = Channel
        .fromPath( params.primer_bed )

	ch_refseq = Channel
		.fromPath( params.reference )

	ch_refgff = Channel
		.fromPath( params.gff )

	ch_sc2_genes = Channel
		.fromPath( params.sc2_genes )

	ch_r_config = Channel
		.fromPath( params.gene_to_bed_config )

	// Workflow steps
	RESPLICE_PRIMERS (
		ch_primer_bed
	)

	GET_GENE_BED (
		ch_r_config,
		ch_sc2_genes
	)

	CROSS_REF_WITH_GENES (
		RESPLICE_PRIMERS.out,
		GET_GENE_BED.out
	)

    SUBSET_BED_FILE (
        CROSS_REF_WITH_GENES.out
    )

	SPLIT_PRIMER_COMBOS (
		SUBSET_BED_FILE.out.bed
	)

	SPLIT_PRIMER_COMBOS.out.count().view()


	GET_PRIMER_SEQS (
		SPLIT_PRIMER_COMBOS.out.flatten()
	)

	if ( params.illumina_pe == true ) {

		MERGE_PAIRS (
			ch_reads
		)

		CLUMP_READS (
			MERGE_PAIRS.out
		)

	} else {

		CLUMP_READS (
			ch_reads
		)

	}

    FIND_ADAPTER_SEQS (
        CLUMP_READS.out
    )

	TRIM_ADAPTERS (
        FIND_ADAPTER_SEQS.out
	)

    FIND_COMPLETE_AMPLICONS (
        TRIM_ADAPTERS.out,
        GET_PRIMER_SEQS.out
    )

	if ( params.illumina_pe == true || params.single_short_reads == true ) {

		REMOVE_OPTICAL_DUPLICATES (
			FIND_COMPLETE_AMPLICONS.out
		)

		REMOVE_LOW_QUALITY_REGIONS (
			REMOVE_OPTICAL_DUPLICATES.out
		)

		REMOVE_ARTIFACTS (
			REMOVE_LOW_QUALITY_REGIONS.out
		)

		ERROR_CORRECT_PHASE_ONE (
			REMOVE_ARTIFACTS.out
		)

		ERROR_CORRECT_PHASE_TWO (
			ERROR_CORRECT_PHASE_ONE.out
		)

		ERROR_CORRECT_PHASE_THREE (
			ERROR_CORRECT_PHASE_TWO.out
		)

		QUALITY_TRIM (
			ERROR_CORRECT_PHASE_THREE.out
		)

	} else {

		REMOVE_ARTIFACTS (
			FIND_COMPLETE_AMPLICONS.out
		)

		QUALITY_TRIM (
			REMOVE_ARTIFACTS.out
		)

	}

    VALIDATE_SEQS (
        QUALITY_TRIM.out
			.map { name, combo, reads -> tuple( name, combo, file(reads), file(reads).countFastq() ) }
			.filter { it[3] >= params.min_reads }
			.map { name, combo, reads, count -> tuple( name, combo, file(reads) ) }
    )

	FASTQC (
		VALIDATE_SEQS.out
	)

	MULTIQC (
		FASTQC.out.multiqc_data.collect()
	)

	IDENTIFY_HAPLOTYPES (
		VALIDATE_SEQS.out
	)

	NAME_HAPLOTYPES (
		IDENTIFY_HAPLOTYPES.out.deduped_fasta
	)

    MAP_HAPLOTYPE_TO_REF (
        NAME_HAPLOTYPES.out
			.splitFasta( by: 1, file: true ),
		ch_refseq
    )

    CALL_HAPLOTYPE_VARIANTS (
        MAP_HAPLOTYPE_TO_REF.out,
		ch_refseq
    )

	RUN_SNPEFF (
		CALL_HAPLOTYPE_VARIANTS.out,
		ch_refseq
	)

	GENERATE_TIDY_VCF (
		RUN_SNPEFF.out
	)

	GENERATE_IVAR_TABLE (
        MAP_HAPLOTYPE_TO_REF.out,
		ch_refseq,
		ch_refgff
	)

	// GENERATE_FINAL_REPORT (
	// 	GENERATE_TIDY_VCF.out.collect(),
	// 	GENERATE_IVAR_TABLE.out.collect(),
	// 	CROSS_REF_WITH_GENES.out
	// )

}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config

// whether to terminate when unrecoverable error occurs
if ( params.debugmode ) {
	params.errorMode = 'terminate'
	params.pubMode = 'copy'
} else {
	params.errorMode = 'ignore'
	params.pubMode = 'symlink'
}

// overarching first level in results file hierarchy
params.amplicon_results = params.results + "/amplicon_${params.desired_amplicon}"

// Preprocessing results subdirectories
params.preprocessing = params.amplicon_results + "/01_preprocessing"
params.merged_reads = params.preprocessing + "/01_merged_pairs"
params.clumped_reads = params.preprocessing + "/02_clumped_reads"
params.trim_adapters = params.preprocessing + "/03_trim_adapters"
params.amplicon_reads = params.preprocessing + "/04_amplicon_reads"
params.optical_dedupe = params.preprocessing + "/05_optical_dedup"
params.low_quality = params.preprocessing + "/06_remove_low_quality"
params.remove_artifacts = params.preprocessing + "/07_remove_artifacts"
params.error_correct = params.preprocessing + "/08_error_correct"
params.qtrim = params.preprocessing + "/09_quality_trim"
params.clipped = params.preprocessing + "/10_clipped_reads"
params.complete_amplicon = params.preprocessing + "/11_complete amplicons"
params.fastqc_results = params.preprocessing + "/12_FastQC_results"

// haplotyping results
params.haplotyping_results = params.amplicon_results + "/02_haplotyping_results"
params.haplotypes = params.haplotyping_results + "/01_haplotypes"
params.variants = params.haplotyping_results + "/02_contig_VCFs"
params.ivar_tables = params.haplotyping_results + "/03_ivar_tables"


// --------------------------------------------------------------- //




// PROCESS SPECIFICATION
// --------------------------------------------------------------- //

process RESPLICE_PRIMERS {

    /*
    */

	tag "${params.desired_amplicon}"
    label "general"

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

    input:
    path bed_file

    output:
    path "respliced.bed"

    script:
    """
	resplice_primers.py -i ${bed_file}
    """

}

process GET_GENE_BED {

	/*
    */

	tag "${params.desired_amplicon}"
    label "general"

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

    input:
    path config
	path fasta

    output:
    path "*.bed"

    script:
    """
	ncbi_fasta_to_bed.R
    """

}

process CROSS_REF_WITH_GENES {

	/*
    */

	tag "${params.desired_amplicon}"
    label "general"

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

    input:
    path primer_bed
	path gene_bed

    output:
    path "${bed_name}_with_genes.bed"

    script:
	bed_name = file(primer_bed.toString()).getSimpleName()
    """
	ref=`csvtk replace -t ${primer_bed} -f 1 -p " " -r "" | cut -f 1 | uniq` && \
	unmatched=`csvtk replace -t ${gene_bed} -f 1 -p " " -r "" | cut -f 1 | uniq` && \
	csvtk replace -H -t ${gene_bed} -p "\$unmatched" -r "\$ref" -o corrected.bed && \
	bedtools intersect -a ${primer_bed} -b corrected.bed -wb | \
	csvtk cut -f 1,2,3,4,5,6,11 -t > ${bed_name}_with_genes.bed
    """

}

process SUBSET_BED_FILE {

    /*
    */

	tag "${params.desired_amplicon}"
    label "general"

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

    input:
    path bed_file

    output:
    path "all_primer_combos.bed", emit: bed

    script:
    """
    grep ${params.desired_amplicon} ${bed_file} > all_primer_combos.bed
    """

}

process SPLIT_PRIMER_COMBOS {

    /*
    */

	tag "${params.desired_amplicon}"
    label "general"

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	path all_primer_combos

	output:
	path "${params.desired_amplicon}*.bed"

	script:
	"""
	split_primer_combos.py -i ${all_primer_combos}
	"""

}

process GET_PRIMER_SEQS {

    /*
    */

	tag "${params.desired_amplicon}"
    label "general"

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	path bed

	output:
	path "${primer_combo}_patterns.txt"

	script:
	primer_combo = file(bed.toString()).getSimpleName()
	"""
	bedtools getfasta -fi ${params.reference} -bed ${bed} | \
	grep -v "^>" > ${primer_combo}_patterns.txt
	"""
}

process MERGE_PAIRS {

    /* */

	tag "${sample_id}"
    label "general"
	publishDir params.merged_reads, mode: params.pubMode, overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 4

    input:
	tuple val(sample_id), path(reads1), path(reads2)

    output:
    tuple val(sample_id), path("${sample_id}_merged.fastq.gz")

    script:
    """
    bbmerge-auto.sh in1=`realpath ${reads1}` \
	in2=`realpath ${reads2}` \
    out=${sample_id}_merged.fastq.gz \
    outu=${sample_id}_unmerged.fastq.gz \
    strict k=93 extend2=80 rem ordered \
    ihist=${sample_id}_ihist_merge.txt \
    threads=${task.cpus}
    """

}

process CLUMP_READS {

    /*
    */

	tag "${sample_id}"
    label "general"
	publishDir params.clumped_reads, mode: params.pubMode, overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 4

	input:
	tuple val(sample_id), path(reads)

	output:
    tuple val(sample_id), path("${sample_id}_clumped.fastq.gz")

	script:
	"""
	clumpify.sh in=${reads} out=${sample_id}_clumped.fastq.gz t=${task.cpus} reorder
	"""
}

process FIND_ADAPTER_SEQS {

	/* */

	tag "${sample_id}"
    label "general"

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 2

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path(reads), path("${sample_id}_adapters.fasta")

	script:
	"""
    bbmerge.sh in=`realpath ${reads}` outa="${sample_id}_adapters.fasta" # ow qin=33
	"""

}

process TRIM_ADAPTERS {

	/* */

	tag "${sample_id}"
    label "general"
	publishDir params.trim_adapters, mode: params.pubMode, overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 2

	input:
	tuple val(sample_id), path(reads), path(adapters)

    output:
    tuple val(sample_id), path("${sample_id}_no_adapters.fastq.gz")

    script:
    """
	reformat.sh in=`realpath ${reads}` \
	out=${sample_id}_no_adapters.fastq.gz \
	ref=`realpath ${adapters}` \
	uniquenames=t overwrite=true t=${task.cpus} -Xmx8g
    """

}

process FIND_COMPLETE_AMPLICONS {

    /*
    */

    tag "${sample_id}"
    label "general"
	publishDir params.amplicon_reads, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 4

    input:
	tuple val(sample_id), path(reads)
    each path(search_patterns)

    output:
    tuple val(sample_id), val(primer_combo), path("${sample_id}_${primer_combo}_amplicons.fastq.gz")

    script:
	primer_combo = file(search_patterns.toString()).getSimpleName().replace("_patterns", "")
    """
	cat ${reads} | \
    seqkit grep \
	--threads 2 \
	--max-mismatch 3 \
	--by-seq \
	--pattern `head -n 1 ${search_patterns}` | \
	seqkit grep \
	--threads 2 \
	--max-mismatch 3 \
	--by-seq \
	--pattern `tail -n 1 ${search_patterns}` \
    -o ${sample_id}_${primer_combo}_amplicons.fastq.gz
    """

}

process REMOVE_OPTICAL_DUPLICATES {

	/*
	This process removes optical duplicates from the Illumina flow cell.
	*/

	tag "${sample_id}"
    label "general"
	publishDir params.optical_dedupe, pattern: "*.fastq.gz", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 1
	memory 3.GB

	input:
	tuple val(sample_id), val(primer_combo), path(reads)

	output:
	tuple val(sample_id), val(primer_combo), path("${sample_id}_${primer_combo}_deduped.fastq.gz")

	script:
	if ( sample_id.toString().contains("SRR") )
		"""
		rename.sh \
		in=${reads} \
		out=${sample_id}_${primer_combo}_deduped.fastq.gz \
		addpairnum=t
		"""
	else
		"""
		clumpify.sh -Xmx2g in=`realpath ${reads}` \
		out=${sample_id}_${primer_combo}_deduped.fastq.gz \
		threads=${task.cpus} \
		optical tossbrokenreads reorder
		"""

}

process REMOVE_LOW_QUALITY_REGIONS {

	/*
	Low quality regions of each read are removed in this process.
	*/

	tag "${sample_id}"
    label "general"
	publishDir params.low_quality, pattern: "*.fastq.gz", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 1
	memory 3.GB

	input:
	tuple val(sample_id), val(primer_combo), path(reads)

	output:
	tuple val(sample_id), val(primer_combo), path("${sample_id}_${primer_combo}_filtered.fastq.gz")

	script:
	if ( sample_id.toString().contains("SRR") )
		"""
		clumpify.sh -Xmx2g in=`realpath ${reads}` \
		out=${sample_id}_${primer_combo}_filtered.fastq.gz \
		threads=${task.cpus} \
		reorder markduplicates
		"""
	else
		"""
		filterbytile.sh -Xmx2g in=`realpath ${reads}` \
		out=${sample_id}_${primer_combo}_filtered.fastq.gz \
		threads=${task.cpus}
		"""

}

process REMOVE_ARTIFACTS {

	/*
	Here we remove various contantimants that may have ended up in the reads,
	such as PhiX sequences that are often used as a sequencing control.
	*/

	tag "${sample_id}"
    label "general"
	publishDir params.remove_artifacts, pattern: "*.fastq.gz", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 1
	memory 3.GB

	input:
	tuple val(sample_id), val(primer_combo), path(reads)

	output:
	tuple val(sample_id), val(primer_combo), path("${sample_id}_${primer_combo}_remove_artifacts.fastq.gz")

	script:
	"""
	bbduk.sh -Xmx2g in=`realpath ${reads}` \
	out=${sample_id}_${primer_combo}_remove_artifacts.fastq.gz \
	k=31 ref=artifacts,phix ordered cardinality \
	threads=${task.cpus}
	"""

}

process ERROR_CORRECT_PHASE_ONE {

	/*
	Bbmap recommends three phases of read error correction, the first of which
	goes through BBMerge.
	*/

	tag "${sample_id}"
    label "general"
	publishDir params.error_correct, pattern: "*.fastq.gz", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 1
	memory 3.GB

	input:
	tuple val(sample_id), val(primer_combo), path(reads)

	output:
	tuple val(sample_id), val(primer_combo), path("${sample_id}_${primer_combo}_error_correct1.fastq.gz")

	script:
	"""
	bbmerge.sh -Xmx2g in=`realpath ${reads}` \
	out=${sample_id}_${primer_combo}_error_correct1.fastq.gz \
	ecco mix vstrict ordered \
	ihist=${sample_id}_${primer_combo}_ihist_merge1.txt \
	threads=${task.cpus}
	"""

}

process ERROR_CORRECT_PHASE_TWO {

	/*
	The second phase of error correction goes through clumpify.sh
	*/

	tag "${sample_id}"
    label "general"
	publishDir params.error_correct, pattern: "*.fastq.gz", mode: 'copy', overwrite: true

	cpus 1
	memory 3.GB

	input:
	tuple val(sample_id), val(primer_combo), path(reads)

	output:
	tuple val(sample_id), val(primer_combo), path("${sample_id}_${primer_combo}_error_correct2.fastq.gz")

	script:
	"""
	clumpify.sh -Xmx2g in=`realpath ${reads}` \
	out=${sample_id}_${primer_combo}_error_correct2.fastq.gz \
	ecc passes=4 reorder \
	threads=${task.cpus} \
	-Xmx2g
	"""

}

process ERROR_CORRECT_PHASE_THREE {

	/*
	The third phase of error correction uses tadpole.sh.
	*/

	tag "${sample_id}"
    label "general"
	publishDir params.error_correct, pattern: "*.fastq.gz", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 1
	memory 3.GB

	input:
	tuple val(sample_id), val(primer_combo), path(reads)

	output:
	tuple val(sample_id), val(primer_combo), path("${sample_id}_${primer_combo}_error_correct3.fastq.gz")

	script:
	"""
	tadpole.sh -Xmx2g in=`realpath ${reads}` \
	out=${sample_id}_${primer_combo}_error_correct3.fastq.gz \
	ecc k=62 ordered \
	threads=${task.cpus}
	"""

}

process QUALITY_TRIM {

	/*
	Here we quality trim reads from both ends to a minimum Phred quality of 10,
	and enforce a minimum read length of 70 bases.
	*/

	tag "${sample_id}"
    label "general"
	publishDir params.qtrim, pattern: "*.fastq.gz", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 1
	memory 3.GB

	input:
	tuple val(sample_id), val(primer_combo), path(reads)

	output:
	tuple val(sample_id), val(primer_combo), path("${sample_id}_${primer_combo}_qtrimmed.fastq.gz")

	script:
	"""
	bbduk.sh -Xmx2g in=`realpath ${reads}` \
	out=${sample_id}_${primer_combo}_qtrimmed.fastq.gz \
	qtrim=rl trimq=10 minlen=70 ordered \
	threads=${task.cpus}
	"""

}

process VALIDATE_SEQS {

    /*
    */

	tag "${sample_id}"
    label "general"
	publishDir params.complete_amplicon, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 1

	input:
	tuple val(sample_id), val(primer_combo), path(reads)

	output:
	tuple val(sample_id), val(primer_combo), path("*.fastq.gz")

	script:
	"""
	cat ${reads} | \
    seqkit seq \
	--threads ${task.cpus} \
	--remove-gaps \
	--validate-seq \
    -o ${sample_id}_${primer_combo}_filtered.fastq.gz
	"""
}

process FASTQC {

    /*
    */

	tag "${sample_id}"
    label "general"
	publishDir params.fastqc_results, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 1

	input:
	tuple val(sample_id), val(primer_combo), path(reads)

	output:
	path "${sample_id}_${primer_combo}_qc.html", emit: html
	path "${sample_id}_${primer_combo}/", emit: multiqc_data

	script:
	"""
	fqc -q ${reads} -s . > ${sample_id}_${primer_combo}_qc.html
	mkdir ${sample_id}_${primer_combo}
	mv fastqc_data.txt ${sample_id}_${primer_combo}/fastqc_data.txt
	"""

}

process MULTIQC {

    /*
    */

	tag "${params.desired_amplicon}"
    label "multiqc"
	publishDir params.preprocessing, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 1

	input:
	path fastqc_files

	output:
	path("*.html")

	script:
	"""
	multiqc ${fastqc_files}
	"""

}

process IDENTIFY_HAPLOTYPES {

	/*
    */

	tag "${sample_id}"
    label "general"

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 1

	input:
	tuple val(sample_id), val(primer_combo), path(reads)

	output:
	tuple val(sample_id), val(primer_combo), path("${sample_id}_deduped.fasta"), emit: deduped_fasta
	path "${sample_id}_${primer_combo}_haplotype_metadata.tsv", emit: metadata

	script:
	"""
	vsearch \
	--fastx_uniques ${reads} \
	--fastaout tmp.fasta \
	--sizeout \
	--minuniquesize ${params.min_reads} \
	--tabbedout tmp.tsv \
	--strand both

	seqkit replace -p ";" -r " " tmp.fasta -o ${sample_id}_deduped.fasta && \
	rm tmp.fasta

	csvtk add-header -t \
	--names orig_label,clust_label,clust_index,seq_index_in_clust,clust_abundance,first_seq_label \
	tmp.tsv | \
	csvtk filter -t --filter "clust_abundance>=${params.min_reads}" \
	> ${sample_id}_${primer_combo}_haplotype_metadata.tsv

	rm tmp.tsv
	"""

}

process NAME_HAPLOTYPES {

	/*
    */

	tag "${sample_id}"
    label "general"
	publishDir "${params.haplotypes}/${sample_id}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 1

	input:
	tuple val(sample_id), val(primer_combo), path(fasta)

	output:
	path "${sample_id}_${primer_combo}_haplotypes.fasta"

	script:
	"""
	seqkit replace -p ";" -r " " ${fasta} | \
	seqkit replace -p "^." -r "${sample_id}_${primer_combo}_haplotype " \
	| seqkit rename --rename-1st-rec \
	-o ${sample_id}_${primer_combo}_haplotypes.fasta
	"""

}

process MAP_HAPLOTYPE_TO_REF {

    /*
    */

    label "general"

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

    cpus 2

	input:
	path haplotype
	each path(refseq)

	output:
	tuple env(hap_name), path("*.bam")

	script:
	"""
	hap_name=`seqkit seq --name --only-id ${haplotype}` && \
	bbmap.sh int=f ref=${refseq} in=`realpath ${haplotype}` out=stdout.sam maxindel=200 | \
	reformat.sh in=stdin.sam out="\${hap_name}.bam"
	"""

}

process CALL_HAPLOTYPE_VARIANTS {

    /*
    */

	tag "${name}"
    label "general"
	publishDir "${params.variants}/${sample_id}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

    cpus 2

	input:
	tuple val(name), path(bam)
	each path(refseq)

	output:
	tuple val(name), path("${name}.vcf"), val(sample_id)

	script:
	sample_id = name.toString().split("_")[0]
	"""
    bcftools mpileup -Ou -f ${refseq} ${bam} | bcftools call --ploidy 1 -mv -Ov -o ${name}.vcf
	"""

}

process RUN_SNPEFF {

    /*
    */

	tag "${name}"
    label "general"
	publishDir "${params.variants}/${sample_id}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

    cpus 1

	input:
	tuple val(name), path(vcf), val(sample_id)
	each path(refseq)

	output:
	tuple val(name), path("${name}_annotated.vcf"), val(sample_id)

	script:
	"""
	REF=`seqkit seq --name --only-id ${refseq}` && \
	snpeff -v \$REF ${vcf} > ${name}_annotated.vcf
	"""

}

process GENERATE_TIDY_VCF {

    /*
    */

	tag "${name}"
    label "general"
	publishDir "${params.variants}/${sample_id}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

    cpus 1

	input:
	tuple val(name), path(vcf), val(sample_id)

	output:
	tuple val(name), path("${name}_annotated.tvcf.tsv"), val(sample_id)

	script:
	"""
	tidyvcf -i ${name}_annotated.vcf -s -o ${name}_annotated.tvcf.tsv
	"""

}

process GENERATE_IVAR_TABLE {

	/* */


	tag "${name}"
    label "iVar"
	publishDir "${params.ivar_tables}/${sample_id}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	input:
	tuple val(name), path(bam)
	each path(refseq)
	each path(refgff)

	output:
	path "*.tsv"

	script:
	sample_id = name.toString().split("_")[0]
	"""
	samtools mpileup -aa -A -d 0 -B -Q 0 --reference ${refseq} ${bam} | \
	ivar variants -p ${name} -t 0 -r ${refseq} -g ${refgff}
	"""

}

process GENERATE_FINAL_REPORT {

	/* */

	tag "${sample_id}"
	publishDir params.amplicon_results, mode: 'copy'

	input:
	path tvcf_files
	path ivar_tables
	each path(gene_bed)

	output:
	path "final_report.xlsx", emit: report_xlsx
	path "*.arrow", emit: arrow_data

	script:
	"""
	read_zap_report.py \
	--results_dir ${params.amplicon_results} \
	--gene_bed ${} \
	--config ${params.reporting_config}
	"""

}

// --------------------------------------------------------------- //
