#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	
	// input channels
	ch_reads = Channel
        .fromFilePairs( "${params.fastq_dir}/*{_R1,_R2}_001.fastq.gz", flat: true )
	
	// Workflow steps
    TRIM_TO_AMPLICONS (
        ch_reads
    )

    CLUMP_READS (
        TRIM_TO_AMPLICONS.out
    )

    // MERGE_PAIRS (
    //     CLUMP_READS.out
    // )

    // MAP_TO_REF (
    //     MERGE_PAIRS.out
    // )

    // EXTRACT_AMPLICON (
    //     MAP_TO_REF.out
    // )

    // BAM_TO_FASTQ (
    //     EXTRACT_AMPLICON.out
    // )

    // FILTER_BY_LENGTH (
    //     BAM_TO_FASTQ.out
    // )

    // HAPLOTYPE_ASSEMBLY (
    //     FILTER_BY_LENGTH.out
    // )

    // // RECORD_FREQUENCIES (
    // //     HAPLOTYPE_ASSEMBLY.out
    // // )

    // DOWNSAMPLE_ASSEMBLIES (
    //     HAPLOTYPE_ASSEMBLY.out
    // )

    // MAP_TO_REF2 (
    //     DOWNSAMPLE_ASSEMBLIES.out
    // )

    // CALL_CONSENSUS_SEQS (
    //     MAP_TO_REF2.out
    // )

    // CALL_VARIANTS (
    //     MAP_TO_REF2.out
    // )

    // GENERATE_REPORT (
    //     DOWNSAMPLE_ASSEMBLIES.out,
    //     MAP_TO_REF2.out,
    //     CALL_CONSENSUS_SEQS.out,
    //     CALL_VARIANTS.out
    // )
	
	
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process TRIM_TO_AMPLICONS {
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'

	cpus 1
	
	input:
	tuple val(sample_id), path(reads)
	
	output:
    tuple val(sample_id), path("*.fastq.gz")
	
	script:
	"""
	ampl-BBDuk.py 
	"""
}

process CLUMP_READS {
	
	tag "${sample_id}"
	publishDir params.results, mode: 'copy'

    cpus 3
	
	input:
	tuple val(sample_id), path(reads1)
	
	output:
    tuple val(sample_id), path("*.fastq.gz")
	
	script:
	"""
	clumpify.sh in=${reads1} out=${sample_id}_clumped.fastq.gz reorder
	"""
}

// process PROCESS_NAME {
	
// 	// This process does something described here
	
// 	tag "${tag}"
// 	publishDir params.results, mode: 'copy'
	
// 	memory 1.GB
// 	cpus 1
// 	time '10minutes'
	
// 	input:
	
	
// 	output:
	
	
// 	when:
	
	
// 	script:
// 	"""
	
// 	"""
// }

// --------------------------------------------------------------- //