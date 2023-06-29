#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	
	// input channels
	
	
	// Workflow steps
    CLUMP_READS ()

    TRIM_TO_AMPLICONS (
        CLUMP_READS.out
    )

    MERGE_PAIRS (
        TRIM_TO_AMPLICONS.out
    )

    MAP_TO_REF (
        MERGE_PAIRS.out
    )

    EXTRACT_AMPLICON (
        MAP_TO_REF.out
    )

    BAM_TO_FASTQ (
        EXTRACT_AMPLICON.out
    )

    FILTER_BY_LENGTH (
        BAM_TO_FASTQ.out
    )

    HAPLOTYPE_ASSEMBLY (
        FILTER_BY_LENGTH.out
    )

    // RECORD_FREQUENCIES (
    //     HAPLOTYPE_ASSEMBLY.out
    // )

    DOWNSAMPLE_ASSEMBLIES (
        HAPLOTYPE_ASSEMBLY.out
    )

    MAP_TO_REF2 (
        DOWNSAMPLE_ASSEMBLIES.out
    )

    CALL_CONSENSUS_SEQS (
        MAP_TO_REF2.out
    )

    CALL_VARIANTS (
        MAP_TO_REF2.out
    )

    GENERATE_REPORT (
        DOWNSAMPLE_ASSEMBLIES.out,
        MAP_TO_REF2.out,
        CALL_CONSENSUS_SEQS.out,
        CALL_VARIANTS.out
    )
	
	
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process CLUMP_READS {
	
	// This process does something described here
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	
	
	output:
	
	
	when:
	
	
	script:
	"""
	
	"""
}

process PROCESS_NAME {
	
	// This process does something described here
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	
	
	output:
	
	
	when:
	
	
	script:
	"""
	
	"""
}

process PROCESS_NAME {
	
	// This process does something described here
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	
	
	output:
	
	
	when:
	
	
	script:
	"""
	
	"""
}

process PROCESS_NAME {
	
	// This process does something described here
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	
	
	output:
	
	
	when:
	
	
	script:
	"""
	
	"""
}

process PROCESS_NAME {
	
	// This process does something described here
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	
	
	output:
	
	
	when:
	
	
	script:
	"""
	
	"""
}

process PROCESS_NAME {
	
	// This process does something described here
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	
	
	output:
	
	
	when:
	
	
	script:
	"""
	
	"""
}

process PROCESS_NAME {
	
	// This process does something described here
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	
	
	output:
	
	
	when:
	
	
	script:
	"""
	
	"""
}

process PROCESS_NAME {
	
	// This process does something described here
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	
	
	output:
	
	
	when:
	
	
	script:
	"""
	
	"""
}

process PROCESS_NAME {
	
	// This process does something described here
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	
	
	output:
	
	
	when:
	
	
	script:
	"""
	
	"""
}

process PROCESS_NAME {
	
	// This process does something described here
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	
	
	output:
	
	
	when:
	
	
	script:
	"""
	
	"""
}

process PROCESS_NAME {
	
	// This process does something described here
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	
	
	output:
	
	
	when:
	
	
	script:
	"""
	
	"""
}

process PROCESS_NAME {
	
	// This process does something described here
	
	tag "${tag}"
	publishDir params.results, mode: 'copy'
	
	memory 1.GB
	cpus 1
	time '10minutes'
	
	input:
	
	
	output:
	
	
	when:
	
	
	script:
	"""
	
	"""
}

// --------------------------------------------------------------- //