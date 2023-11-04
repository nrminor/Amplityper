#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {

    println("Note: This workflow currently only supports amplicons sequenced on a Illumina paired-end platform.")
    println("Support for longer reads will be added in the future.")
	
	// input channels
	ch_reads = Channel
        .fromFilePairs( "${params.fastq_dir}/*{_R1,_R2}_001.fastq.gz", flat: true )

    ch_primer_bed = Channel
        .fromPath( params.primer_bed )
	
	// Workflow steps
    MERGE_PAIRS (
        ch_reads
    )

    FIND_ADAPTER_SEQS (
        MERGE_PAIRS.out
    )

	TRIM_ADAPTERS (
        FIND_ADAPTER_SEQS.out
	)

    CLUMP_READS (
        TRIM_ADAPTERS.out
    )

    GET_PRIMER_SEQS (
        ch_primer_bed
    )

    FIND_COMPLETE_AMPLICONS (
        CLUMP_READS.out,
        GET_PRIMER_SEQS.out.txt
    )

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

    EXTRACT_REF_AMPLICON (
        GET_PRIMER_SEQS.out.txt
    )

    MAP_TO_AMPLICON (
        QUALITY_TRIM.out,
        EXTRACT_REF_AMPLICON.out.seq
    )

    CLIP_AMPLICONS (
        MAP_TO_AMPLICON.out,
        GET_PRIMER_SEQS.out.bed
    )

    BAM_TO_FASTQ (
        EXTRACT_AMPLICON.out
    )

    FILTER_BY_LENGTH (
        EXTRACT_REF_AMPLICON.out.length,
        BAM_TO_FASTQ.out
    )

    // HAPLOTYPE_ASSEMBLY (
    //     FILTER_BY_LENGTH.out
    // )

    // RECORD_FREQUENCIES (
    //     HAPLOTYPE_ASSEMBLY.out
    // )

    DOWNSAMPLE_ASSEMBLIES (
        HAPLOTYPE_ASSEMBLY.out
    )

    MAP_CONTIGS_TO_REF (
        DOWNSAMPLE_ASSEMBLIES.out
    )

    CALL_VARIANTS (
        MAP_CONTIGS_TO_REF.out
    )

    // GENERATE_REPORT (
    //     DOWNSAMPLE_ASSEMBLIES.out,
    //     MAP_CONTIGS_TO_REF.out,
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

process MERGE_PAIRS {
    /* */

    label "general"

    input:
	tuple val(sample_id), path(reads1), path(reads2)

    output:
    tuple val(sample_id), path("${sample_id}_merged.fastq.gz")

    script:
    """
    bbmerge-auto.sh in=`realpath ${reads}` \
    out=${sample}_merged.fastq.gz \
    outu=${sample}_unmerged.fastq.gz \
    strict k=93 extend2=80 rem ordered \
    ihist=${sample}_ihist_merge.txt \
    threads=${task.cpus}
    """

}

process FIND_ADAPTER_SEQS {
	
	/* */
	
	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2
	
	input:
	tuple val(sample_id), path(reads)
	
	output:
	tuple val(sample_id), path(reads), path("${sample_id}_adapters.fasta")
	
	script:
	"""
    bbmerge.sh in=`realpath ${reads}` outa="${sample_id}_adapters.fasta" ow qin=33
	"""

}

process TRIM_ADAPTERS {
	
	/* */
	
	tag "${sample_id}"

	errorStrategy { task.attempt < 3 ? 'retry' : errorMode }
	maxRetries 2
	
	input:
	tuple val(sample_id), path(reads), path(adapters)

    output:
    tuple val(sample_id), path("${sample_id}_no_adapters.fastq.gz")

    script:
    """
	reformat.sh in=`realpath ${reads}` \
	out=${sample_id}_no_adapters.fastq.gz \
	ref=`realpath ${adapters}` \
	uniquenames=t overwrite=true t=${task.cpus}
    """

}

process CLUMP_READS {

    /*
    */
	
	tag "${sample_id}"
    label "general"
	publishDir params.results, mode: 'copy'

    cpus params.available_cpus
	
	input:
	tuple val(sample_id), path(reads)
	
	output:
    tuple val(sample_id), path("${sample_id}_clumped.fastq.gz")
	
	script:
	"""
	clumpify.sh in=${reads} out=${sample_id}_clumped.fastq.gz t=${task.cpus} reorder
	"""
}

process GET_PRIMER_SEQS {

    /*
    */

    label "general"

    input:
    path bed_file

    output:
    path "primer_seqs.bed", emit: bed
    path "primer_seqs.txt", emit: txt

    script:
    """
    grep ${params.desired_amplicon} ${params.primer_bed} > amplicon-primers.bed && \
    bedtools getfasta -fi ${params.reference} -bed amplicon-primers.bed | \
    seqkit fx2tab --no-qual -o primer_seqs.txt
    """

}

process FIND_COMPLETE_AMPLICONS {

    /*
    */

    tag "${sample_id}"
    label "general"

    input:
	tuple val(sample_id), path(reads)
    path primer_seqs
    
    output:
    tuple val(sample_id), path("${sample_id}_amplicons.fastq.gz")

    script:
    """
    seqkit amplicon \
    --primer-file ${primer_seqs} \
    --max-mismatch 2 \
    ${reads} \
    -o ${sample_id}_amplicons.fastq.gz
    """

}

process REMOVE_OPTICAL_DUPLICATES {

	/* 
	This process removes optical duplicates from the Illumina flow cell.
	*/

	tag "${sample_id}"
    label "general"
	// publishDir params.optical_dedupe, pattern: "*.fastq.gz", mode: params.publishMode, overwrite: true

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id) path("${sample_id}_deduped.fastq.gz")

	script:
	"""
	clumpify.sh in=`realpath ${reads}` \
	out=${sample_id}_deduped.fastq.gz \
	threads=${task.cpus} \
	dedupe optical tossbrokenreads
	"""

}

process REMOVE_LOW_QUALITY_REGIONS {

	/* 
	Low quality regions of each read are removed in this process.
	*/

	tag "${sample_id}"
    label "general"
	// publishDir params.low_quality, pattern: "*.fastq.gz", mode: params.publishMode, overwrite: true

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_filtered_by_tile.fastq.gz")

	script:
	"""
	filterbytile.sh in=`realpath ${reads}` \
	out=${sample_id}_filtered_by_tile.fastq.gz \
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
	// publishDir params.remove_artifacts, pattern: "*.fastq.gz", mode: params.publishMode, overwrite: true

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample}_remove_artifacts.fastq.g")

	script:
	"""
	bbduk.sh in=`realpath ${reads}` \
	out=${sample_id}_remove_artifacts.fastq.gz \
	k=31 ref=artifacts,phix ordered cardinality \
	threads=${task.cpus}
	"""

}

process ERROR_CORRECT_PHASE_ONE {

	/* 
	Bbmap recommends three phases of read error correction, the first of which
	goes through BBMerge.
	*/

	tag "${sample}"
    label "general"
	// publishDir params.error_correct, pattern: "*.fastq.gz", mode: params.publishMode, overwrite: true

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_error_correct1.fastq.gz")

	script:
	"""
	bbmerge.sh in=`realpath ${reads}` \
	out=${sample_id}_error_correct1.fastq.gz \
	ecco mix vstrict ordered \
	ihist=${sample_id}_ihist_merge1.txt \
	threads=${task.cpus}
	"""

}

process ERROR_CORRECT_PHASE_TWO {

	/* 
	The second phase of error correction goes through clumpify.sh
	*/

	tag "${sample_id}"
    label "general"
	// publishDir params.error_correct, pattern: "*.fastq.gz", mode: params.publishMode, overwrite: true

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_error_correct2.fastq.gz")

	script:
	"""
	clumpify.sh in=`realpath ${reads}` \
	out=${sample_id}_error_correct2.fastq.gz \
	ecc passes=4 reorder \
	threads=${task.cpus}
	"""

}

process ERROR_CORRECT_PHASE_THREE {

	/* 
	The third phase of error correction uses tadpole.sh.
	*/

	tag "${sample_id}"
    label "general"
	// publishDir params.error_correct, pattern: "*.fastq.gz", mode: params.publishMode, overwrite: true

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_error_correct3.fastq.gz")

	script:
	"""
	tadpole.sh in=`realpath ${reads}` \
	out=${sample_id}_error_correct3.fastq.gz \
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
	// publishDir params.qtrim, pattern: "*.fastq.gz", mode: 'copy', overwrite: true

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_qtrimmed.fastq.gz")

	script:
	"""
	bbduk.sh in=`realpath ${reads}` \
	out=${sample_id}_qtrimmed.fastq.gz \
	qtrim=rl trimq=10 minlen=70 ordered \
	threads=${task.cpus}
	"""

}

process EXTRACT_REF_AMPLICON {

    /* */

    input:
    path primer_file

    output:
    path "amplicon.fasta", emit: seq
    env len, emit: length

    script:
    """
    seqkit amplicon \
    --max-mismatch 0 \
    --primer-file ${primer_file}
    ${params.reference} \
    > amplicon.fasta && \
    seqkit fx2tab --no-qual --length amplicon.fasta -o amplicon.stats && \
    len=`cat amplicon.stats | tail -n 1 | awk '{print $2}'`
    """

}

process MAP_TO_AMPLICON {

    /*
    */
	
	tag "${sample_id}"
    label "general"
	publishDir params.results, mode: 'copy'

    cpus params.available_cpus
	
	input:
	tuple val(sample_id), path(reads)
    path amplicon_seq
	
	output:
	tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai")
	
	script:
    if ( params.geneious_mode == true )
        """
        geneious -i ${amplicon_seq} ${reads} -operation Map_to_Reference -o ${sample_id}.bam && \
        samtools sort -o ${sample_id}_sorted.bam ${sample_id}.bam && \
        samtools index -o ${sample_id}_sorted.bam.bai ${sample_id}_sorted.bam
        """
    else
        """
        bbmap.sh ref=${amplicon_seq} in=${reads} out=stdout.sam t=${task.cpus} maxindel=200 | \
        samtools sort -o ${sample_id}_sorted.bam - && \
        samtools index -o ${sample_id}_sorted.bam.bai ${sample_id}_sorted.bam
        """
}

process CLIP_AMPLICONS {

    /*
    */
	
	tag "${sample_id}"
    label "general"
	publishDir params.results, mode: 'copy'

    cpus 3
	
	input:
	tuple val(sample_id), path(bam), path(index)
    path amplicon_bed
	
	output:
    tuple val(sample_id), path("${sample_id}_clipped.bam")
	
	script:
	"""
	samtools ampliconclip -b ${amplicon_bed} ${sample_id}_sorted.bam -o ${sample_id}_clipped.bam
	"""
}

process BAM_TO_FASTQ {

    /*
    */
	
	tag "${sample_id}"
    label "general"
	publishDir params.results, mode: 'copy'

    cpus params.available_cpus
	
	input:
	tuple val(sample_id), path(bam)
	
	output:
	tuple val(sample_id), path("*.fastq.gz")
	
	script:
	"""
	reformat.sh in=${bam} out=stdout.fastq.gz t=${task.cpus} | \
	clumpify.sh in=stdin.fastq.gz out=${sample_id}_amplicon_reads.fastq.gz t=${task.cpus} reorder
	"""
}

process FILTER_BY_LENGTH {

    /*
    */
	
	tag "${sample_id}"
    label "general"
	publishDir params.results, mode: 'copy'

    cpus 3
	
	input:
    val max_len
	tuple val(sample_id), path(reads)
	
	output:
	tuple val(sample_id), path("*.fastq.gz")
	
	script:
	"""
    seqkit seqkit \
    --max-len ${max_len} \
    ${reads} \
    -o ${sample_id}_filtered.fastq.gz
	"""
}

process HAPLOTYPE_ASSEMBLY {

    /*
    */
	
	tag "${sample_id}"
    label "general"
	publishDir params.results, mode: 'copy'

    cpus 8
	
	input:
	tuple val(sample_id), path(reads)
	
	output:
	tuple val(sample_id), path(bam)
	
	script:
    if ( params.geneious_mode == true )
        """
        geneious -i ${reads} –operation de_novo_assemble -x ${params.assembly_profile} -o ${sample_id}.bam
        """
    else
        """
        spades.py --merged ${reads} --corona --threads ${task.cpus} -o .
        """
}

// process RECORD_FREQUENCIES {

//     /*
//     */
	
// 	tag "${sample_id}"
//    label "general"
// 	publishDir params.results, mode: 'copy'
	
// 	input:
	
	
// 	output:
	
	
// 	script:
// 	"""
	
// 	"""
// }

process DOWNSAMPLE_ASSEMBLIES {

    /*
    */
	
	tag "${sample_id}"
    label "general"
	publishDir params.results, mode: 'copy'
	
	input:
	tuple val(sample_id), path(bam)
	
	output:
	tuple val(sample_id), path("*.bam")
	
	script:
	"""
	samtools view -b -s 0.1 ${bam} > ${sample_id}_downsampled.bam
	"""
}

process MAP_CONTIGS_TO_REF {

    /*
    */
	
	tag "${sample_id}"
    label "general"
	publishDir params.results, mode: 'copy'

    cpus 3
	
	input:
	tuple val(sample_id), path(bam)
	
	output:
	tuple val(sample_id), path("*.bam")
	
	script:
	if ( params.geneious_mode == true )
        """
        echo in development > false.bam
        """
    else
        """
        reformat.sh in=${bam} out=stdout.fastq.gz | \
        bbmap.sh ref=${params.reference} in=${reads} out=stdout.sam | \
        reformat.sh in=stdin.sam out=${sample_id}.bam 
        """
}

process TRIM_PRIMERS {

    /*
    */
	
	tag "${sample_id}"
    label "iVar"
	publishDir params.results, mode: 'copy'

    cpus 3
	
	input:
	tuple val(sample_id), path(bam)
	
	output:
	tuple val(sample_id), path("*.bam")
	
	script:
	"""
    ivar trim -b ${params.primer_bed} -p sample_id_trimmed -i ${bam} -q 15 -m 50 -s 4
	"""
}

process CALL_CONSENSUS_SEQS {

    /*
    */
	
	tag "${sample_id}"
    label "general"
	publishDir params.results, mode: 'copy'

    cpus 3
	
	input:
	tuple val(sample_id), path(bam)
	
	output:
	tuple val(sample_id), path("*.fasta")
	
	script:
	"""
    pileup.sh in=${bam} out=${sample_id}_pileup.txt consensus=${sample_id}_consensus.fasta
	"""
}

process CALL_VARIANTS {

    /*
    */
	
	tag "${sample_id}"
    label "general"
	publishDir params.results, mode: 'copy'

    cpus 3
	
	input:
	tuple val(sample_id), path(bam)
	
	output:
	tuple val(sample_id), path("*.vcf")
	
	script:
	"""
    bcftools mpileup -Ou -f ${params.reference} ${bam} | bcftools call -mv -Ov -o ${sample_id}.vcf
	"""

}

process RUN_IVAR {

    /*
    */
	
	tag "${sample_id}"
    label "iVar"
	publishDir params.results, mode: 'copy'

    cpus 3
	
	input:
	tuple val(sample_id), path(bam)
	
	output:
	tuple val(sample_id), path("*.tsv"), path("*.fasta")
	
	script:
	"""
    samtools mpileup -aa -A -d 600000 -B -Q 0 ${bam} | \
    tee >(ivar variants -p ${sample_id} -q 20 -t 0.03 -r ${params.reference} -g ${params.gff}) | \
    ivar consensus -p ${sample_id}_ivar_consensus -q 20 -t 0
	"""

}

// process GENERATE_REPORT {
	
// 	// This process does something described here
	
// 	tag "${sample_id}"
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
