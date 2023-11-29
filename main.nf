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
	
	// Workflow steps
    GET_PRIMER_SEQS (
        ch_primer_bed
    )

    EXTRACT_REF_AMPLICON (
		ch_refseq,
        GET_PRIMER_SEQS.out.amplicon_coords
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
        GET_PRIMER_SEQS.out.patterns
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

    MAP_TO_AMPLICON (
        QUALITY_TRIM.out
			.map { name, reads -> tuple( name, file(reads), file(reads).countFastq() ) }
			.filter { it[2] >= params.min_reads }
			.map { name, reads, count -> tuple( name, file(reads) ) },
        EXTRACT_REF_AMPLICON.out.seq
    )

    CLIP_AMPLICONS (
        MAP_TO_AMPLICON.out,
        GET_PRIMER_SEQS.out.bed
    )

	if ( params.geneious_mode != true ) {

			CALL_AMPLICON_VARIANTS (
				CLIP_AMPLICONS.out
			)

	}

    BAM_TO_FASTQ (
        CLIP_AMPLICONS.out
    )

    VALIDATE_SEQS (
        BAM_TO_FASTQ.out
    )

	if ( params.geneious_mode == true ) {

		ASSEMBLE_WITH_GENEIOUS (
			VALIDATE_SEQS.out
				.map { name, reads -> tuple( name, file(reads), file(reads).countFastq() ) }
				.filter { it[2] >= params.min_reads }
				.map { name, reads, count -> tuple( name, file(reads) ) }
		)

		FILTER_ASSEMBLIES (
			ASSEMBLE_WITH_GENEIOUS.out.flatten()
		)

	} else {

		ASSEMBLE_WITH_WHATSHAP (
			VALIDATE_SEQS.out
				.map { name, reads -> tuple( name, file(reads), file(reads).countFastq() ) }
				.filter { it[2] >= params.min_reads }
				.map { name, reads, count -> tuple( name, file(reads) ) },
			CALL_AMPLICON_VARIANTS.out.vcf
		)

		if ( params.long_reads == true ) {

			RUN_HAIRSPLITTER (
				ASSEMBLE_WITH_WHATSHAP.out.flatten()
			)

			FILTER_ASSEMBLIES (
				RUN_HAIRSPLITTER.out
			)

		} else {

			FILTER_ASSEMBLIES (
				ASSEMBLE_WITH_WHATSHAP.out.flatten()
			)

		}

	}

    MAP_ASSEMBLY_TO_REF (
        FILTER_ASSEMBLIES.out
			.map { filename, reads, id -> tuple( filename, file(reads), id, file(reads).countFastq() ) }
			.filter { it[3] >= params.min_reads },
		ch_refseq
    )

	CALL_CONSENSUS_SEQS (
		MAP_ASSEMBLY_TO_REF.out
	)

	EXTRACT_AMPLICON_CONSENSUS (
		CALL_CONSENSUS_SEQS.out
	)

    CALL_HAPLOTYPE_VARIANTS (
        MAP_ASSEMBLY_TO_REF.out,
		ch_refseq
    )

	GENERATE_IVAR_TABLE (
        MAP_ASSEMBLY_TO_REF.out,
		ch_refseq,
		ch_refgff
	)

    // GENERATE_REPORT (
    //     DOWNSAMPLE_ASSEMBLIES.out,
    //     MAP_ASSEMBLY_TO_REF.out,
    //     CALL_CONSENSUS_SEQS.out,
    //     CALL_HAPLOTYPE_VARIANTS.out
    // )
	
	
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config

// whether to terminate when unrecoverable error occurs
if ( params.debugmode ) {
	params.errorMode = 'terminate'
} else {
	params.errorMode = 'ignore'
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

// assembly results
params.assembly_results = params.amplicon_results + "/02_assembly_results"
params.assembly_reads = params.assembly_results + "/01_assembly_reads"
params.aligned_assembly = params.assembly_results + "/02_aligned_assemblies"
params.consensus = params.assembly_results + "/03_contigs"
params.variants = params.assembly_results + "/04_contig_VCFs"
params.ivar_tables = params.assembly_results + "/05_ivar_tables"


// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process GET_PRIMER_SEQS {

    /*
    */

	tag "${params.desired_amplicon}"
    label "general"

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

    input:
    path bed_file

    output:
    path "primer_seqs.bed", emit: bed
    path "primer_seqs.txt", emit: txt
	path "amplicon_coords.bed", emit: amplicon_coords
	path "patterns.txt", emit: patterns

    script:
    """

	# get the actual primer sequences and make sure the amplicon's reverse primer
	# is complementary to the reference sequence it's based on
    grep ${params.desired_amplicon} ${params.primer_bed} > primer_seqs.bed && \
    bedtools getfasta -fi ${params.reference} -bed primer_seqs.bed > tmp.fasta && \
	grep -v "^>" tmp.fasta > patterns.txt && \
	seqkit head -n 1 tmp.fasta -o primer_seqs.fasta && \
	seqkit range -r -1:-1 tmp.fasta | \
	seqkit seq --complement --validate-seq --seq-type DNA >> primer_seqs.fasta && \
	rm tmp.fasta

	# determine the amplicon coordinates
	cat primer_seqs.bed | \
	awk 'NR==1{start=\$2} NR==2{end=\$3} END{print \$1, start, end}' OFS="\t" \
	> amplicon_coords.bed

	# convert to a text file that can be read by `seqkit amplicon`
    seqkit fx2tab --no-qual primer_seqs.fasta | \
	awk '{print \$2}' | paste -sd \$'\t' - - | awk -v amplicon="${params.desired_amplicon}" \
	'BEGIN {OFS="\t"} {print amplicon, \$0}' | head -n 1 > primer_seqs.txt
    """

}

process EXTRACT_REF_AMPLICON {

    /* */

	tag "${params.desired_amplicon}"
    label "general"
	publishDir params.amplicon_results, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

    input:
	path refseq
    path amplicon_coords

    output:
    path "${params.desired_amplicon}.fasta", emit: seq
    env len, emit: length

    script:
    """
	cat ${refseq} | \
    seqkit subseq \
	--bed ${amplicon_coords} \
    -o ${params.desired_amplicon}.fasta && \
    seqkit fx2tab --no-qual --length ${params.desired_amplicon}.fasta \
	-o ${params.desired_amplicon}.stats && \
    len=`cat ${params.desired_amplicon}.stats | tail -n 1 | awk '{print \$3}'`
    """

}

process MERGE_PAIRS {

    /* */
	
	tag "${sample_id}"
    label "general"
	publishDir params.merged_reads, mode: 'copy', overwrite: true

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
	publishDir params.clumped_reads, mode: 'copy', overwrite: true

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

	cpus 4
	
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
	publishDir params.trim_adapters, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 4
	
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
    tuple val(sample_id), path("${sample_id}_amplicons.fastq.gz")

    script:
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
    -o ${sample_id}_amplicons.fastq.gz
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

	cpus 4
	memory 3.GB

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_deduped.fastq.gz")

	script:
	if ( sample_id.toString().contains("SRR") )
		"""
		rename.sh \
		in=${reads} \
		out=${sample_id}_deduped.fastq.gz \
		addpairnum=t
		"""
	else
		"""
		clumpify.sh -Xmx2g in=`realpath ${reads}` \
		out=${sample_id}_deduped.fastq.gz \
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

	cpus 4
	memory 3.GB

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_filtered.fastq.gz")

	script:
	if ( sample_id.toString().contains("SRR") )
		"""
		clumpify.sh -Xmx2g in=`realpath ${reads}` \
		out=${sample_id}_filtered.fastq.gz \
		threads=${task.cpus} \
		reorder markduplicates
		"""
	else
		"""
		filterbytile.sh -Xmx2g in=`realpath ${reads}` \
		out=${sample_id}_filtered.fastq.gz \
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

	cpus 4
	memory 3.GB

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_remove_artifacts.fastq.gz")

	script:
	"""
	bbduk.sh -Xmx2g in=`realpath ${reads}` \
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

	tag "${sample_id}"
    label "general"
	publishDir params.error_correct, pattern: "*.fastq.gz", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 4
	memory 3.GB

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_error_correct1.fastq.gz")

	script:
	"""
	bbmerge.sh -Xmx2g in=`realpath ${reads}` \
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
	publishDir params.error_correct, pattern: "*.fastq.gz", mode: 'copy', overwrite: true

	cpus 4
	memory 3.GB

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_error_correct2.fastq.gz")

	script:
	"""
	clumpify.sh -Xmx2g in=`realpath ${reads}` \
	out=${sample_id}_error_correct2.fastq.gz \
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

	cpus 4
	memory 3.GB

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_error_correct3.fastq.gz")

	script:
	"""
	tadpole.sh -Xmx2g in=`realpath ${reads}` \
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
	publishDir params.qtrim, pattern: "*.fastq.gz", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 4
	memory 3.GB

	input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_qtrimmed.fastq.gz")

	script:
	"""
	bbduk.sh -Xmx2g in=`realpath ${reads}` \
	out=${sample_id}_qtrimmed.fastq.gz \
	qtrim=rl trimq=10 minlen=70 ordered \
	threads=${task.cpus}
	"""

}

process MAP_TO_AMPLICON {

    /*
    */
	
	tag "${sample_id}"
    label "general"

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 4
	
	input:
	tuple val(sample_id), path(reads)
    each path(amplicon_seq)
	
	output:
	tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai")
	
	script:
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
	publishDir params.clipped, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 4
	
	input:
	tuple val(sample_id), path(bam), path(index)
    each path(amplicon_bed)
	
	output:
    tuple val(sample_id), path("${sample_id}_clipped.bam")
	
	script:
	"""
	samtools ampliconclip \
	-b ${amplicon_bed} \
	--hard-clip \
	--both-ends \
	--clipped \
	${sample_id}_sorted.bam \
	-o ${sample_id}_clipped.bam
	"""
}

// process CALL_AMPLICON_VARIANTS {}

process BAM_TO_FASTQ {

    /*
    */
	
	tag "${sample_id}"
    label "general"
	publishDir params.clipped, mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 4
	
	input:
	tuple val(sample_id), path(bam)
	
	output:
	tuple val(sample_id), path("*.fastq.gz")
	
	script:
	"""
	reformat.sh in=${bam} out=stdout.fastq t=${task.cpus} | \
	clumpify.sh in=stdin.fastq out=${sample_id}_amplicon_reads.fastq.gz t=${task.cpus} reorder
	"""
}

process VALIDATE_SEQS {

    /*
    */
	
	tag "${sample_id}"
    label "general"
	publishDir params.complete_amplicon, mode: 'copy'

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

	cpus 4
	
	input:
	tuple val(sample_id), path(reads)
	
	output:
	tuple val(sample_id), path("*.fastq.gz")
	
	script:
	"""
	cat ${reads} | \
    seqkit seq \
	--threads ${task.cpus} \
	--remove-gaps \
	--validate-seq \
    -o ${sample_id}_filtered.fastq.gz
	"""
}

process ASSEMBLE_WITH_GENEIOUS {

    /*
    */
	
	tag "${sample_id}"
    // label "general"

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

    cpus 8
	
	input:
	tuple val(sample_id), path(reads)
	
	output:
	path "*.fastq.gz"

	when:
	params.geneious_mode == true
	
	script:
	"""
	geneious -i ${reads} -x ${params.assembly_profile} -o ${sample_id}.fastq.gz --multi-file && \
	for file in ${sample_id}*Assembly*.fastq.gz; do
		mv "\$file" "\${file// /_}"
	done
	"""
}

// process ASSEMBLE_WITH_WHATSHAP {}

// process RUN_HAIRSPLITTER {}

process FILTER_ASSEMBLIES {

	/* */

	tag "${file_name}"
    label "general"
	publishDir "${params.assembly_reads}/${sample_id}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

    cpus 4

	input:
	path assembly_reads

	output:
	tuple val(file_name), path("*.fastq.gz"), val(sample_id)

	when:
	assembly_reads.getSimpleName().contains("Contig")

	script:
	sample_id = assembly_reads.getSimpleName().split("_L00")[0]
	contig_num = assembly_reads.getName().split(".fastq.gz")[0].replace("Unpaired", "").split("_")[-1]
	file_name = sample_id + "_contig" + contig_num
	"""
	clumpify.sh in="${assembly_reads}" out="${file_name}.fastq.gz" reorder
	"""

}

process MAP_ASSEMBLY_TO_REF {

    /*
    */
	
	tag "${name}"
    label "general"
	publishDir "${params.aligned_assembly}/${sample_id}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

    cpus 4
	
	input:
	tuple val(name), path(assembly_reads), val(sample_id), val(read_support)
	each path(refseq)
	
	output:
	tuple val(name), path("*.bam"), val(sample_id), val(read_support)
	
	script:
	"""
	bbmap.sh int=f ref=${refseq} in=${assembly_reads} out=stdout.sam maxindel=200 | \
	reformat.sh in=stdin.sam out="${name}.bam"
	"""

}

process CALL_CONSENSUS_SEQS {

    /*
    */
	
	tag "${name}"
    label "iVar"

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2

    cpus 4
	
	input:
	tuple val(name), path(bam), val(sample_id), val(read_support)
	
	output:
	tuple val(name), path("${name}_consensus.fa*"), val(sample_id), val(read_support)
	
	script:
	"""
	samtools sort "${bam}" | \
	samtools mpileup \
	-aa -A -d 0 -Q 0 \
	- | ivar consensus \
	-t 0 -m 1 -q 0 -k \
	-p ${name}_consensus
	"""

}

process EXTRACT_AMPLICON_CONSENSUS {

	/* */
	
	tag "${name}"
    label "general"
	publishDir "${params.consensus}/${sample_id}", mode: 'copy', overwrite: true

	errorStrategy { task.attempt < 3 ? 'retry' : params.errorMode }
	maxRetries 2
	
	input:
	tuple val(name), path(fasta), val(sample_id), val(read_support)
	
	output:
	tuple val(name), path("${name}_${params.desired_amplicon}_consensus.fa*"), val(sample_id), val(read_support)
	
	script:
	"""
	cat ${fasta} | \
	seqkit replace \
	-p \$ -r "_readsupport_${read_support}" \
	-o "${name}_${params.desired_amplicon}_consensus.fa"
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

    cpus 4
	
	input:
	tuple val(name), path(bam), val(sample_id), val(read_support)
	each path(refseq)
	
	output:
	tuple val(name), path("${name}.vcf"), val(sample_id), val(read_support)
	
	script:
	"""
    bcftools mpileup -Ou -f ${refseq} ${bam} | bcftools call --ploidy 1 -mv -Ov -o ${name}.vcf
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
	tuple val(name), path(bam), val(sample_id), val(read_support)
	each path(refseq)
	each path(refgff)

	output:
	path "*.tsv"

	script:
	"""
	samtools mpileup -aa -A -d 0 -B -Q 0 --reference ${refseq} ${bam} | \
	ivar variants -p ${name} -t 0 -r ${refseq} -g ${refgff}
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
