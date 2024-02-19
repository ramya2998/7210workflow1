// Defining the input parameters.
params.fastQ_reads = "$projectDir/data/training_data/ANN0831_R{1,2}.fastq.gz"
params.output_directory = "Outputs"

// Log pipeline information (GitHub)
log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    reads: ${params.fastQ_reads}
    outdir: ${params.output_directory}
"""
    .stripIndent(true)

// Defining the fastQ trimming process
process fastQ_trimming {
    // Identification tag in logs
    tag "Data Being Trimmed: $sample_id"
    // Output directory --> Results
    publishDir params.output_directory, mode: 'copy'

    // Input tuple of sample ID and fastq reads
    input:
    tuple val(sample_id), path(fastQ_reads)

    // Output tuple of sample IDs and paths to paired/unpaired trimmed fastq files
    output:
    tuple val(sample_id), path("${sample_id}_r1.paired.fastq.gz"), path("${sample_id}_r2.paired.fastq.gz")

    // Actual trimming...
    script:

    """
    trimmomatic PE -phred33 \\
        ${fastQ_reads[0]} \\
        ${fastQ_reads[1]} \\
        ${sample_id}_r1.paired.fastq.gz \\
        ${sample_id}_r1.unpaired.fastq.gz \\
        ${sample_id}_r2.paired.fastq.gz \\
        ${sample_id}_r2.unpaired.fastq.gz \\
        SLIDINGWINDOW:5:30 \\
        1> trimmed_stdout.log \\
        2> trimmed_stderr.log 
    """
}

// Defining fasta assembly process
process fasta_assembly {
    // Identification tags, again. 
    tag "Sample Being Assembled: $sample_id"
    // Specified output directory...
    publishDir params.output_directory, mode: 'copy'

    // Input tuple of sample ID and paths to paired fastq files.
    input:
    tuple val(sample_id), path(r1_paired), path(r2_paired)

    // Output path for assembled fasta file
    output:
    path("${sample_id}_assembled")

    // Actual assembly process
    """
    spades.py -1 $r1_paired -2 $r2_paired -o ${sample_id}_assembled
    """
}

// Defining the workflow...
workflow {
    // Channel to read fastq files in paris
    fastQ_pairs = Channel.fromFilePairs(params.fastQ_reads, checkIfExists: true)
    // Perform fastq trimming on the pairs
    trimmed_files = fastQ_trimming(fastQ_pairs)
    // Perform fasta assembly on trimmed files 
    assembled_files = fasta_assembly(trimmed_files)
}