#!/usr/bin/env nextflow


/**
* run bowtie2 to align the reads to the human genome build hg19
*/


process run_bowtie2_aligner {
    maxForks 1

    

    tag { srr }

    publishDir "data/processed", mode: 'copy'
    container "community.wave.seqera.io/library/bowtie2_samtools:a08292672802bd5b"
    

    input:
        tuple val(srr), path(read1), path(read2) // reads is a list: [R1, R2]
        path index_base

    output:
        tuple val(srr), path("*.bam")

    script:

    """
    bowtie2 -x ${index_base}/hg19 -1 ${read1} -2 ${read2} | samtools view -bS - > "${srr}.bam"
    """
}



