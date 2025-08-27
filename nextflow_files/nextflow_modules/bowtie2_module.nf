#!/usr/bin/env nextflow


/**
* run bowtie2 to align the reads to the human genome build hg19
*/


process run_bowtie2_aligner {

    tag "${sample_id}"

    publishDir "data/processed", mode: 'copy'
    container "community.wave.seqera.io/library/bowtie2_samtools:a08292672802bd5b"
    

    input:
    tuple val(sample_id), path(reads) // reads is a list: [R1, R2]
    val index_base

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    bowtie2 -x ${index_base} \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        | samtools view -bS - > ${sample_id}.bam
    """
}



