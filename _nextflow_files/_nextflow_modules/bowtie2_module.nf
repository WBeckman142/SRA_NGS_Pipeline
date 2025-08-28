#!/usr/bin/env nextflow

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Run Bowtie2 and align to the human genome hg19
* Run samtools to convert sam files to bam files then sort and index
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

process run_bowtie2_aligner {
    // limited parallelisation to run on local machine
    maxForks 1

    tag { srr }

    publishDir "data/processed", mode: 'copy'
    container "community.wave.seqera.io/library/bowtie2_samtools:a08292672802bd5b"
    

    input:
        tuple val(srr), path(read1), path(read2)
        path index_base

    output:
        tuple val(srr), path("*.sorted.bam"), path("*.sorted.bam.bai")

    script:

    """
    bowtie2 -x ${index_base}/hg19 -1 ${read1} -2 ${read2} \
      | samtools view -bS - \
      | samtools sort -o ${srr}.sorted.bam
    
    samtools index -b ${srr}.sorted.bam
    """
}



