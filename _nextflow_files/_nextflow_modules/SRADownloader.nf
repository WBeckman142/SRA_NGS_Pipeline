#!/usr/bin/env nextflow

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Download 100,000 reads from target files using fastq-dump. Split files into read pairs
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

 process sradownloader {

   tag { srr }

    publishDir 'data/raw', mode: 'copy'
    container 'community.wave.seqera.io/library/sra-tools:3.2.1--2063130dadd340c5'

    input:
        val srr

    output:
       tuple val(srr), path("*_1.fastq"), path("*_2.fastq")

    script:
    
    // If `params.max_reads` is set, add `-X N`, otherwise leave it empty
    def max_reads_opt = params.max_reads ? "-X ${params.max_reads}" : ""

    """
    fastq-dump --split-files ${max_reads_opt} --outdir ./ ${srr}
    """
   }
