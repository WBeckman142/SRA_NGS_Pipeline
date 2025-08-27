#!/usr/bin/env nextflow


/**
* download 50,000 reads from target files using fastq-dump. Split files into read pairs
*/


 process sradownloader {

   tag { srr }

    publishDir 'data/raw', mode: 'copy'
    container 'community.wave.seqera.io/library/sra-tools:3.2.1--2063130dadd340c5'

    input:
        val srr

    output:
       tuple val(srr), path("*_1.fastq"), path("*_2.fastq")

    script:
    """
    fastq-dump  --split-files -X 5 --outdir ./ ${srr}  
    """
   }
