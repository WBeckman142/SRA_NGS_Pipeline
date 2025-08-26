#!/usr/bin/env nextflow


/**
* download 50,000 reads from target files using fastq-dump. Split files into read pairs
*/


 process sradownloader {

    publishDir 'data/raw', mode: 'copy'
    container 'community.wave.seqera.io/library/sra-tools:3.2.1--2063130dadd340c5'

    input:
        val srr

    output:
       file("*.fastq")

    script:
    """
    fastq-dump  --split-files -X 50000 --outdir ./ ${srr}  
    """
   }
