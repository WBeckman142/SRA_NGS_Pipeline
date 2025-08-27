#!/usr/bin/env nextflow


/**
* run fastqc from docker image and save the output html to reports folder
*/


process run_fastqc {

    publishDir 'reports', mode: 'copy'
    container 'community.wave.seqera.io/library/fastqc:0.12.1--af7a5314d5015c29'

    input:
        tuple val(srr), path(read1), path(read2)  
    
    output:
        file("*.html")

    script:
    """
    fastqc ${read1} ${read2}
    """
}

