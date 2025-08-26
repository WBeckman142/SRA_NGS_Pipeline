#!/usr/bin/env nextflow


/**
* run fastqc from docker image and save the output html to reports folder
*/


process run_fastqc {

    publishDir 'reports', mode: 'copy'
    container 'community.wave.seqera.io/library/fastqc:0.12.1--af7a5314d5015c29'

    input:
        path read  
    
    output:
        file("*.html")

    script:
    """
    fastqc ${read}
    """
}

