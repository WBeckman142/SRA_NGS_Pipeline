#!/usr/bin/env nextflow

/**
* Import srr_codes.txt containing SRR codes of files for download
*
*/

input_sra_codes = 'txt_inputs/srr_codes.txt'
params.reads = "${launchDir}/data/raw/*.fastq"

/**
* Create channel from imported file and split to text
*/

srr_ch = Channel.fromPath(input_sra_codes)
                .splitText()
                .view()

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


workflow{

    // download fastq files from SRA using fastq-dump
    fastq_ch = sradownloader(srr_ch)

    run_fastqc(sradownloader.out)
}
