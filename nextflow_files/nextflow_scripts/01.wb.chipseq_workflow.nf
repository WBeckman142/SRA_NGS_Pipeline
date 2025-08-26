#!/usr/bin/env nextflow

/**
* Import srr_codes.txt containing SRR codes of files for download
*
*/

input_sra_codes = 'txt_inputs/srr_codes.txt'

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

    """
    fastq-dump  --split-files -X 50000 --outdir ./ ${srr}  
    """
   }


workflow{

    // download fastq files from SRA using fastq-dump
    fastq_ch = sradownloader(srr_ch)
    
}
