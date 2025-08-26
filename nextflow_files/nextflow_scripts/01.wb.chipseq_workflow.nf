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


/**
* Import modules from nextflow_modules folder
*/
include { sradownloader } from "${launchDir}/nextflow_files/nextflow_modules/SRADownloader.nf"
include { run_fastqc } from "${launchDir}/nextflow_files/nextflow_modules/FASTQC_module.nf"


/**
* run workflow using imported modules
*/
workflow{

    // download fastq files from SRA using fastq-dump
    fastq_ch = sradownloader(srr_ch)

    // run fastqc
    run_fastqc(fastq_ch)
}
