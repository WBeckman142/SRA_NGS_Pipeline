#!/usr/bin/env nextflow

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Import srr_codes.txt containing SRR codes of files for download
* Import raw fastq files
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

input_sra_codes = 'txt_inputs/srr_codes.txt'
params.reads    = "${launchDir}/data/raw/*.fastq"



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Create channels
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

srr_ch          = Channel.fromPath(input_sra_codes)
                         .splitText()

index_ch        = Channel.value("/Users/willbeckman/Documents/Nextflow/reference_genome/hg19")


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Import modules from nextflow_files/nextflow_modules folder
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

include { sradownloader          } from "${launchDir}/nextflow_files/nextflow_modules/SRADownloader.nf"
include { run_fastqc             } from "${launchDir}/nextflow_files/nextflow_modules/FASTQC_module.nf"
include { run_bowtie2_aligner    } from "${launchDir}/nextflow_files/nextflow_modules/bowtie2_module.nf"




/**
* run workflow using imported modules
*/
workflow{

    // download fastq files from SRA using fastq-dump
    fastq_ch = sradownloader(srr_ch)

    // run fastqc
    run_fastqc(fastq_ch)

    read_pairs_ch   = Channel.fromFilePairs("data/raw/*_{1,2}.fastq")                    
                             .map { id, reads -> tuple(id, reads) }
                             .view()

    // run bowtie2
    run_bowtie2_aligner(read_pairs_ch, index_ch)
}
