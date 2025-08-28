# 🧬 NGS Processing Pipeline

![Nextflow](https://img.shields.io/badge/nextflow-%2300A388.svg?style=flat&logo=nextflow&logoColor=white)
![Docker](https://img.shields.io/badge/docker-%230db7ed.svg?style=flat&logo=docker&logoColor=white)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![GitHub last commit](https://img.shields.io/github/last-commit/WBeckman142/SRA_NGS_Pipeline/main)


This repository contains a **Nextflow pipeline** for processing NGS data from SRA.  
It downloads FASTQ files, runs **FastQC** for quality control, and aligns reads to the **hg19 human genome** using **Bowtie2**, producing sorted and indexed BAM files ready for visualization in **IGV**.

---

## 📖 Overview
- Automates NGS data preprocessing with reproducibility in mind.
- Designed to work with SRA accession codes.
- Allows you to download sepcified number of reads per sample.
- Generates alignment-ready BAM files with indexes for downstream analysis.

---

## 🚀 Features
- 🔽 Download FASTQ files from SRA
- 📊 Run FastQC for quality assessment
- 🎯 Align paired-end reads with Bowtie2
- 📂 Output sorted & indexed BAM files
- 👀 Ready-to-visualize in IGV

---

## 📂 Repository Structure
```
├── data/ # Raw and processed data
│ ├── raw/ # Input FASTQ files
│ └── processed/ # Aligned BAM outputs
├── nextflow_files/ # Nextflow modules
│ ├── FASTQC_module.nf
│ ├── SRADownloader.nf
│ └── bowtie2_module.nf
├── reference_genome/ # Bowtie2 index files for hg19
├── txt_inputs/ # SRA accession list (srr_codes.txt)
└── main.nf # Main Nextflow pipeline
```

---

## ⚡ Installation
Make sure you have:
- [Nextflow](https://www.nextflow.io/)  
- [Docker](https://www.docker.com/) or [Apptainer](https://apptainer.org/)  

Clone this repository:
```bash
git clone https://github.com/WBeckman142/SRA_NGS_Pipeline.git
cd SRA_NGS_Pipeline
```


▶️ Usage

Add your SRA accession codes to txt_inputs/srr_codes.txt.<br>
Add your indexed hg19 reference genome to the folder reference_genome in the main directory<br>
<br>
Run the pipeline:
```
nextflow run main.nf
```
By default, the pipeline downloads all of the reads from each SRA accession. You can override this with the --max_reads parameter:
```
nextflow run main.nf --max_reads 50000
```
- --max_reads: Maximum number of reads to download from each SRA accession.

📊 Example Output

After running the pipeline, you’ll find:<br>
- FastQC HTML reports in reports/<br>
- Aligned BAM files + indexes in data/processed/

Example:
```
data/processed/
├── SRR30621347.sorted.bam
├── SRR30621347.sorted.bam.bai
```
In IGV:
![IGV Tracks](https://github.com/WBeckman142/SRA_NGS_Pipeline/blob/main/images/IGV_tracks.png?raw=true)

🤝 Contributing

Contributions welcome! 
