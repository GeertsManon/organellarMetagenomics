# organellarMetagenomics

A bioinformatics pipeline designed for the recovery of complete organellar genomes from metagenomic data. This workflow efficiently extracts high-quality organellar genomes from environmental DNA (short-read sequencing), uncovering previously hidden eukaryotic diversity.

For detailed instructions on how to use this pipeline, please refer to my Wiki pages:

- [Installation and Dependencies](https://github.com/GeertsManon/organellarMetagenomics/wiki/A:-Installation-&-Dependencies) - Overview of the system and tool requirements
- [Upstream Pipeline Workflow using HPC](https://github.com/GeertsManon/organellarMetagenomics/wiki/B:-Upstream-pipeline-workflow-using-HPC) - Instructions for running the upstream pipeline (data retrieval, read trimming and assembly) on high-performance computing environments
- [Manual Disentanglement of Chloroplast Genomes from Metagenomic Data](https://github.com/GeertsManon/organellarMetagenomics/wiki/C:-Manual-Disentanglement-of-Chloroplast-Genomes-from-Metagenomic-Data-in-Bandage) - Step-by-step guide for manually reconstructing complete chloroplast genomes from assembly graphs
- [Quality Assessment of Assembled Organellar Genomes](https://github.com/GeertsManon/organellarMetagenomics/wiki/D:-Quality-Assessment-of-Assembled-Organellar-Genomes) - Description of the quality assessment process (mapping and variant calling)

[Scripts](https://github.com/GeertsManon/organellarMetagenomics/tree/main/scripts) are developed for:

- Data downloading from NCBI SRA using SRA-Toolkit
- Quality trimming with Trimmomatic
- _De novo_ assembly with [GetOrganelle](https://github.com/Kinggerm/GetOrganelle)
- Quality assessment through 1) mapping with SMALT and SAMtools and 2) variant calling with BCFtools
