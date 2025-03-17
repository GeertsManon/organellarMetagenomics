# OrganellarMetagenomics

A bioinformatics pipeline for recovering complete organellar genomes from metagenomic data. Enables extraction of high-quality organellar genomes from environmental DNA (short reads), revealing hidden eukaryotic diversity.

[Scripts](https://github.com/GeertsManon/organellarMetagenomics/tree/main/scripts) for:

- Data download from NCBI SRA using SRA-Toolkit
- Quality trimming with Trimmomatic
- _De novo_ assembly with [GetOrganelle](https://github.com/Kinggerm/GetOrganelle)
- Quality assessment through 1) mapping with SMALT and SAMtools and 2) variant calling with BCFtools

For detailed instructions on how to use this pipeline, please refer to my Wiki pages:

- [Upstream Pipeline Workflow using HPC](https://github.com/GeertsManon/organellarMetagenomics/wiki/Upstream-pipeline-workflow-using-HPC) - Instructions for running the pipeline on high-performance computing environments
- [Manual Disentanglement of Chloroplast Genomes from Metagenomic Data](https://github.com/GeertsManon/organellarMetagenomics/wiki/Manual-Disentanglement-of-Chloroplast-Genomes-from-Metagenomic-Data) - Step-by-step guide for manually reconstructing complete chloroplast genomes from assembly graphs
