# OrganellarMetagenomics

A bioinformatics pipeline for recovering complete organellar genomes from metagenomic data. Scripts for:

- Data download from NCBI SRA using SRA-Toolkit
- Quality trimming with Trimmomatic
- De novo assembly with GetOrganelle
- Quality assessment through 1) mapping with SMALT and SAMtools and 2) variant calling with BCFtools

Enables extraction of high-quality organellar genomes from environmental DNA (short reads), revealing hidden eukaryotic diversity.
