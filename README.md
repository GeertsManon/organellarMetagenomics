# organellarMetagenomics

A bioinformatic toolkit for extracting and analyzing complete organellar genomes from metagenomic datasets. This workflow enables the recovery of high-quality mitochondrial and plastid genomes from environmental DNA, revealing hidden eukaryotic diversity that traditional metagenomic approaches often miss.

## Key Features

- Targeted _de novo_ assembly of organellar sequences from complex metagenomic data
- Improved assembly using a graph-based approach for resolving organellar genome structures
- Quality assessment of assembly and manual disentanglement of assembly graphs
- Support for both high and low abundance taxa
- Taxonomic identification through multi-faceted analysis

## Documentation

For detailed instructions, please refer to our Wiki pages:

- [Installation and Dependencies](https://github.com/GeertsManon/organellarMetagenomics/wiki/A:-Installation-&-Dependencies) - System requirements, dependencies, and installation guide
- [Upstream Pipeline Workflow using HPC](https://github.com/GeertsManon/organellarMetagenomics/wiki/B:-Upstream-pipeline-workflow-using-HPC) - Instructions for running the upstream pipeline (data retrieval, read trimming and assembly) on high-performance computing environments
- [Manual Disentanglement of Chloroplast Genomes from Metagenomic Data](https://github.com/GeertsManon/organellarMetagenomics/wiki/C:-Manual-Disentanglement-of-Chloroplast-Genomes-from-Metagenomic-Data-in-Bandage) - Step-by-step guide for manually reconstructing complete chloroplast genomes from assembly graphs
- [Quality Assessment of Assembled Organellar Genomes](https://github.com/GeertsManon/organellarMetagenomics/wiki/D:-Quality-Assessment-of-Assembled-Organellar-Genomes) - Description of the quality assessment process (mapping and variant calling)

## Core Components

The toolkit includes [scripts](https://github.com/GeertsManon/organellarMetagenomics/tree/main/scripts) for:

* High-throughput data acquisition from NCBI SRA
* Quality control and adapter trimming
* Targeted _de novo_ assembly with optimized parameters
* Assembly graph visualization and disentanglement
* Quality assessment through read mapping, coverage analysis and variant calling
* Circular genome visualization
