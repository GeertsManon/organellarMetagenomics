#!/bin/bash

source /data/leuven/347/vsc34774/miniconda3/etc/profile.d/conda.sh
conda activate GetOrganelle

# Check if the script has exactly two positional arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file.gb> <input_file.fasta>"
    exit 1
fi


# Extract annotated regions from the .gb file
get_annotated_regions_from_gb.py "$1" -o "${1%.gb}_CDS" -t CDS --mix

# Move the gene.fasta file to a new location with a new name
mv "${1%.gb}_CDS/gene/gene.fasta" "${1%.gb}.label.fasta"

# Check if the file was moved successfully and echo a message
if [ -f "${1%.gb}.label.fasta" ]; then
    echo "File has been successfully created: ${1%.gb}.label.fasta"
else
    echo "Error: ${1%.gb}.label.fasta was not created."
fi

# Move the second file and rename it
mv "$2" "${2%.fasta}.seed.fasta"

# Check if the file was moved successfully and echo a message
if [ -f "${2%.fasta}.seed.fasta" ]; then
    echo "File has been successfully created: ${2%.fasta}.seed.fasta"
else
    echo "Error: ${2%.fasta}.seed.fasta was not created."
fi

conda deactivate
