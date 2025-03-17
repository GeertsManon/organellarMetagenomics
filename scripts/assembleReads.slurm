#!/bin/bash

#SBATCH --cluster=wice
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --account=l_fishgenomics_at_kul
#SBATCH --job-name=GetOrganelle


############################################
#
# Assemble chloroplast genomes
# ============================
#
#       1. De novo assembly of cpDNA using GetOrganelle
#
#       HOWTO: sbatch --export=SAMPLE=AveRiv scripts/assembleReads.slurm
#
############################################



###########	1	Download modules


module purge
module load cluster/wice/batch
module load SPAdes
module load BLAST+
module load Python

source /data/leuven/347/vsc34774/miniconda3/etc/profile.d/conda.sh
conda activate GetOrganelle




###########     2	Set Variables



LOG_DATE=$(date +%Y%m%d)
export SLURM_OUTPUT="${LOG_DATE}_${SLURM_JOB_NAME}_${SAMPLE}.out"
export SLURM_ERROR="${LOG_DATE}_${SLURM_JOB_NAME}_${SAMPLE}.err"
exec > "$SLURM_OUTPUT" 2>> "$SLURM_ERROR"

if [ -z "$SAMPLE" ]; then
    echo "Sample not set, abort script."
    exit 1
fi

DIR=3_Assembly/$SAMPLE
READS1=2_TrimmedData/${SAMPLE}_R1_paired.fastq.gz
READS2=2_TrimmedData/${SAMPLE}_R2_paired.fastq.gz

[ -d 3_Assembly ] || mkdir 3_Assembly
[ -d $DIR ] || mkdir $DIR


###########     3	Launch command


echo ">>>>>>>   Assembling data for $SAMPLE...   <<<<<<<"
get_organelle_from_reads.py -1 $READS1 -2 $READS2 \
	-o $DIR --overwrite -R 40 \
	-k 21,33,55,77,99,127 \
	-F other_pt --memory-save -t 32  2>> "$SLURM_ERROR"


conda deactivate
