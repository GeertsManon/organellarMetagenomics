#!/bin/bash

#SBATCH --cluster=wice
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --account=l_fishgenomics_at_kul
#SBATCH --job-name=finetune
#SBATCH --mem-per-cpu=5G


############################################
#
# Assemble chloroplast genomes
# ============================
#
#       1. De novo assembly cpDNA using GetOrganelle
#
#       HOWTO: sbatch --export=SAMPLE=AveRiv,WS=60 scripts/finetuneAssembly.slurm
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
export SLURM_OUTPUT="${LOG_DATE}_${SLURM_JOB_NAME}_${SAMPLE}_${WS}.out"
export SLURM_ERROR="${LOG_DATE}_${SLURM_JOB_NAME}_${SAMPLE}_${WS}.err"
exec > "$SLURM_OUTPUT" 2>> "$SLURM_ERROR"

SEED=customDB/Cyclotella-n8.seed.fasta
LABEL=customDB/Cyclotella-n8.label.fasta

READS1=2_TrimmedData/${SAMPLE}_R1_paired.fastq.gz
READS2=2_TrimmedData/${SAMPLE}_R2_paired.fastq.gz




###########     3	Launch command


R=40
KMER=21,33,55,77,99,127
# default kmers: 21,55,85,115

if [ -z "$R" ]; then
	echo "R is not set, using default value"
	if [ -z "$WS" ]; then
		echo "WS is not set, using default value"
		DIR="3_Assembly/${SAMPLE}_Rdef_DBn8_WSdef"
	else
		echo "WS is set to: $WS"
		DIR="3_Assembly/${SAMPLE}_Rdef_DBn8_WS${WS}"
	fi
else
	echo "R is set to: $R"
	if [ -z "$WS" ]; then
		echo "WS is not set, using default value"
		DIR="3_Assembly/${SAMPLE}_R${R}_DBn8_WSdef"
	else
		echo "WS is set to: $WS"
		DIR="3_Assembly/${SAMPLE}_R${R}_DBn8_WS${WS}"
	fi
fi


[ ! -d "$DIR" ] && mkdir "$DIR"

echo ">>>>>>>   Assembling data for $SAMPLE...   <<<<<<<"
get_organelle_from_reads.py -1 $READS1 -2 $READS2 \
	-o $DIR --overwrite -t 32 --memory-save \
	-s $SEED --genes $LABEL -F other_pt \
	-k $KMER ${R:+-R $R} ${WS:+-w $WS} \
	--reduce-reads-for-coverage inf \
       	--max-reads inf



conda deactivate
