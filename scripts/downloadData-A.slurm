#!/bin/bash

#SBATCH --cluster=wice
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --account=l_fishgenomics_at_kul
#SBATCH --job-name=downloadData-A

############################################
#
# Downloading Data from NCBI SRA
# ==============================
#
#	  1. Download raw sequencing data in its native .sra format (single core)
#	  2. Launch script to retrieve .sra files and store them locally (multiple cores)
#	  3. Launch script for trimming reads (multiple cores)
#
#  	HOWTO: sbatch --array=1 scripts/downloadData-A.slurm
#
############################################


###########	1	Download modules


module purge
module load cluster/wice/batch
module load SRA-Toolkit


###########     2	Set Variables


CONFIG=datasets.txt
SRR=$(awk -v ID=$SLURM_ARRAY_TASK_ID '$1==ID {print $2}' $CONFIG)

LOG_DATE=$(date +%Y%m%d)
export SLURM_OUTPUT="${LOG_DATE}_${SLURM_ARRAY_TASK_ID}_${SRR}.out"
export SLURM_ERROR="${LOG_DATE}_${SLURM_ARRAY_TASK_ID}_${SRR}.err"
exec > "$SLURM_OUTPUT" 2> "$SLURM_ERROR"

[ -d 1_RawData ] || mkdir 1_RawData


###########     3	 Launch command (single core)


echo ">>>>>>>   Prefetching $SRR...   <<<<<<<"
prefetch $SRR --progress yes -O 1_RawData


###########     4       Launch fastq-dump (multiple cores)


echo ">>>>>>>   Launched downloadData-B.slurm   <<<<<<<"
sbatch --array=$SLURM_ARRAY_TASK_ID scripts/downloadData-B.slurm
