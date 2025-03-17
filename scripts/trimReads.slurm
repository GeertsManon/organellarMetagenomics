#!/bin/bash

#SBATCH --cluster=wice
#SBATCH --cpus-per-task=16
#SBATCH --time=5:00:00
#SBATCH -A l_fishgenomics_at_kul
#SBATCH -J trimmomatic
#SBATCH --mem-per-cpu=10G


############################################
#
# Trim Reads
# ==========
#
#	1. Visualise read quality before trimming
#	2. Trim reads
#	3. Visualise read quality after trimming
#
#	HOWTO: sbatch --array=1 scripts/downloadData-A.slurm
#
#	POST script: 
#		cat 2_TrimmedData/*R1_paired* > 2_TrimmedData/AveRiv_R1_paired.fastq.gz
#		cat 2_TrimmedData/*R2_paired* > 2_TrimmedData/AveRiv_R2_paired.fastq.gz
#
#		mv 2_TrimmedData/15_SRR10955151_R1_paired.fastq.gz > 2_TrimmedData/LakeMan_R1_paired.fastq.gz
#		mv 2_TrimmedData/15_SRR10955151_R2_paired.fastq.gz > 2_TrimmedData/LakeMan_R2_paired.fastq.gz
#
#		mv 2_TrimmedData/16_ERR9631049_R1_paired.fastq.gz > 2_TrimmedData/ResCRep_R1_paired.fastq.gz
#		mv 2_TrimmedData/16_ERR9631049_R2_paired.fastq.gz > 2_TrimmedData/ResCRep_R2_paired.fastq.gz
#		
############################################


###########	1	Download modules

module purge
module load cluster/wice/oldhierarchy
module load Trimmomatic
module load FastQC


###########     2	Set Variables


LOG_DATE=$(date +%Y%m%d)
export SLURM_OUTPUT="${LOG_DATE}_${SLURM_JOB_NAME}_${SLURM_ARRAY_TASK_ID}.out"
export SLURM_ERROR="${LOG_DATE}_${SLURM_JOB_NAME}_${SLURM_ARRAY_TASK_ID}.err"
exec > "$SLURM_OUTPUT" 2> "$SLURM_ERROR"


CONFIG=datasets.txt
SRR=$(awk -v ID=$SLURM_ARRAY_TASK_ID '$1==ID {print $2}' $CONFIG)
READS1=1_RawData/$SLURM_ARRAY_TASK_ID'_'$SRR'_'R1.fastq.gz
READS2=1_RawData/$SLURM_ARRAY_TASK_ID'_'$SRR'_'R2.fastq.gz
RESULTS_FILE=2_TrimmedData/results.log


###########     3       Launch command


echo ">>>>>>>   Visualising raw reads...   <<<<<<<"
fastqc --threads $SLURM_CPUS_PER_TASK $READS1 $READS2

echo ">>>>>>>   Trimming reads...   <<<<<<<"
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
	-threads $SLURM_CPUS_PER_TASK \
	-phred33 \
	$READS1 $READS2 \
	2_TrimmedData/$SLURM_ARRAY_TASK_ID'_'$SRR'_R1_paired.fastq.gz' 2_TrimmedData/$SLURM_ARRAY_TASK_ID'_'$SRR'_R1_unpaired.fastq.gz' \
	2_TrimmedData/$SLURM_ARRAY_TASK_ID'_'$SRR'_R2_paired.fastq.gz' 2_TrimmedData/$SLURM_ARRAY_TASK_ID'_'$SRR'_R2_unpaired.fastq.gz' \
	ILLUMINACLIP:adapters.fasta:2:30:10 \
	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

echo ">>>>>>>   Visualising trimmed reads...   <<<<<<<"
fastqc --threads $SLURM_CPUS_PER_TASK  2_TrimmedData/$SLURM_ARRAY_TASK_ID'_'$SRR'_R1_paired.fastq.gz' \
	2_TrimmedData/$SLURM_ARRAY_TASK_ID'_'$SRR'_R2_paired.fastq.gz'



###########     4       Results


echo '>>>>>>>   Calculating and writing result for $SRR...   <<<<<<<'
INPUT_PAIRS=$(grep -oP "Input Read Pairs: \K\d+" "$SLURM_ERROR")
BOTH_SURVIVING=$(grep -oP "Both Surviving: \K\d+" "$SLURM_ERROR")
FORWARD_ONLY=$(grep -oP "Forward Only Surviving: \K\d+" "$SLURM_ERROR")
REVERSE_ONLY=$(grep -oP "Reverse Only Surviving: \K\d+" "$SLURM_ERROR")
DROPPED=$(grep -oP "Dropped: \K\d+" "$SLURM_ERROR")

# Update results file
if [ ! -f "$RESULTS_FILE" ]; then
    echo -e "SLURM_ARRAY_TASK_ID\tINPUT_PAIRS\tBOTH_SURVIVING\tFORWARD_ONLY\tREVERSE_ONLY\tDROPPED" > "$RESULTS_FILE"
fi

if grep -q -P "^$SLURM_ARRAY_TASK_ID\t" "$RESULTS_FILE"; then
    # Update existing entry
    awk -v ID="$SLURM_ARRAY_TASK_ID" -v INPUT_PAIRS="$INPUT_PAIRS" -v BOTH_SURVIVING="$BOTH_SURVIVING" -v FORWARD_ONLY="$FORWARD_ONLY" -v REVERSE_ONLY="$REVERSE_ONLY" -v DROPPED="$DROPPED" '
        BEGIN { FS=OFS="\t" }
        $1 == ID {
            $2 = INPUT_PAIRS;
            $3 = BOTH_SURVIVING;
            $4 = FORWARD_ONLY;
            $5 = REVERSE_ONLY;
            $6 = DROPPED;
        }
        { print }
    ' "$RESULTS_FILE" > tmp && mv tmp "$RESULTS_FILE"
else
    # Append new entry
    echo -e "${SLURM_ARRAY_TASK_ID}\t${INPUT_PAIRS}\t${BOTH_SURVIVING}\t${FORWARD_ONLY}\t${REVERSE_ONLY}\t$DROPPED" >> "$RESULTS_FILE"
fi
