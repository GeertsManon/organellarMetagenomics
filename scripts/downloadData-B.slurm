#!/bin/bash

#SBATCH --cluster=wice
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --account=l_fishgenomics_at_kul
#SBATCH --job-name=downloadData-B


############################################
#
# Downloading Data from NCBI SRA
# ==============================
#
#       1. Download raw sequencing data in its native .sra format (single core)
#       2. Launch script to retrieve .sra files and store them locally (multiple cores)
#	      3. Launch script for trimming reads (multiple cores)
#
#       HOWTO: sbatch --array=1 scripts/downloadData-A.slurm
#
############################################




###########     1       Download modules

module purge
module load cluster/wice/batch
module load SRA-Toolkit



###########     2       Set Variables


LOG_DATE=$(date +%Y%m%d)
export SLURM_OUTPUT="${LOG_DATE}_${SLURM_JOB_NAME}_${SLURM_ARRAY_TASK_ID}.out"
export SLURM_ERROR="${LOG_DATE}_${SLURM_JOB_NAME}_${SLURM_ARRAY_TASK_ID}.err"
exec > "$SLURM_OUTPUT" 2> "$SLURM_ERROR"


CONFIG=samples
SRR=$(awk -v ID=$SLURM_ARRAY_TASK_ID '$1==ID {print $2}' $CONFIG)

[ -d 2_TrimmedData ] || mkdir 2_TrimmedData



###########     3       Launch command


echo ">>>>>>>   Downloading data for $SRR...   <<<<<<<"
fasterq-dump $SRR -O 1_RawData -o $SLURM_ARRAY_TASK_ID"_"$SRR \
	--threads $SLURM_NTASKS_PER_NODE --progress

mv 1_RawData/$SLURM_ARRAY_TASK_ID"_"$SRR"_1.fastq" 1_RawData/$SLURM_ARRAY_TASK_ID"_"$SRR"_R1.fastq"
mv 1_RawData/$SLURM_ARRAY_TASK_ID"_"$SRR"_2.fastq" 1_RawData/$SLURM_ARRAY_TASK_ID"_"$SRR"_R2.fastq"

echo ">>>>>>>   gzipping fastq files...   <<<<<<<"
pigz -c -p $SLURM_NTASKS_PER_NODE -v 1_RawData/$SLURM_ARRAY_TASK_ID"_"$SRR"_R1.fastq" > 1_RawData/$SLURM_ARRAY_TASK_ID"_"$SRR"_R1.fastq.gz"
pigz -c -p $SLURM_NTASKS_PER_NODE -v 1_RawData/$SLURM_ARRAY_TASK_ID"_"$SRR"_R2.fastq" > 1_RawData/$SLURM_ARRAY_TASK_ID"_"$SRR"_R2.fastq.gz"

if [ -f 1_RawData/$SLURM_ARRAY_TASK_ID"_"$SRR"_R1.fastq.gz" ]; then
        rm -f 1_RawData/$SLURM_ARRAY_TASK_ID"_"$SRR"_R1.fastq"
fi

if [ -f 1_RawData/$SLURM_ARRAY_TASK_ID"_"$SRR"_R2.fastq.gz" ]; then
	rm -f 1_RawData/$SLURM_ARRAY_TASK_ID"_"$SRR"_R2.fastq"
fi


###########     4       Results


echo ">>>>>>>   Calculating and writing result for $SRR...   <<<<<<<"
# Extract the number of reads read and written from the SLURM_ERROR file using awk
READS_READ=$(awk '/reads read/ {print $NF}' "$SLURM_ERROR")
READS_WRITTEN=$(awk '/reads written/ {print $NF}' "$SLURM_ERROR")

# Check if the file exists
if [ ! -f 1_RawData/results.log ]; then
    # Create the file and add a header (optional)
    echo -e "SLURM_ARRAY_TASK_ID\tSRR\tREADS_READ\tREADS_WRITTEN" > 1_RawData/results.log
fi

if grep -q -P "^$SLURM_ARRAY_TASK_ID\t" 1_RawData/results.log; then
    # Update the existing line if ID is found
    awk -v ID="$SLURM_ARRAY_TASK_ID" -v READS_READ="$READS_READ" -v READS_WRITTEN="$READS_WRITTEN" '
        BEGIN { FS=OFS="\t" }
        $1 == ID {
            $3 = READS_READ;
            $4 = READS_WRITTEN;
        }
        { print }
    ' 1_RawData/results.log > tmp && mv tmp 1_RawData/results.log
else
    # Append a new line if ID is not found
    echo -e "$SLURM_ARRAY_TASK_ID\t$SRR\t$READS_READ\t$READS_WRITTEN" >> 1_RawData/results.log
fi




###########     5	Launch trimming


echo ">>>>>>>   Launched trimReads.slurm   <<<<<<<"
sbatch --array=$SLURM_ARRAY_TASK_ID scripts/trimReads.slurm
