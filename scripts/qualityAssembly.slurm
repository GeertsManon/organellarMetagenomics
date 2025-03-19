#!/bin/bash

#SBATCH --cluster=wice
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --account=l_fishgenomics_at_kul
#SBATCH --job-name=mapping


############################################
#
# Assess quality assembly
# =======================
#
#	      1. Generate hash index for assembled scaffold
#       2. Map trimmed paired reads against assembled scaffold
#	      3. Calculate mappings stats
#	      4. Calculate per-site read depths
#	      5. Calculate depth stats
#	      6. Call variants
#
#       HOWTO: sbatch --export=SAMPLE=AveRiv scripts/qualityAssembly.slurm
#
############################################


###########	1	Download modules


module purge
module load cluster/genius/login
module load libdeflate/1.18-GCCcore-12.3.0
module load SAMtools/1.18-GCC-12.3.0
module load SMALT/0.7.6-intel-2021a
module load BCFtools/1.18-GCC-12.3.0



###########     2	Set Variables


LOG_DATE=$(date +%Y%m%d)
export SLURM_OUTPUT="${LOG_DATE}_${SLURM_JOB_NAME}_${SAMPLE}.out"
export SLURM_ERROR="${LOG_DATE}_${SLURM_JOB_NAME}_${SAMPLE}.err"
exec > "$SLURM_OUTPUT" 2> "$SLURM_ERROR"

if [ -z "$SAMPLE" ]; then
    echo "Sample not set, abort script."
    exit 1
fi

DIR=4_QualityAssembly/$SAMPLE
RES_FILE=$DIR/results.log
READS1=2_TrimmedData/${SAMPLE}_R1_paired.fastq.gz
READS2=2_TrimmedData/${SAMPLE}_R2_paired.fastq.gz
REF=$(ls "$DIR"/*fasta)
INDEX="${REF%.fasta}"

[ -d 4_QualityAssembly ] || mkdir 4_QualityAssembly
[ -d $DIR ] || mkdir $DIR




###########     3	Generate hash index for ref


smalt index -k 13 -s 1 $INDEX $REF
samtools faidx $REF


IDS=(80 85 90 95 96 97 98 99)

for ID in "${IDS[@]}"; do


###########     4       Map reads & convert to bam


	echo "Running SMALT mapping with identity threshold: $ID"

	smalt map -y 0.$ID -x -n $SLURM_NTASKS_PER_NODE -o ${INDEX}_id${ID}.sam $INDEX $READS1 $READS2
	samtools sort -O bam -@ $SLURM_NTASKS_PER_NODE -o ${INDEX}_id${ID}_sorted.bam ${INDEX}_id${ID}.sam
	samtools view -@ $SLURM_NTASKS_PER_NODE -b -F 4 ${INDEX}_id${ID}_sorted.bam > ${INDEX}_id${ID}_mapped.bam
	samtools index -@ $SLURM_NTASKS_PER_NODE ${INDEX}_id${ID}_mapped.bam


###########     5       Calculate mappings stats


	BAM=${INDEX}_id${ID}_sorted.bam
	N_READS=$(samtools view -@ $SLURM_NTASKS_PER_NODE -c $BAM)
	N_READS_UNMAPPED=$(samtools view -@ $SLURM_NTASKS_PER_NODE -c -f 4 $BAM)
	N_READS_MAPPED=$(samtools view -@ $SLURM_NTASKS_PER_NODE -c -F 4 $BAM)
	rm ${INDEX}_id${ID}.sam ${INDEX}_id${ID}_sorted.bam


###########     6       Calculate read depth


	DEPTH=${INDEX}_id${ID}_depth.txt
	samtools depth -a ${INDEX}_id${ID}_mapped.bam > $DEPTH


###########     7       Caluclate depth stats


	MEDIAN=$(cat $DEPTH | awk '{print $3}' | sort -n | awk '{a[i++]=$1} END {print a[int(i/2)]}')
	awk '{sum+=$3; sumsq+=$3^2} END {mean=sum/NR; print mean, sqrt(sumsq/NR - mean^2)}' $DEPTH > 4_QualityAssembly/$SAMPLE/DEPTH_STATS
	MEAN=$(awk '{print $1}' 4_QualityAssembly/$SAMPLE/DEPTH_STATS)
	SD=$(awk '{print $2}' 4_QualityAssembly/$SAMPLE/DEPTH_STATS)
	rm 4_QualityAssembly/$SAMPLE/DEPTH_STATS

	covered_bases=$(samtools mpileup ${INDEX}_id${ID}_mapped.bam | awk -v X=1 '$4>=X' | wc -l)
        ref_length=$(awk '{L+=$2} END {print L}' $REF.fai)
        BREADTH=$(echo "scale=2; ($covered_bases / $ref_length) * 100" | bc)


###########     8	Call variants



	bcftools mpileup -Ou --skip-indels -f $REF -a "INFO/AD" --threads "$SLURM_NTASKS_PER_NODE" ${INDEX}_id${ID}_mapped.bam | \
    		bcftools call -Oz -mv -o ${INDEX}_id${ID}.vcf.gz --threads "$SLURM_NTASKS_PER_NODE"

	bcftools view -i 'DP>=10 && QUAL>=20' ${INDEX}_id${ID}.vcf.gz -o ${INDEX}_id${ID}_filtered.vcf.gz
	#rm  ${INDEX}_id${ID}.vcf.gz

	# Extract heterozygous (0/1) and homozygous (1/1) variant counts
	HET=$(gzip -cd ${INDEX}_id${ID}_filtered.vcf.gz | grep -v '^#' | grep -c "0/1")
	HOM=$(gzip -cd ${INDEX}_id${ID}_filtered.vcf.gz | grep -v '^#' | grep -c "1/1")




###########     9       Output results


	# Check if the file exists
	if [ ! -f $RES_FILE ]; then
        	echo -e "SAMPLE\tID\tN_READS\tN_READS_UNMAPED\tN_READS_MAPPED\tMEAN\tMEDIAN\tSD\tBREADTH\tHET\tHOM" > $RES_FILE
	fi

	echo -e "$SAMPLE\t$ID\t${N_READS}\t${N_READS_UNMAPPED}\t${N_READS_MAPPED}\t$MEAN\t$MEDIAN\t$SD\t$BREADTH\t$HET\t$HOM" >> $RES_FILE


done
