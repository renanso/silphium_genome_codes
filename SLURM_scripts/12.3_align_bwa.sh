#!/bin/bash
#SBATCH --job-name=map-reference	# Job name
#SBATCH --nodes=1		            # N of nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem="64G"			# Memory per node; by default using M as unit
#SBATCH --time=48:00:00     # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --export=ALL        #export users explicit environment variables to compute node
#SBATCH --output=./logs/%x_%j.out	# Standard output log
#SBATCH --error=./logs/%x_%j.err		# Standard error log
#SBATCH --partition=normal
#SBATCH --array 1-270  #count files with specific names for the array number: ls -lR ./*_R1.fq.gz | wc -l


# Script to map paired-end reads to a reference genome using BWA MEM,
# convert to BAM, sort and index with samtools. Designed to run as a
# SLURM array where each task maps one sample (paired R1/R2 files).
# Usage: Set `REF` and BED files as needed, run with `sbatch
# 12.3_align_bwa.sh`. The script expects R1/R2 files in the working
# directory and uses `SLURM_ARRAY_TASK_ID` to select the sample.
# Dependencies: BWA, samtools; the script loads modules via `ml`.


# Load software modules 

ml cluster/bwa/0.7.17
ml samtools/1.19.2-gcc-13.1.0

#directories and inputs
cd /cluster/lab/clevenger/RSouza/011_flexprep_silphium/updated_name/merged_scans/raw_data

REF=/cluster/home/rsouza/silphium/ref_index/bwa_index/Silphium_integrifolium_var_Bad_Astra.mainGenome.fasta
#REF=/cluster/lab/clevenger/RSouza/007_silphium50k/diagnostics/blotch1.fasta
#REF=/cluster/lab/clevenger/RSouza/007_silphium50k/diagnostics/DCMV-KS4.fasta
#REF=/cluster/lab/clevenger/RSouza/007_silphium50k/diagnostics/rust.fasta
#REF=/cluster/lab/clevenger/RSouza/007_silphium50k/diagnostics/Silphium_integrifolium_var_Bad_Astra.chloroplast.fasta

BED1=/cluster/lab/clevenger/RSouza/007_silphium50k/scripts/aux_data/bed_50k_target.txt
#BED2=/cluster/lab/clevenger/RSouza/007_silphium50k/scripts/aux_data/bed_50k_bind.txt
#BED3=/cluster/lab/clevenger/RSouza/007_silphium50k/scripts/aux_data/bed_50k_offtarget.txt
#BED4=/cluster/lab/clevenger/RSouza/007_silphium50k/scripts/blotch_target.bed

FILE=$(ls *R1.* | sed -n ${SLURM_ARRAY_TASK_ID}p)

R1=$(ls *R1.* | sed -n ${SLURM_ARRAY_TASK_ID}p)
echo ${R1}

R2=$(echo ${R1} | sed "s/R1/R2/")
echo ${R2}

ID=$(echo ${FILE}| sed 's/_R[1,2].fq.gz//')
echo ${ID}

RG="@RG\tID:silphium\tSM:${ID}\tPL:${ILLUMINA}\tLB:${ID}_1"
# running

echo "Starting mapping"

bwa mem -M -t 8 -R ${RG} ${REF} ${R1} ${R2} | samtools view -bS -@ 8 -o ${ID}.bam

echo "Finishing mapping"

echo "Starting samtools sort and index"

samtools sort -@ 8 ${ID}.bam -o ${ID}_srt.bam

samtools index -c -@ 8 ${ID}_srt.bam

echo "Finish samtools sort and index"

echo "Starting samtools depth calculation"
#samtools depth -b ${BED1} ${ID}_srt.bam > ${ID}_srt.target.depth
#samtools depth -b ${BED1} ${ID}_srt.bam | wc -l | awk '{print $1}' > ${ID}_srt.target_count.depth
#samtools depth -b ${BED1} ${ID}_srt.bam |  awk '{sum+=$3} END { print sum/NR}' > ${ID}_srt.target_mean.depth

#samtools depth -b ${BED2} ${ID}_srt.bam > ${ID}_srt.bind.depth
#samtools depth -b ${BED2} ${ID}_srt.bam | wc -l | awk '{print $1}' > ${ID}_srt.bind_count.depth
#samtools depth -b ${BED2} ${ID}_srt.bam |  awk '{sum+=$3} END { print sum/NR}' > ${ID}_srt.bind_mean.depth

#samtools depth -b ${BED3} ${ID}_srt.bam > ${ID}_srt.offtarget.depth
#samtools depth -b ${BED3} ${ID}_srt.bam | wc -l | awk '{print $1}' > ${ID}_srt.offtarget_count.depth
#samtools depth -b ${BED3} ${ID}_srt.bam |  awk '{sum+=$3} END { print sum/NR}' > ${ID}_srt.offtarget_mean.depth

#samtools depth -b ${BED4} ${ID}_srt.bam > ${ID}_srt.target.depth
#samtools depth -b ${BED4} ${ID}_srt.bam | wc -l | awk '{print $1}' > ${ID}_srt.bind_count.depth
#samtools depth -b ${BED4} ${ID}_srt.bam |  awk '{sum+=$3} END { print sum/NR}' > ${ID}_srt.bind_mean.depth

#samtools depth ${ID}_srt.bam -o ${ID}_srt.depth
#samtools coverage ${ID}_srt.bam -o ${ID}_srt.coverage

echo "finish depth calculation"

echo "Deleting original bam"

rm ${ID}.bam

#mkdir ./blotch_aln
#mv *srt* blotch_aln