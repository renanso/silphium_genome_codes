#!/bin/bash
#SBATCH --job-name=bcftools	# Job name
#SBATCH --nodes=1		            # N of nodes
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem="48G"			# Memory per node; by default using M as unit
#SBATCH --time=48:00:00     # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --export=ALL        #export users explicit environment variables to compute node
#SBATCH --output=%x_%j.out		# Standard output log
#SBATCH --error=%x_%j.err		# Standard error log
#SBATCH --partition=plant
#SBATCH --array 1-270 #count files with specific names for the array number: ls -lR ./*bam | wc -l


# Script to call variants from BAM files using `bcftools mpileup` and
# `bcftools call`. Intended to run as a SLURM array where each task
# processes one BAM file and produces a sorted/indexed BCF.
#
# Usage: Set the `REF` variable to the correct reference FASTA and
# run `sbatch 17.1_var_call.sh`. Adjust annotations/filters as needed.
# Dependencies: bcftools, SLURM.



# Load software modules 

ml bcftools/1.19-gcc-13.1.0
#ml cluster/vcflib/1.0.1
#ml cluster/vcftools/0.1.16

#directories and inputs
#cd /cluster/home/rsouza/003_D192_project/raw_fastqs/int_bam_files
#cd /cluster/home/rsouza/003_D192_project/raw_fastqs/per_bam_files
#cd /cluster/lab/clevenger/RSouza/004_silphium_diversity/target_seq/mapped_ba_fwd
#cd /cluster/lab/clevenger/RSouza/004_silphium_diversity/target_seq/mapped_ba_rev
#cd /cluster/home/rsouza/002_S99_project/raw_fastqs/bam_files/sorted_bam
#cd /cluster/lab/clevenger/RSouza/004_silphium_diversity/merged_1x_ts/combined/mapped
cd /cluster/lab/clevenger/RSouza/011_flexprep_silphium/updated_name/merged_scans

#REF=/cluster/home/rsouza/silphium/dual_genome/dual_genome.fasta
#REF=/cluster/home/rsouza/silphium/Silphium_perfoliatum_var_1JHW_V1_release/Silphium_perfoliatum_var_1JHW/sequences/Silphium_perfoliatum_var_1JHW.mainGenome.fasta
REF=/cluster/home/rsouza/silphium/Silphium_integrifolium_var_Bad_Astra_V1_release/Silphium_integrifolium_var_Bad_Astra/sequences/Silphium_integrifolium_var_Bad_Astra.mainGenome.fasta
#REF=/cluster/home/rsouza/silphium/Silphium_integrifolium_var_Bad_Astra_V1_release/Silphium_integrifolium_var_Bad_Astra/sequences/reverse_genome/Silphium_integrifolium_var_Bad_Astra.mainGenome_reverse.fasta
#REF=/cluster/lab/clevenger/RSouza/007_silphium50k/diagnostics/Silphium_integrifolium_var_Bad_Astra.chloroplast.fasta


BAM=$(ls *bam | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo ${BAM}

ID=$(echo ${BAM}| sed 's/.bam//')

echo ${ID}

echo "Part 1 - mpileup and call"
# part 1 - mpileup and call
#bcftools mpileup -Ou -f ${REF} ${BAM} | bcftools call -vm -Ob --format-fields GQ,GP -o ${ID}.bcf

bcftools mpileup -Ou --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR -f ${REF} ${BAM} | bcftools call -vm -Ob --format-fields GQ,GP -o ${ID}.bcf

##1b call only
#bcftools call out11.vcf.gz -vmO z --format-fields GQ,GP -o out11_called.vcf.gz
echo "Part 1 complete"

echo "Part 2 - sort and index"
# part 2 - sort and index

bcftools sort ${ID}.bcf -Ob -o  ${ID}_srt.bcf
bcftools index -c ${ID}_srt.bcf

echo "Finish sort and index"

echo "Deleting original bcf"

rm ${ID}.bcf