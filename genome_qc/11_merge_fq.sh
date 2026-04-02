#!/bin/bash
#SBATCH --job-name=merge_fq	# Job name
#SBATCH --nodes=1		            # N of nodes
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem="48G"			# Memory per node; by default using M as unit
#SBATCH --time=48:00:00     # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --export=ALL        #export users explicit environment variables to compute node
#SBATCH --output=%x_%j.out	# Standard output log
#SBATCH --error=%x_%j.err		# Standard error log
#SBATCH --partition=plant

file1=/cluster/lab/clevenger/RSouza/007_silphium50k/fastqs/stax/raw_data/og_name
file2=/cluster/lab/clevenger/RSouza/011_flexprep_silphium/updated_name/fastqs
combined=/cluster/lab/clevenger/RSouza/011_flexprep_silphium/updated_name/phylogeny

for i in $(ls $file1); do cat $file1/$i $file2/$i > $combined/$i; done


#cat /cluster/lab/clevenger/RSouza/004_silphium_diversity/merged_1x_ts/combined/fastqs/E2.63.1_R1.fq.gz  /cluster/lab/clevenger/RSouza/004_silphium_diversity/target_seq/fastqs/E2.63.1_R1.fq.gz > E2.63.1_R1.fq.gz
#cat /cluster/lab/clevenger/RSouza/004_silphium_diversity/merged_1x_ts/combined/fastqs/E2.63.1_R2.fq.gz  /cluster/lab/clevenger/RSouza/004_silphium_diversity/target_seq/fastqs/E2.63.1_R2.fq.gz > E2.63.1_R2.fq.gz 