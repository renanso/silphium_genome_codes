#!/bin/bash
#SBATCH --job-name=FASTQC		# Job name
#SBATCH --nodes=1		        # N of nodes
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem="16G"			# Memory per node; by default using M as unit
#SBATCH --time=24:00:00     # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --export=ALL        #export user’s explicit environment variables to compute node
#SBATCH -o /cluster/lab/clevenger/RSouza/004_silphium_diversity/fastqc/target_seq_qc/slurm-%A_%a.out
#SBATCH --output=%x_%j.out		# Standard output log
#SBATCH --error=%x_%j.err		# Standard error log
#SBATCH --partition=plant
#SBATCH --array 1-540 #count files with specific names for the array number: ls -lR ./*fq.gz | wc -l
 
# Run FastQC on sequencing FASTQ files to perform per-sample
# quality control and produce per-file HTML reports. Intended to be
# submitted as a SLURM array job where each task processes one file.
#
# Usage: Edit `INPUT` and `OUTPUT` directories below, then submit with
# `sbatch 11_fastqc.sh`. This script activates a `fastqc` conda env.
# Notes: MultiQC can be run afterwards to aggregate reports.




#ml /cluster/home/rsouza/FastQC/fastqc		# Load software module
#ml cluster/multiqc/1.17 

eval "$(conda shell.bash hook)"
conda activate fastqc

#######

#input and output directories
INPUT=/cluster/lab/clevenger/RSouza/011_flexprep_silphium/updated_name/merged_scans/raw_data
OUTPUT=/cluster/lab/clevenger/RSouza/011_flexprep_silphium/updated_name/merged_scans/qc

echo "Starting FastQC Job"
#directory where FastQC is:
#cd /cluster/home/rsouza/apps/FastQC/

file=$(ls $INPUT | sed -n ${SLURM_ARRAY_TASK_ID}p)

#/cluster/home/rsouza/apps/FastQC/fastqc -o $OUTPUT $INPUT/${file}
fastqc -o $OUTPUT $INPUT/${file}
echo "Finishing FastQC Job"

#enter interactive to run
#ml cluster/multiqc/1.17
#multiqc . --interactive