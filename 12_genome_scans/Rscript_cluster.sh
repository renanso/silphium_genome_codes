#!/bin/bash
#SBATCH --job-name=gwas		# Job name
#SBATCH --nodes=1		        # N of nodes
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem="128G"			# Memory per node; by default using M as unit
#SBATCH --time=72:00:00              	# Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --export=NONE                   # Do not export any user’s explicit environment variables to compute node
#SBATCH --output=%x_%j.out		# Standard output log
#SBATCH --error=%x_%j.err		# Standard error log
#SBATCH --partition=normal

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

eval "$(conda shell.bash hook)"
conda activate r_env

R CMD BATCH gwas_gapit.R

# Need to run a conversion dos2unix before running the script if editing on a windows machine. sbatch: error: Batch script contains DOS line breaks (\r\n) instead of expected UNIX line breaks (\n)
#cp ./gwas_gapit_v4.R ./gwas_gapit_v4.bkp; dos2unix ./gwas_gapit_v4.R



