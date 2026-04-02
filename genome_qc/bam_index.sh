#!/bin/bash
#SBATCH --job-name=b_index	# Job name
#SBATCH --nodes=1		    # N of nodes
#SBATCH --ntasks=8
#SBATCH --mem="96G"		# Memory per node; by default using M as unit
#SBATCH --time=72:00:00    # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --export=ALL        #export users explicit environment variables to compute node
#SBATCH --output=%x_%j.out	# Standard output log
#SBATCH --error=%x_%j.err		# Standard error log
#SBATCH --partition=normal

ml cluster/samtools/1.16.1 

#samtools index -c -@ 8 /cluster/home/rsouza/silphium/genome_scripts/tad/1jhw/mapped.PT.bam
samtools index -c -@ 8 /cluster/home/rsouza/silphium/genome_scripts/tad/bad_astra/mapped.PT.bam