#!/bin/bash
#SBATCH --job-name=oatk	# Job name
#SBATCH --nodes=1		    # N of nodes
#SBATCH --ntasks=8
#SBATCH --mem="256G"		# Memory per node; by default using M as unit
#SBATCH --time=240:00:00    # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --export=ALL        #export users explicit environment variables to compute node
#SBATCH --output=%x_%j.out	# Standard output log
#SBATCH --error=%x_%j.err		# Standard error log
#SBATCH --partition=normal

#conda
eval "$(conda shell.bash hook)"
conda activate oatk


oatk -k 1001 -c 30 -t 8 --nhmmscan /cluster/home/rsouza/anaconda3/envs/oatk/bin/nhmmscan -m angiosperm_mito.fam -p angiosperm_pltd.fam -o bad_astra /cluster/home/rsouza/silphium/trio_bin/drafts/S3ZZP.asm.dip.hap2.p_ctg.fa

#oatk -k 1001 -c 30 -t 8 --nhmmscan /cluster/home/rsouza/anaconda3/envs/oatk/bin/nhmmscan -m angiosperm_mito.fam -p angiosperm_pltd.fam -o 1jhw /cluster/home/rsouza/silphium/trio_bin/drafts/S3ZZP.asm.dip.hap1.p_ctg.fa
