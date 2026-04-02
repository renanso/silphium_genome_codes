#!/bin/bash
#SBATCH --job-name=hifiasm	# Job name
#SBATCH --nodes=1		    # N of nodes
#SBATCH --ntasks=1
#SBATCH --mem="128G"		# Memory per node; by default using M as unit
#SBATCH --time=24:00:00    # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --export=ALL        #export users explicit environment variables to compute node
#SBATCH --output=%x_%j.out	# Standard output log
#SBATCH --error=%x_%j.err		# Standard error log
#SBATCH --partition=normal

#module load cluster/hifiasm/0.19.8

#data input
#hifiseq="/cluster/home/rsouza/silphium/PacBio/combined/all_data.fastq.gz"
#read1="/cluster/home/rsouza/silphium/OmniC/JQIP_OmniC_NA_NA_CTTGTCGA_Silphium_OmniC-Silphium_OmniC_I1177_L3_R1.fastq"
#read2="/cluster/home/rsouza/silphium/OmniC/JQIP_OmniC_NA_NA_CTTGTCGA_Silphium_OmniC-Silphium_OmniC_I1177_L3_R2.fastq"

## hifi data only #####
#hifiasm -o S3ZZP_hifiasm  -l0  -t 32 $hifiseq

##### with Hi-C data #####
#hifiasm -o S3ZZP_hifi.asm -t32 --h1 $read1 --h2 $read2 $hifiseq

## gfa to fasta for output
cd /cluster/home/rsouza/silphium/trio_bin
awk '/^S/{print ">"$2;print $3}' S3ZZP.asm.dip.hap1.p_ctg.gfa > S3ZZP.asm.dip.hap1.p_ctg.fa
