#!/bin/bash
#SBATCH --job-name=homer	# Job name
#SBATCH --nodes=1		    # N of nodes
#SBATCH --ntasks=8
#SBATCH --mem="256G"		# Memory per node; by default using M as unit
#SBATCH --time=144:00:00    # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --export=ALL        #export users explicit environment variables to compute node
#SBATCH --output=%x_%j.out	# Standard output log
#SBATCH --error=%x_%j.err		# Standard error log
#SBATCH --partition=normal

ml cluster/bowtie/2.3.5.1
ml samtools/1.19.2-gcc-13.1.0

eval "$(conda shell.bash hook)"
conda activate base

maketag=/cluster/home/rsouza/apps/homer/bin/makeTagDirectory

#using bams
#outputdir="/cluster/home/rsouza/silphium/genome_scripts/tad/bad_astra/homer"
#alignment="/cluster/home/rsouza/silphium/genome_scripts/tad/bad_astra/mapped.PT.bam"
#makeTagDirectory <OutputTagDirectory> <alignment file-read1>,<alignment file-read -tbp 1
#$maketag $outputdir $alignment -tbp 1


### individual alignments

#raw reads
REF=/cluster/home/rsouza/silphium/ref_index/bowtie_index/silphium_bad_astra
REF2=/cluster/home/rsouza/silphium/ref_index/bowtie_index/silphium_1jhw
R1="/cluster/home/rsouza/silphium/OmniC/JQIP_OmniC_NA_NA_CTTGTCGA_Silphium_OmniC-Silphium_OmniC_I1177_L3_R1.fastq"
R2="/cluster/home/rsouza/silphium/OmniC/JQIP_OmniC_NA_NA_CTTGTCGA_Silphium_OmniC-Silphium_OmniC_I1177_L3_R2.fastq"
outputdir="/cluster/home/rsouza/silphium/genome_scripts/tad/bad_astra/homer3"
outputdir2="/cluster/home/rsouza/silphium/genome_scripts/tad/1jhw/homer"

#bowtie2 -p 20 -x $REF -U $R1 > omniC_R1.ba.sam
#bowtie2 -p 20 -x $REF -U $R2 > omniC_R2.ba.sam

#bowtie2 -p 20 -x $REF2 -U $R1 | samtools view -bS -@ 8 -o omniC_R1.1jhw.bam
#bowtie2 -p 20 -x $REF2 -U $R2 | samtools view -bS -@ 8 -o omniC_R2.1jhw.bam

$maketag $outputdir omniC_R1.ba.bam,omniC_R2.ba.bam -tbp 1
#$maketag $outputdir2 omniC_R1.1jhw.bam,omniC_R2.1jhw.bam -tbp 1
