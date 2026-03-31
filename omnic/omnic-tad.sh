#!/bin/bash
#SBATCH --job-name=omnic-tad	# Job name
#SBATCH --nodes=1		    # N of nodes
#SBATCH --ntasks=32
#SBATCH --mem="128G"		# Memory per node; by default using M as unit
#SBATCH --time=200:00:00    # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --export=ALL        #export users explicit environment variables to compute node
#SBATCH --output=%x_%j.out	# Standard output log
#SBATCH --error=%x_%j.err		# Standard error log
#SBATCH --partition=normal

module load cluster/bedtools/2.28.0
module load cluster/bwa/0.7.17
module load cluster/pairtools
module load cluster/samtools/1.16.1
module load cluster/preseq/2.0
module load cluster/python/3.8.10

eval "$(conda shell.bash hook)"
conda activate omnic

#genome input
maternal="/cluster/home/rsouza/silphium/Silphium_integrifolium_var_Bad_Astra_V1_release/Silphium_integrifolium_var_Bad_Astra/sequences/Silphium_integrifolium_var_Bad_Astra.mainGenome.fasta"
paternal="/cluster/home/rsouza/silphium/Silphium_perfoliatum_var_1JHW_V1_release/Silphium_perfoliatum_var_1JHW/sequences/Silphium_perfoliatum_var_1JHW.mainGenome.fasta"

# If too many contigs, filter the genome to remove contigs < 10kb
#module load cluster/bbmap/38.42
#reformat.sh in=S3ZZP.asm.dip.hap1.p_ctg.fa out=S3ZZP.asm.dip.hap1.p_ctg_1kb_fltr.fa minlength=10000 -Xmx64g

#omnic input

R1="/cluster/home/rsouza/silphium/OmniC/JQIP_OmniC_NA_NA_CTTGTCGA_Silphium_OmniC-Silphium_OmniC_I1177_L3_R1.fastq"
R2="/cluster/home/rsouza/silphium/OmniC/JQIP_OmniC_NA_NA_CTTGTCGA_Silphium_OmniC-Silphium_OmniC_I1177_L3_R2.fastq"

#1. Generate an index file for your reference

#samtools faidx $maternal
#samtools faidx $paternal

#2. Use the index file to generate the genome file by printing the first two columns into a new file.

#cut -f1,2 ${maternal}.fai > ${maternal}.genome
#cut -f1,2 ${paternal}.fai > ${paternal}.genome

# 3. Generate a bwa index file for the chosen reference.

#bwa index $maternal
#bwa index $paternal

# 4. Making a temp folder to store intermediate files

#mkdir ./temp

#5. From fastq to final valid pairs bam file

bwa mem -5SP -T0 -t16 $paternal $R1 $R2 | \
pairtools parse --min-mapq 40 --walks-policy 5unique \
--max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path $paternal | \
pairtools sort --tmpdir=/cluster/home/rsouza/silphium/genome_scripts/tad/temp --nproc 16|pairtools dedup --nproc-in 8 \
--nproc-out 8 --mark-dups --output-stats stats.txt|pairtools split --nproc-in 8 \
--nproc-out 8 --output-pairs mapped.pairs --output-sam -|samtools view -bS -@16 | \
samtools sort -@16 -o mapped.PT.bam;samtools index -c mapped.PT.bam

