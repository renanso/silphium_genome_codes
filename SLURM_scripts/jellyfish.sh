#!/bin/bash
#SBATCH --job-name=kmer_fq # Job name
#SBATCH --nodes=1		     # N of nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 10
#SBATCH --mem="1000G"		# Memory per node; by default using M as unit
#SBATCH --time=48:00:00    # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --export=ALL        #export users explicit environment variables to compute node
#SBATCH --output=%x_%j.out	# Standard output log
#SBATCH --error=%x_%j.err		# Standard error log
#SBATCH --partition=highmem


# Purpose: Run Jellyfish k-mer counting and produce a k-mer histogram
# for downstream analyses (e.g., genome size estimation, heterozygosity
# assessment). Adjust `-m` (k-mer size), `-s` (hash size) and thread
# counts to fit your data and memory limits.
# Usage: Update `jellyfish` binary path and `reads` location, then
# submit with `sbatch jellyfish.sh`. The script runs `jellyfish histo`
# to output a histogram file for plotting or further processing.



jellyfish="/cluster/home/rsouza/apps/jellyfish/jellyfish-2.3.1/bin/jellyfish"
reads="/cluster/home/rsouza/silphium/illumina_parents"

#$jellyfish count -C -m 21 -s 1000000000 -t 10 $reads/JTTN_PCRfree_1_IJHW_ATGTAAGT_Silphium_I1208_L4_R*.fastq -o reads_ijhw.jf

#Plot the histogram:

cd /cluster/home/rsouza/silphium/genome_scripts
$jellyfish histo -t 10 reads_ijhw.jf > reads_ijhw.histo

#Plot the histogram with an x axis of one million instead of default 10,000:
#jellyfish histo -t 10 --high=1000000 mer_counts.jf > reads.histo
