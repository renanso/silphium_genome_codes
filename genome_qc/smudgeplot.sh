#!/bin/bash
#SBATCH --job-name=smudge	# Job name
#SBATCH --nodes=1		    # N of nodes
#SBATCH --ntasks=16
#SBATCH --mem="128G"		# Memory per node; by default using M as unit
#SBATCH --time=200:00:00    # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --export=ALL        #export users explicit environment variables to compute node
#SBATCH --output=%x_%j.out	# Standard output log
#SBATCH --error=%x_%j.err		# Standard error log
#SBATCH --partition=normal

export TMPDIR=$JOB_TMP_DIR

eval "$(conda shell.bash hook)"
conda activate smudgeplot

# run FastK to create a k-mer database
#path1="/cluster/home/rsouza/silphium/illumina_parents/JTTM_PCRfree_1_Bad_Astra_CGGCGTGA_Silphium_I1208_L4_R[12].fastq"
#path2="/cluster/home/rsouza/silphium/genome_scripts/smudge_data/Sint"

path1="/cluster/home/rsouza/silphium/illumina_parents/JTTN_PCRfree_1_IJHW_ATGTAAGT_Silphium_I1208_L4_R[12].fastq"
path2="/cluster/home/rsouza/silphium/genome_scripts/smudge_data/Sper"

FastK -v -t4 -k31 -M16 -T4 $path1 -N$path2/FastK_Table -P$JOB_TMP_DIR

# Find all k-mer pairs in the dataset using hetmer module
smudgeplot.py hetmers -L 12 -t 16 -o $path2/kmerpairs --verbose $path2/FastK_Table

# this now generated `data/Scer/kmerpairs_text.smu` file;
# it's a flat file with three columns; covB, covA and freq (the number of k-mer pairs with these respective coverages)

# use the .smu file to infer ploidy and create smudgeplot
smudgeplot.py all -o $path2/trial_run $path2/kmerpairs_text.smu


# check that bunch files are generated (3 pdfs; some summary tables and logs)

#ls $path2/trial_run_*


######################
##example 
# download data
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR326/001/SRR3265401/SRR3265401_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR326/001/SRR3265401/SRR3265401_2.fastq.gz

# sort them in a reasonable place
#mkdir ./data/Scer
#mv *fastq.gz data/Scer/

#FastK -v -t4 -k31 -M16 -T4 data/Scer/SRR3265401_[12].fastq.gz -Ndata/Scer/FastK_Table

# Find all k-mer pairs in the dataset using hetmer module
#smudgeplot.py hetmers -L 12 -t 4 -o data/Scer/kmerpairs --verbose data/Scer/FastK_Table
# this now generated `data/Scer/kmerpairs_text.smu` file;
# it's a flat file with three columns; covB, covA and freq (the number of k-mer pairs with these respective coverages)

# use the .smu file to infer ploidy and create smudgeplot
#smudgeplot.py all -o data/Scer/trial_run data/Scer/kmerpairs_text.smu

# check that bunch files are generated (3 pdfs; some summary tables and logs)
#ls data/Scer/trial_run_*