################################################################################
# Genome-Wide Association Study (GWAS) using GAPIT3
################################################################################
# Purpose: Performs genome-wide association analysis on SNP data and phenotypes
#          using two complementary GWAS models (FarmCPU and BLINK)
# Input: 
#   - SNP data: 011_silphium_smiss1_miss03_maf01_449K_258n_updated_2n_std.Ihapmap
#   - Phenotypes: rfc.txt (Ray Floret Count and other traits)
# Output: GAPIT results, Manhattan plots, Q-Q plots
# Reference: GAPIT3 package (https://github.com/jiabowang/GAPIT3)
################################################################################

# Clear workspace and reset options
rm(list = ls())

# Note: Installation code below (commented out after initial setup)
# devtools::install_github("jiabowang/GAPIT3", force=TRUE)
# options(download.file.method="libcurl", url.method="libcurl")
# source("http://www.zzlab.net/GAPIT/GAPIT.library.R")
# source("http://www.zzlab.net/GAPIT/gapit_functions.txt")

# ============================================================================
# LOAD REQUIRED LIBRARIES
# ============================================================================
library("data.table")           # Fast reading and manipulation of large dataframes
library("GAPIT")                # Genome Association and Prediction Integrated Tool
library("BiocParallel")         # Parallel computation for Bioconductor

# Test parallel processing with 8 cores
bplapply(1:10, print, BPPARAM = MulticoreParam(workers = 8))

# ============================================================================
# SECTION 1: READ AND QUALITY CONTROL GENOTYPING DATA
# ============================================================================
# Read SNP genotype data in Hapmap format
# File contains 449K SNPs for 258 samples (after QC filtering)
# QC filters applied: missing data <1%, <3%, and MAF >0.01
myG <- fread(file = "011_silphium_smiss1_miss03_maf01_449K_258n_updated_2n_std.Ihapmap", 
             sep = "\t")

# Remove header row (contains SNP metadata, not actual genotypes)
myG <- myG[2:nrow(myG), ]

# Display structure of genotype matrix
str(myG)

# ============================================================================
# SECTION 2: DATA FILTERING AND FORMATTING
# ============================================================================
# Filter SNPs: Keep only biallelic markers (alleles field has ≤3 characters)
# Removes multiallelic sites which can complicate analysis
myG <- subset(myG, nchar(as.character(alleles)) <= 3)

# Split genotype data into metadata and genotype calls
# Columns 1-11: SNP metadata (rs#, alleles, chromosome, position, etc.)
myG1 <- myG[, 1:11]
# Columns 12+: Genotype calls for each sample (0/1/2)
myG2 <- myG[, 12:ncol(myG)]

# Examine sample of data
myG2[1:20, 1:30]

# ============================================================================
# SECTION 3: READ AND PROCESS PHENOTYPIC DATA
# ============================================================================
# Load phenotype file (RF = Ray Floret count, etc.)
myY <- fread("rfc.txt", sep = "\t")

# Keep only genotyped samples (those present in SNP data)
# This ensures genotype-phenotype matching
myY2 <- data.frame(myY[myY$Taxa %in% colnames(myG2), ])

# ============================================================================
# SECTION 4: MATCH GENOTYPES AND PHENOTYPES
# ============================================================================
# Extract genotypes for samples with phenotype data only
myG3 <- myG2[, colnames(myG2) %in% myY2$Taxa, with = FALSE]

# Reorder phenotypes to match genotype column order
col_order <- c(colnames(myG3))
myY3 <- myY2[match(col_order, myY2$Taxa), ]

# Verify names match between genotypes and phenotypes
all(myY3$Taxa == colnames(myG3))

# Combine SNP metadata with matched genotype calls
myG4 <- data.frame(myG1, myG3)

# Display data structure for verification
str(myY3)
str(myG4)

# Add sample names as first row (required GAPIT format)
myG4[1, ] <- names(myG4)
names(myG4) <- NULL

# ============================================================================
# SECTION 5: SUBSET DATA FOR TESTING (OPTIONAL)
# ============================================================================
# Use smaller dataset for initial testing/validation
# Comment out if running full analysis on all ~449K SNPs
geno_test <- myG4[1:50000, 1:ncol(myG4)]

# ============================================================================
# SECTION 6: RUN GWAS ANALYSIS
# ============================================================================
# Perform GWAS using two complementary models:
#   - FarmCPU: Fast, computationally efficient, good for large datasets
#   - BLINK: Uses penalized regression, controls false positives well
# Parameters:
#   Y: Phenotype data
#   G: Genotype data
#   PCA.total=3: Use first 3 principal components as covariates
#   NJtree.group=3: Use kinship matrix based on 3 groups (population structure)
#   Multiple_analysis=TRUE: Perform analysis with both models
myGAPIT <- GAPIT(
  Y = myY3,
  G = myG4,
  PCA.total = 3,
  model = c("FarmCPU", "Blink"),
  NJtree.group = 3,
  Multiple_analysis = TRUE
)
