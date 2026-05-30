################################################################################
# Genetic Diversity Analysis: Heterozygosity Estimation
################################################################################
# Purpose: Calculate observed and expected heterozygosity per locus to assess
#          genetic variation and allele frequency patterns across populations
# Input: test.hapmap - Genotype data in HapMap format with SNP metadata
#        info_nj_groups2.csv - Sample population/group assignments
# Output: Heterozygosity statistics, distribution plots, RDS objects for
#         downstream analysis (observed_het1.jpeg, het-distribution1.jpeg)
# References: adegenet and pegas packages for population genetic analysis
################################################################################

# ============================================================================
# LOAD REQUIRED LIBRARIES
# ============================================================================
library("adegenet")     # Genetic data analysis and population genetics
library("hierfstat")    # Hierarchical F-statistics and genetic diversity
library("pegas")        # Population and evolutionary genetics analysis
library("data.table")   # Fast data manipulation

# ============================================================================
# SECTION 1: DATA INPUT AND FORMATTING
# ============================================================================
# Load genotype data in HapMap format
# Format: rs#, alleles, chromosome, position, ... followed by genotype calls
data1 <- read.table("test.hapmap", header = TRUE)

# Create SNP identifiers combining chromosome and position
data1$snp <- paste0(data1$chr, "_", data1$pos)
rownames(data1) <- data1$snp

# Count number of SNPs
n <- ncol(data1) - 1

# ============================================================================
# SECTION 2: TRANSPOSE AND PREPARE GENOTYPE MATRIX
# ============================================================================
# Extract genotype calls (columns 3 through n)
# This matrix has SNPs as columns and samples as rows
data2 <- t(data1[, 3:n])

# Extract sample identifiers from row names
tree_id <- row.names(data2)

# Create data frame with sample IDs and genotypes
data3 <- cbind(tree_id, data2)

# ============================================================================
# SECTION 3: MERGE WITH POPULATION/LOCATION METADATA
# ============================================================================
# Load population assignment information
# Expected columns: sample_id, country, state, and/or subpopulation group
info <- read.csv("info_nj_groups2.csv", header = TRUE)

# Create paired metadata for analysis
info2 <- info[, c(1, 2, 3, 2)]
colnames(info2) <- c("tree_id", "country", "state", "country_state")

# Merge genotype data with population information
data4 <- merge(data3, info2, by = "tree_id")

# ============================================================================
# SECTION 4: REORDER COLUMNS FOR ADEGENET
# ============================================================================
# Calculate number of genotype columns
ncols <- ncol(data4) - 3

# Reorder to put metadata columns first: tree_id, country, state, country_state
# Followed by all SNP columns
new_cols <- c(1, (ncols + 1), (ncols + 2), (ncols + 3), 2:ncols)
data5 <- data4[, new_cols]

# ============================================================================
# SECTION 5: CREATE GENIND OBJECT FOR ADEGENET
# ============================================================================
# Extract SNP genotypes (exclude metadata columns)
locus <- data5[, -c(1, 2, 3, 4)]

# Extract individual labels (sample IDs)
ind <- as.character(data5$tree_id)

# Extract population assignments (geographic state)
population <- as.character(data5$state)

# Create genind object (adegenet format)
# df2genind: dataframe to genind conversion
#   ploidy=2: diploid organism
#   sep="": genotypes are already separated into individual characters
Mydata1 <- df2genind(locus, ploidy = 2, ind.names = ind, 
                     pop = population, sep = "")
Mydata1

# ============================================================================
# SECTION 6: CALCULATE GENETIC DIVERSITY STATISTICS
# ============================================================================
# Count number of alleles per locus
nAll(Mydata1)

# Save genind object for future use
saveRDS(Mydata1, "Mydata1")

# Calculate summary diversity statistics
# This includes observed heterozygosity (Hobs), expected heterozygosity (He),
# and allele frequency distributions
div <- summary(Mydata1)
saveRDS(div, "div1")

# Calculate and display mean observed heterozygosity
mean(div$Hobs)

# Display available summary statistics
names(div)

# ============================================================================
# SECTION 7: VISUALIZE HETEROZYGOSITY
# ============================================================================
# Plot observed heterozygosity per locus
jpeg(file = "observed_het1.jpeg", width = 640, height = 640, units = "px", quality = 100)
plot(div$Hobs, xlab = "Loci number", ylab = "Observed Heterozygosity",
     main = "Observed heterozygosity per locus")
dev.off()

# Plot histogram of heterozygosity distribution across loci
jpeg(file = "het-distribution1.jpeg", width = 640, height = 640, units = "px", quality = 100)
hist(div$Hobs)
dev.off()









