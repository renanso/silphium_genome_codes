################################################################################
# Linkage Disequilibrium Analysis for Candidate SNPs
################################################################################
# Purpose: Calculate pairwise linkage disequilibrium (LD) between candidate SNPs
#          associated with quantitative traits; create correlation matrix and
#          visualization to assess SNP independence in GWAS results
# Input: alleles.csv - Genotype data for candidate SNPs
# Output: LD r² matrix, corrplot visualization (ld_r2_matrix.csv)
# SNPs analyzed: 6 candidate SNPs associated with Ray Floret Count (RFC),
#                Receptacle Diameter (RD), Seed Number Count (SNC),
#                and Agricultural Index (AI), Heat Index (HI)
################################################################################

# Clear workspace
rm(list = ls())

# Load required library for LD calculation
library(SNPRelate)

# ============================================================================
# SECTION 1: DATA INPUT AND GENOTYPE ENCODING
# ============================================================================
# Load genotype data for candidate SNPs
# Format: each column is a SNP, each row is a sample
# Genotypes are stored as two-character strings (e.g., "AA", "AT", "TT")
data <- read.csv("alleles.csv", header = TRUE, stringsAsFactors = FALSE)

# ============================================================================
# SECTION 2: GENOTYPE CONVERSION FUNCTION
# ============================================================================
# Convert genotype strings (diploid format) to numeric codes (0, 1, 2)
# 0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternate
convert_geno <- function(geno, alleles) {
  # Extract individual alleles from genotype string
  a1 <- substr(geno, 1, 1)   # First allele
  a2 <- substr(geno, 2, 2)   # Second allele
  
  # Define reference and alternate alleles
  ref <- alleles[1]
  alt <- alleles[2]
  
  # Return numeric code based on allele content
  if (a1 == ref && a2 == ref) return(0)        # Homozygous reference
  if (a1 == alt && a2 == alt) return(2)        # Homozygous alternate
  return(1)                                     # Heterozygous
}

# ============================================================================
# SECTION 3: EXTRACT ALLELES AND ENCODE GENOTYPES
# ============================================================================
# Identify all unique alleles at each locus
# This ensures correct reference/alternate assignment
alleles_list <- lapply(data, function(x) sort(unique(unlist(strsplit(x, "")))))

# Initialize genotype matrix (rows: samples, columns: SNPs)
geno_matrix <- matrix(0, nrow = nrow(data), ncol = ncol(data))

# Convert each genotype string to numeric code
for (i in 1:ncol(data)) {
  alleles <- alleles_list[[i]]
  geno_matrix[, i] <- sapply(data[, i], convert_geno, alleles = alleles)
}

# ============================================================================
# SECTION 4: CREATE GDS FILE FOR SNPRelate
# ============================================================================
# SNPRelate requires data in GDS (Genomic Data Structure) format
# SNPs are rows, samples are columns in GDS format (hence the transpose)
snpgdsCreateGeno("alleles.gds",
  genmat = t(geno_matrix),                                    # Transpose for correct format
  sample.id = paste0("sample", 1:nrow(data)),                # Sample identifiers
  snp.id = colnames(data),                                   # SNP identifiers
  snp.position = 1:ncol(data),                               # Arbitrary position (no physical map)
  snp.chromosome = rep(1, ncol(data)),                       # All SNPs assigned to chromosome 1
  snp.allele = sapply(alleles_list, function(a) paste(a, collapse = "/"))  # Allele pairs
)

# ============================================================================
# SECTION 5: CALCULATE LINKAGE DISEQUILIBRIUM
# ============================================================================
# Open the GDS file for analysis
genofile <- snpgdsOpen("alleles.gds")

# Calculate LD matrix using Pearson correlation
# Parameters:
#   slide = -1: Compare all SNP pairs (no sliding window)
#   method = "r": Use Pearson correlation coefficient
ld <- snpgdsLDMat(genofile, slide = -1, method = "r")

# Extract r² values (correlation coefficient squared = r²)
# r² values range from 0 (no LD) to 1 (complete LD)
r2_matrix <- ld$LD^2

# Assign SNP names to rows and columns
rownames(r2_matrix) <- colnames(data)
colnames(r2_matrix) <- colnames(data)

# Display r² matrix
r2_matrix

# ============================================================================
# SECTION 6: CUSTOMIZE SNP LABELS FOR VISUALIZATION
# ============================================================================
# Create informative row labels including chromosome position and associated traits
row.names(r2_matrix) <- c("Chr01_695173749\nRFC/SNC", 
                          "Chr01_1583902655\nRFC/RD/SNC",
                          "Chr03_721440854\nRFC/RD",
                          "Chr04_222292225\nRFC", 
                          "Chr06_274305900\nAI",
                          "Chr03_78530134\nHI")

# Create informative column labels (with extra spacing for readability)
colnames(r2_matrix) <- c("Chr01_695173749\n  RFC/SNC",
                         "Chr01_1583902655\n  RFC/RD/SNC",
                         "Chr03_721440854\n  RFC/RD",
                         "Chr04_222292225\n RFC", 
                         "Chr06_274305900\n  AI",
                         "Chr03_78530134\n  HI")

# ============================================================================
# SECTION 7: VISUALIZE LINKAGE DISEQUILIBRIUM
# ============================================================================
# Load visualization library
library(corrplot)

# Define color palette: red (negative/low values) → white (0) → blue (positive/high values)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

# Create correlation plot
# method="color": Use color intensity to represent r² values
# type="lower": Show only lower triangle (symmetric matrix)
# diag=FALSE: Hide diagonal (r² = 1 for self-comparison)
corrplot(r2_matrix, method = "color", col = col(200),
         type = "lower", order = "original",
         # Add r² values as text on the plot
         addCoef.col = "black", number.cex = 0.8,
         # Format text labels
         tl.col = "black", tl.srt = 45, tl.cex = 0.8,
         # Hide diagonal values (all 1.0)
         diag = FALSE
)

# ============================================================================
# SECTION 8: SAVE RESULTS
# ============================================================================
# Export r² matrix as CSV for downstream analysis
write.csv(r2_matrix, "ld_r2_matrix.csv")


