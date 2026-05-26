################################################################################
# Manhattan Plot and SNP Density Visualization for GWAS Results
################################################################################
# Purpose: Create Manhattan plots showing genome-wide association scan results
#          with threshold lines for significant associations; highlights SNPs
#          associated with quantitative traits
# Input: snp_data_449K.txt - GWAS results with LOD scores per SNP
# Output: Manhattan plot visualization highlighting significant associations
#         at custom thresholds (CMplot PNG output)
# Reference: CMplot package for publication-quality Manhattan plots
################################################################################

# ============================================================================
# LOAD REQUIRED LIBRARY
# ============================================================================
library("CMplot")   # Publication-quality Manhattan plots for GWAS

# ============================================================================
# SECTION 1: DATA INPUT
# ============================================================================
# Load GWAS results table
# Expected columns: SNP identifiers, chromosome, position, LOD scores
snp_data <- read.table("snp_data_449K.txt", header = TRUE)

# ============================================================================
# SECTION 2: CREATE MANHATTAN PLOT - DENSITIES ACROSS ALL CHROMOSOMES
# ============================================================================
# Generate density Manhattan plot showing LOD scores genome-wide
# Parameters explained:
#   type="h": Histogram/bar plot style
#   plot.type="d": Density plot showing SNP significance across chromosomes
#   LOG10=FALSE: LOD scores already in -log10 scale
#   threshold=3.7: Significance threshold line (genome-wide significant)
#   highlight.type="p": Highlight as point scatter
#   file="jpg": Output format
#   dpi=300: High-resolution output for publication
CMplot(snp_data, type = "h", plot.type = "d", LOG10 = FALSE,
       highlight = SNPs,                               # Highlight specific SNPs (if defined)
       highlight.type = "p",                           # Plot as points
       ylab = "LOD scores",                            # Y-axis label
       threshold = 3.7,                                # Significance threshold
       threshold.col = "black",                        # Color of threshold line
       threshold.lwd = 1,                              # Line width of threshold
       main = "",                                      # Title (empty)
       highlight.col = NULL,                          # Highlight colors (auto-assigned)
       highlight.cex = 1.2,                           # Size of highlighted points
       highlight.pch = 19,                            # Point character (filled circle)
       file = "jpg",                                   # Output file type
       dpi = 300,                                      # Resolution (dots per inch)
       file.output = TRUE,                            # Save to file
       verbose = TRUE,                                 # Print progress to console
       width = 10,                                     # Plot width (inches)
       height = 4,                                     # Plot height (inches)
       band = 1                                        # Band thickness between chromosomes
)

# Display first few rows of data for verification
snp_data[1:2, ]

# ============================================================================
# SECTION 3: ALTERNATIVE MANHATTAN PLOT - MULTI-PANEL BY TRAIT (COMMENTED)
# ============================================================================
# Alternative visualization: Multiple Manhattan plots arranged by phenotype
# This would show results for different traits side-by-side
# Uncomment to use with multi-trait GWAS results:

# CMplot(snp_data, type="h", plot.type="m", LOG10=F,
#        highlight=SNPs, highlight.type="p",
#        ylab="LOD scores",
#        threshold=3.7, threshold.col="black", threshold.lwd=1,
#        main="Seed size QTL",
#        highlight.col=NULL, highlight.cex=1.2, highlight.pch=19,
#        file="jpg", memo="", dpi=300,
#        file.output=TRUE, verbose=TRUE,
#        width=10, height=4, band=1)