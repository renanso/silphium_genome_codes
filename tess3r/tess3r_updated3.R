################################################################################
# Population Structure Analysis using TESS3R
################################################################################
# Purpose: Perform spatially-informed ancestry estimation (TESS3) to identify
#          population substructure and admixture patterns; combines genetic and
#          geographic information to assign samples to ancestral populations
# Input: Numeric genotype data (HapMap format), geographic coordinates
# Output: Q-matrix (ancestry proportions), cross-validation plots, structure plots
# Reference: TESS3R package for large-scale population genomics
#            Caye et al. (2016) PNAS 113(48):13603-13608
################################################################################

# Clear workspace
rm(list = ls())

# ============================================================================
# LOAD REQUIRED LIBRARIES
# ============================================================================
library("Rcpp")            # C++ code integration for performance
library("tess3r")          # Spatially-informed ancestry estimation
library("maps")            # Map drawing utilities
library("BiocParallel")    # Parallel computation
library("data.table")      # Fast data frame operations
library("ggplot2")         # Grammar of graphics for visualization
library("ggrepel")         # Better label placement in plots

# Test parallel processing with 8 cores
bplapply(1:10, print, BPPARAM = MulticoreParam(workers = 8))

# ============================================================================
# SECTION 1: LOAD AND PREPARE GENOTYPE DATA
# ============================================================================
# Load numeric genotype data (0/1/2 coding for diploid calls)
# File format: rows are SNPs, columns are samples
# Expected format: first column = SNP identifier
data2 <- fread(file = "test.txt", sep = "\t")
colnames(data2)[1] <- "clone"

# Alternative: Load data from HapMap format if different file available
# data2 <- read.table("RSouza_silphium.Ihapmap.150gps.numeric2.hmp.txt", header = T)
# colnames(data2)[1] <- "clone"

# ============================================================================
# SECTION 2: LOAD GEOGRAPHIC COORDINATES
# ============================================================================
# Load latitude/longitude coordinates for each sample
# Format: sample ID, latitude, longitude columns
coord2 <- read.table("coordinates_integrifolium.txt", header = T)

# Alternative coordinate file (commented):
# coord2 <- read.table("coordinates2.txt", header = T)

# ============================================================================
# SECTION 3: MERGE GENOTYPES WITH COORDINATES
# ============================================================================
# Combine genotype and coordinate data, matching by clone ID
n1 <- ncol(data2) + 1
n2 <- ncol(data2) + 2

# Extract coordinates only for samples with genotypes
coord3 <- (merge(data2, coord2, by = 'clone'))[, n2:n1]
coord3 <- as.matrix(coord3)

# Save coordinate matrix for downstream use
write.table(coord3, file = "coordinates_matrix.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE)

# Extract genotype data matched with coordinates
data3 <- (merge(data2, coord2, by = 'clone'))[, 1:ncol(data2)]

# ============================================================================
# SECTION 4: VISUALIZE SAMPLE COLLECTION LOCATIONS
# ============================================================================
# Generate map showing geographic distribution of samples
png(filename = "plot.png", width = 10000, height = 10000, res = 300, units = "px")

# Convert coordinates to data frame for ggplot
df <- as.data.frame(coord3)
df$label <- labels

# Create map with sample locations and labels
ggplot(coord2, aes(x = long, y = lat)) +
  geom_point(size = 1.5) +
  # Add state boundaries as context
  borders("state") +
  # Add repelling labels to avoid overlaps
  geom_text_repel(aes(label = clone), size = 3, color = "blue",
                  max.overlaps = Inf) +
  labs(x = "Longitude (°E)", y = "Latitude (°N)") +
  # Limit map extent to study region
  coord_quickmap(
    xlim = c(-100, -85),
    ylim = c(30, 45),
    expand = FALSE
  ) +
  theme_minimal(base_size = 14)
dev.off()

# ============================================================================
# SECTION 5: PREPARE DATA FOR TESS3R ANALYSIS
# ============================================================================
# Set row names to sample IDs for downstream analysis
rownames(data3) <- data3$clone
# Remove redundant clone column
data3$clone <- NULL

# ============================================================================
# SECTION 6: RUN TESS3 ANALYSIS
# ============================================================================
# Perform TESS3 with K ranging from 1 to 8 ancestral populations
# This explores different levels of population structure
# Parameters:
#   X: Genotype matrix (rows = samples, columns = SNPs)
#   coord: Geographic coordinates for spatial priors
#   K: Number of ancestral populations to test (1:8)
#   method: "projected.ls" = Least squares with spatial constraints
#   ploidy: 2 = diploid organism
#   openMP.core.num: Number of cores for parallel processing
tess3.obj <- tess3(X = data3, coord = coord3, K = 1:8,
                   method = "projected.ls", ploidy = 2, openMP.core.num = 8)

# Save TESS3 object for future reference
saveRDS(tess3.obj, file = "tess3.rds")

# ============================================================================
# SECTION 7: MODEL SELECTION - CROSS-VALIDATION
# ============================================================================
# Plot cross-validation scores to select optimal K
# Lower score = better fit of the model
png(filename = "plot1.png", width = 480, height = 480, units = "px")
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")
dev.off()

# ============================================================================
# SECTION 8: EXTRACT Q-MATRICES FOR SELECTED K VALUES
# ============================================================================
# Q-matrix: ancestry proportions for each individual in each population
# Rows = samples, Columns = K populations, Values = proportion of ancestry

# Extract Q-matrices for K = 2, 3, 4, 5
# These represent different hypotheses about population number
q.matrix2 <- qmatrix(tess3.obj, K = 2)
q.matrix3 <- qmatrix(tess3.obj, K = 3)
q.matrix4 <- qmatrix(tess3.obj, K = 4)
q.matrix5 <- qmatrix(tess3.obj, K = 5)

# Save Q-matrices for downstream visualization
saveRDS(q.matrix2, file = "q.matrix2.rds")
saveRDS(q.matrix3, file = "q.matrix3.rds")
saveRDS(q.matrix4, file = "q.matrix4.rds")
saveRDS(q.matrix5, file = "q.matrix5.rds")

# ============================================================================
# SECTION 9: DEFINE COLOR PALETTE FOR Q-PLOTS
# ============================================================================
# Create color palette for STRUCTURE-like bar plots
# One color per ancestral population
my.colors <- c("tomato", "orange", "lightblue", "wheat", "olivedrab")

# Create palette with specified colors
my.palette <- CreatePalette(my.colors, 5)

# Note: Use bar.plot() function from tess3r to create STRUCTURE-like plots
# Example:
# bar.plot(q.matrix2, coord2, colorPalette = my.palette,
#          map.polygon = TRUE, cex = 0.4)

#library(RColorBrewer)
#my.colors <- brewer.pal(n = 5, name = "Dark2")
#my.palette <- CreatePalette(my.colors)

png(filename="plot2.png", width = 1800, height = 1000, units = "px")
par(mfrow = c(2, 2))
barplot(q.matrix2, border = NA, space = 0, sort.by.Q =TRUE,
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "K=2", 
        col.palette = my.palette) -> bp
axis(1, at = 1:nrow(q.matrix2), labels = bp$order, las = 3, cex.axis = 1) 

barplot(q.matrix3, border = NA, space = 0, sort.by.Q =TRUE,
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "K=3", 
        col.palette = my.palette) -> bp
axis(1, at = 1:nrow(q.matrix3), labels = bp$order, las = 3, cex.axis = 1)

barplot(q.matrix4, border = NA, space = 0, sort.by.Q =TRUE,
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "K=4", 
        col.palette = my.palette) -> bp
axis(1, at = 1:nrow(q.matrix4), labels = bp$order, las = 3, cex.axis = 1)

barplot(q.matrix5, border = NA, space = 0, sort.by.Q =TRUE,
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "K=5", 
        col.palette = my.palette) -> bp
axis(1, at = 1:nrow(q.matrix5), labels = bp$order, las = 3, cex.axis = 1)
dev.off()

###
#png(filename="plot3.png", width = 1000, height = 1300, units = "px")
#par(mfrow = c(2, 2))

#plot(q.matrix2, coord3, method = "map.max", interpol = FieldsKrigModel(10),  
#     main = "Ancestry model K=2",
#     xlab = "Longitude", ylab = "Latitude", 
#     resolution = c(300,300), cex = .5,
#     col.palette = my.palette)
#map('state', add = T, interior =T)
#map('lakes', add=TRUE, fill=TRUE, col='white', boundary='black')

#plot(q.matrix3, coord3, method = "map.max", interpol = FieldsKrigModel(10),  
#     main = "Ancestry model K=3",
#     xlab = "Longitude", ylab = "Latitude", 
#     resolution = c(300,300), cex = .5,
#     col.palette = my.palette)
#map('state', add = T, interior =T)
#map('lakes', add=TRUE, fill=TRUE, col='white', boundary='black')

#plot(q.matrix4, coord3, method = "map.max", interpol = FieldsKrigModel(10),  
#     main = "Ancestry model K=4",
#     xlab = "Longitude", ylab = "Latitude", 
#     resolution = c(300,300), cex = .5,
#     col.palette = my.palette)
#map('state', add = T, interior =T)
#map('lakes', add=TRUE, fill=TRUE, col='white', boundary='black')

#plot(q.matrix5, coord3, method = "map.max", interpol = FieldsKrigModel(10),  
#     main = "Ancestry model K=5",
#     xlab = "Longitude", ylab = "Latitude", 
#     resolution = c(300,300), cex = .5,
#     col.palette = my.palette)
#map('state', add = T, interior =T)
#map('lakes', add=TRUE, fill=TRUE, col='white', boundary='black')

#dev.off()

# retrieve tess3 results for different K
p.values2 <- pvalue(tess3.obj, K = 2)
p.values3 <- pvalue(tess3.obj, K = 3)
p.values4 <- pvalue(tess3.obj, K = 4)
p.values5 <- pvalue(tess3.obj, K = 5)

saveRDS(p.values2, file = "p.values2.rds")
saveRDS(p.values3, file = "p.values3.rds")
saveRDS(p.values4, file = "p.values4.rds")
saveRDS(p.values5, file = "p.values5.rds")

png(filename="plot4.png", width = 1000, height = 700, units = "px")
par(mfrow = c(2, 2))
hist(p.values2, col = "lightblue", main = "K=2") 
hist(p.values3, col = "lightblue", main = "K=3") 
hist(p.values4, col = "lightblue", main = "K=4") 
hist(p.values5, col = "lightblue", main = "K=5") 
dev.off()

##Fst calculations
# Benjamini-Hochberg algorithm
fdr.level = 1e-4
L = length(p.values2)

w2 = which(sort(p.values2) < fdr.level * (1:L)/L)
w3 = which(sort(p.values3) < fdr.level * (1:L)/L)
w4 = which(sort(p.values4) < fdr.level * (1:L)/L)
w5 = which(sort(p.values5) < fdr.level * (1:L)/L)

candidates2 = order(p.values2)[w2]
candidates3 = order(p.values3)[w3]
candidates4 = order(p.values4)[w4]
candidates5 = order(p.values5)[w5]

length(candidates2)
length(candidates3)
length(candidates4)
length(candidates5)

# Manhattan plot 
png(filename="plot5.png", width = 1000, height = 1300, units = "px")
par(mfrow = c(2, 2))

plot(p.values2, main = "Fst plot K=2", 
     xlab = "Locus id", ylab = "-log10(P-values)",
     cex = .3, col = "grey")
points(candidates2, -log10(p.values2)[candidates2], 
     pch = 19, cex = .2, col = "blue")

plot(p.values3, main = "Fst plot K=3", 
     xlab = "Locus id", ylab = "-log10(P-values)",
     cex = .3, col = "grey")
points(candidates3, -log10(p.values3)[candidates3], 
       pch = 19, cex = .2, col = "blue")

plot(p.values4, main = "Fst plot K=4", 
     xlab = "Locus id", ylab = "-log10(P-values)",
     cex = .3, col = "grey")
points(candidates4, -log10(p.values4)[candidates4], 
       pch = 19, cex = .2, col = "blue")

plot(p.values5, main = "Fst plot K=5", 
     xlab = "Locus id", ylab = "-log10(P-values)",
     cex = .3, col = "grey")
points(candidates5, -log10(p.values5)[candidates5], 
       pch = 19, cex = .2, col = "blue")
dev.off()