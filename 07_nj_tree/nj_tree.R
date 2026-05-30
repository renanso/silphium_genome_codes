################################################################################
# Neighbor-Joining Phylogenetic Tree and Population Structure Visualization
################################################################################
# Purpose: Construct neighbor-joining tree from kinship matrix and visualize
#          population structure by grouping samples based on species/population
# Input: kinship_centered_ibs.txt - Identity-by-state kinship matrix
#        info_nj_groups.csv - Population/species group assignments
# Output: Circular phylogram with species groups color-coded (dendrogram.png)
################################################################################

# Clear workspace
rm(list = ls())

# ============================================================================
# LOAD REQUIRED LIBRARIES
# ============================================================================
library("ggtree")       # Grammar of graphics for phylogenetic trees
library("ggplot2")      # Base graphics system used by ggtree
library("dendextend")   # Dendrogram manipulation and visualization

# ============================================================================
# SECTION 1: LOAD AND PREPARE KINSHIP MATRIX
# ============================================================================
# Load kinship matrix (centered identity-by-state distances)
# Rows and columns are sample names, values are pairwise kinship coefficients
theKin <- read.table("kinship_centered_ibs.txt", header = FALSE, row.names = 1)

# Extract sample clone identifiers from row names
clones <- row.names(theKin)

# ============================================================================
# SECTION 2: CONSTRUCT DISTANCE MATRIX AND CLUSTER TREE
# ============================================================================
# Calculate Euclidean distances from kinship matrix
# These distances are used for hierarchical clustering
distance.matrix <- dist(theKin, upper = TRUE)

# Perform hierarchical clustering using complete linkage
hc <- hclust(distance.matrix)

# ============================================================================
# SECTION 3: LOAD POPULATION/SPECIES GROUP INFORMATION
# ============================================================================
# Load metadata defining which group (species/subpopulation) each sample belongs to
# Expected columns: sample ID, species/group designation, geographic info
info <- read.csv("info_nj_groups.csv", header = TRUE)

# ============================================================================
# SECTION 4: DEFINE SAMPLE GROUPINGS
# ============================================================================
# Create a list mapping each group to its constituent sample IDs
# Groups based on K3 extended species classification
group <- list(
  # S. integrifolium subpopulations
  S.integrifolium_E = info[info$K3_ext_species == "s_integrifolium_east", ][, 1],
  S.integrifolium_W = info[info$K3_ext_species == "s_integrifolium_west", ][, 1],
  S.integrifolium_NE = info[info$K3_ext_species == "s_integrifolium_northeast", ][, 1],
  # Hybrid and related species
  S.integrifolium_x_perfoliatum = info[info$K3_ext_species == "s_integrifolium_x_perfoliatum", ][, 1],
  S.radula = info[info$K3_ext_species == "s_radula", ][, 1],
  S.perfoliatum = info[info$K3_ext_species == "s_perfoliatum", ][, 1],
  # Unassigned samples
  S.integrifolium_NA = info[info$K3_ext_species == "s_integrifolium_NA", ][, 1]
)

# ============================================================================
# SECTION 5: CONSTRUCT AND VISUALIZE PHYLOGENETIC TREE
# ============================================================================
# Open PNG file for high-resolution output
tiff("dendrogram.png", width = 30, height = 25, res = 150, units = "cm")

# Create circular phylogram tree from hierarchical clustering
# Size parameter controls line width
p_iris <- ggtree(hc, layout = 'fan', size = 1.1) +
  # Add taxon labels (sample IDs) at tips, rotated to follow tree angle
  geom_tiplab(size = 1.5, aes(angle = angle))

# Assign colors to groups based on K3 species classification
# groupOTU() collapses groups and assigns them colors
p3 <- groupOTU(p_iris, group, 'group') +
  aes(color = group) +
  # Customize legend appearance
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key.size = unit(1, 'cm')) +
  # Define color palette for each group
  # Colors correspond to: East, West, Northeast, Hybrid, Radula, Perfoliatum, NA
  scale_color_manual(values = c("#EC2D01", "#D3D3D3", "#E7B800", "#6BAED6",
                                 "#9B5FE0", "#42F7C0", "#000000"))

# Display the final tree
p3

# Close PNG file
dev.off()