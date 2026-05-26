################################################################################
# Principal Component Analysis (PCA) Visualization
################################################################################
# Purpose: Create PCA plots showing population structure and geographic patterns
#          across Silphium samples; includes both 3D and 2D projections
# Input: pca_data.csv - Principal component scores and sample metadata
# Output: Interactive 3D plot (Plotly), PNG files for PC1 vs PC2 and PC1 vs PC3
#         colored by species groups and geographic coordinates
################################################################################

# ============================================================================
# LOAD REQUIRED LIBRARIES
# ============================================================================
library("ggplot2")      # Grammar of graphics for static plots
library("plotly")       # Interactive 3D visualization
library("scales")       # Scale transformations and color utilities

# ============================================================================
# SECTION 1: DATA INPUT
# ============================================================================
# Load PCA results including PC1, PC2, PC3 scores and metadata
# PC1 explains 8.45% of variance, PC2 explains 2.06%, PC3 explains 1.79%
dat <- read.csv("pca_data.csv", header = TRUE)

# ============================================================================
# SECTION 2: 3D INTERACTIVE VISUALIZATION
# ============================================================================
# Create interactive 3D scatter plot (requires Plotly.js for interactivity)
# Points colored by species/K3 extended population groups
plot_ly(x = dat$PC1, y = dat$PC2, z = dat$PC3,
        type = "scatter3d", mode = "markers",
        color = dat$K3_ext_species)

# ============================================================================
# SECTION 3: 2D PCA PLOT - PC1 vs PC2 COLORED BY SPECIES GROUPS
# ============================================================================
# Create static 2D plot of first two principal components
p1 <- ggplot(dat, aes(x = PC1, y = PC2, color = K3_ext_species, label = row.names(dat))) +
  geom_point(size = 4) +
  # Customize plot appearance
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  # Define custom colors for each species/group
  scale_color_manual(
    values = c(
      "s_integrifolium_east"  = "#EC2D01",
      "s_integrifolium_west" = "#6BAED6",
      "s_integrifolium_northeast"  = "#E7B800",
      "s_integrifolium_x_perfoliatum" = "#9B5FE0",
      "s_radula" = "#000000",
      "s_perfoliatum"  = "#3BB143",
      "s_integrifolium_NA" = "#dcdcdb"
    ),
    name = "Groups"   # Legend title
  )

# Add axis labels with variance explained percentages
p2 <- p1 + xlab("PC1 (8.45%)") + ylab("PC2 (2.06%)")

# Display plot
p2

# Save high-resolution PNG
ggsave("pca_plot_k3_groups_pc1-2.png", plot = p2, height = 5, width = 8, dpi = 300)

# ============================================================================
# SECTION 4: 2D PCA PLOT - PC1 vs PC2 COLORED BY LONGITUDE
# ============================================================================
# Visualize how geographic location (longitude) correlates with genetic structure
p3 <- ggplot(dat, aes(x = PC1, y = PC2)) +
  # Points colored by longitude (east-west geographic position)
  geom_point(aes(color = long), size = 4) +
  # Samples with missing longitude data shown as white circles
  geom_point(
    data = subset(dat, is.na(long)),
    shape = 21, fill = "white", colour = "black", size = 4, stroke = 0.3
  ) +
  # Customize appearance
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  # Use viridis color palette (inferno) for continuous longitude values
  # Reverse direction so east = warm colors (red), west = cool colors (blue)
  scale_colour_viridis_c(begin = 0, end = 1, direction = -1, na.value = NA,
                         option = "inferno", name = "Longitude") +
  xlab("PC1 (8.45%)") + ylab("PC2 (2.06%)")

# Display and save
p3
ggsave("pca_plot_long_pc1-2.png", plot = p3, height = 5, width = 6, dpi = 300)

dev.off()

# ============================================================================
# SECTION 5: 2D PCA PLOT - PC1 vs PC3 COLORED BY LONGITUDE
# ============================================================================
# Visualize geographic patterns on PC1 vs PC3
p4 <- ggplot(dat, aes(x = PC1, y = PC3)) +
  geom_point(aes(color = long), size = 4) +
  # Samples with missing longitude data shown with open symbols
  geom_point(
    data = subset(dat, is.na(long)),
    shape = 21, fill = "white", colour = "black", size = 4, stroke = 0.3
  ) +
  # Customize appearance
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 12)) +
  # Use viridis color palette
  scale_colour_viridis_c(begin = 0, end = 1, direction = -1, na.value = NA,
                         option = "inferno") +
  xlab("PC1 (8.45%)") + ylab("PC3 (1.79%)")

# Display and save
p4
ggsave("pca_plot_long_pc1-3_new.png", plot = p4, height = 5, width = 6, dpi = 300)

dev.off()


