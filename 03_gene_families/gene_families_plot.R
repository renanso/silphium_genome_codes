################################################################################
# Gene Family Visualization Across Multiple Genomes
################################################################################
# Purpose: Creates bubble plot showing counts of selected gene families across
#          five plant species (2 Silphium spp., 2 Helianthus spp., S. taccada)
# Input: all_genes.csv - Gene family counts per genome
# Output: PNG files with publication-quality gene family comparison plots
################################################################################

library("tidyverse")    # Data manipulation and tidyverse functions
library("ggplot2")      # Grammar of graphics for publication-quality plots

# ============================================================================
# DATA INPUT AND SETUP
# ============================================================================
# Load gene count data for all genomes
df <- read.csv("all_genes.csv")

# ============================================================================
# SECTION 1: DATA PREPARATION AND ORDERING
# ============================================================================
# Define gene families of interest for visualization
gene_families <- c("Histone", "Fatty_acid", "Ethylene", "Cytokinin", 
                   "Auxin", "Squalene", "Self_incompatibility", "LRR")

# Define species order (from most distant to most closely related to S. integrifolium)
genome_order <- c("S.taccada", "H.annuus_XRQ", "H.annuus_ha412", 
                  "S.perfoliatum", "S.integrifolium")

# Define gene family display order
family_order <- c("Histone", "LRR", "Self_incompatibility",  
                  "Fatty_acid", "Squalene", "Auxin", "Cytokinin", "Ethylene")

# Transform data to long format for ggplot
# Select only genomes and gene families of interest
heatmap_df <- df %>%
  select(Genome, all_of(gene_families)) %>%
  # Convert wide to long format (Genome | Gene_family | Count)
  pivot_longer(cols = -Genome, names_to = "Gene_family", values_to = "Count") %>%
  # Apply factor ordering to control visualization order
  mutate(
    Genome = factor(Genome, levels = genome_order),
    Gene_family = factor(Gene_family, levels = family_order)
  )

# ============================================================================
# SECTION 2: CREATE BUBBLE PLOT VISUALIZATION
# ============================================================================
# Bubble size represents gene count, color distinguishes gene families
ggplot(heatmap_df, aes(x = Gene_family, y = Genome)) +
  # Draw bubbles with fill color by gene family
  geom_point(
    aes(size = Count, fill = Gene_family),
    shape = 21, color = "black", alpha = 0.9
  ) +
  # Add count labels on top of bubbles
  geom_text(aes(label = Count), color = "black", size = 5, fontface = "bold") +
  # Scale bubble size (adjust range as needed for optimal visualization)
  scale_size_continuous(range = c(8, 23), guide = "none") +
  # Use Set3 color palette (10 distinct colors) and hide legend
  scale_fill_brewer(palette = "Set3", guide = "none") +
  # Labels and title
  labs(
    title = "Gene family counts",
    x = "Gene family",
    y = "Genome"
  ) +
  # Minimal theme with larger fonts for readability
  theme_minimal(base_size = 20) +
  theme(
    # X-axis: rotate labels to prevent overlap
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    # Make axis titles bold
    axis.title = element_text(size = 20, face = "bold"),
    # Center and enlarge plot title
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    # Remove grid lines for cleaner look
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Reduce spacing between rows
    axis.ticks.length = unit(0, "pt"),
    panel.spacing = unit(0.1, "lines")
  )
