################################################################################
# Transposable Element (TE) Abundance and Distribution Analysis
################################################################################
# Purpose: Visualize abundance of transposable element classes across multiple
#          plant genomes; compare TE counts, base pair content, and genome
#          percentage; show TE type distributions within Silphium species
# Input: TE statistics CSV files (count, bp, percentage by class and type)
# Output: Publication-quality bar plots (JPEG format) showing TE content
#         across species and element classes
################################################################################

# ============================================================================
# LOAD REQUIRED LIBRARIES
# ============================================================================
library("tidyverse")    # Data manipulation (dplyr, tidyr, ggplot2)
library("ggplot2")      # Grammar of graphics for visualization
library("scales")       # Scale transformations for plot aesthetics

# ============================================================================
# SECTION 1: DATA INPUT
# ============================================================================
# Load TE statistics from CSV files
# Files contain counts, base pairs, and percentages for each TE class per species
df1 <- read.csv("te_species_count.csv")       # TE sequence counts
df2 <- read.csv("te_species_bp.csv")          # TE base pair content
df3 <- read.csv("te_species_percentage.csv")  # TE percentage of genome
df4 <- read.csv("gene_count.csv")             # Gene counts for comparison

# ============================================================================
# SECTION 2: DATA RESHAPING AND ORDERING
# ============================================================================
# Define species order (phylogenetic/taxonomic arrangement)
order <- c("S.integrifolium", "S.perfoliatum", "H.annuus_HA412HO",
           "H.annuus_XRQ", "Sc.taccada")

# Convert wide format to long format (one row per species-TE class combination)
# This format is required for faceted ggplot visualization

# TE counts - reshape and order species
df_long1 <- df1 %>%
  pivot_longer(cols = starts_with(c("S", "H.")), 
               names_to = "Species", values_to = "Count") %>%
  mutate(Species_ordered = factor(Species, levels = order))

# TE base pair content - reshape and order
df_long2 <- df2 %>%
  pivot_longer(cols = starts_with(c("S", "H.")), 
               names_to = "Species", values_to = "Base_pairs") %>%
  mutate(Species_ordered = factor(Species, levels = order))

# TE genome percentage - reshape and order
df_long3 <- df3 %>%
  pivot_longer(cols = starts_with(c("S", "H.")), 
               names_to = "Species", values_to = "Percentage") %>%
  mutate(Species_ordered = factor(Species, levels = order))

# Gene counts - reshape and order
df_long4 <- df4 %>%
  pivot_longer(cols = starts_with(c("S", "H.")), 
               names_to = "Species", values_to = "Count") %>%
  mutate(Species_ordered = factor(Species, levels = order))
str(df_long4)

# ============================================================================
# SECTION 3: VISUALIZATION 1 - TE SEQUENCE COUNTS
# ============================================================================
# Bar plot showing number of TE sequences per species
# Fill color represents TE class (different classes stacked)
ggplot(df_long1, aes(x = Species_ordered, y = Count, fill = Class)) +
  geom_col() +                                    # Stacked bar plot
  theme_minimal(base_size = 12) +
  scale_fill_brewer(palette = "Set2", name = "TE Class") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", size = 0.3),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, margin = margin(b = 6))
  ) +
  labs(x = NULL, y = "Sequences Count")

ggsave("TE_counts.jpeg", width = 4, height = 4, dpi = 200)

# ============================================================================
# SECTION 4: VISUALIZATION 2 - TE BASE PAIR CONTENT
# ============================================================================
# Bar plot showing total base pairs contributed by each TE class
# Allows comparison of TE impact on genome size across species
ggplot(df_long2, aes(x = Species_ordered, y = Base_pairs, fill = Class)) +
  geom_col() +
  theme_minimal(base_size = 12) +
  scale_fill_brewer(palette = "Set2", name = "TE Class") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", size = 0.3),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, margin = margin(b = 6))
  ) +
  labs(x = NULL, y = "Base pairs (bp)")

ggsave("TE_basepairs.jpeg", width = 4, height = 4, dpi = 200)

# ============================================================================
# SECTION 5: VISUALIZATION 3 - TE GENOME PERCENTAGE
# ============================================================================
# Bar plot showing TE content as percentage of total genome
# Useful for comparing TE burden across species with different genome sizes
ggplot(df_long3, aes(x = Species_ordered, y = Percentage, fill = Class)) +
  geom_col() +
  theme_minimal(base_size = 12) +
  scale_fill_brewer(palette = "Set2", name = "TE Class") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", size = 0.3),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, margin = margin(b = 6))
  ) +
  labs(x = NULL, y = "Percentage (%)")

ggsave("TE_percentage.jpeg", width = 4, height = 4, dpi = 200)

# ============================================================================
# SECTION 6: VISUALIZATION 4 - TE TYPES IN SILPHIUM SPECIES
# ============================================================================
# Detailed breakdown of TE types (not just classes) for Silphium spp.
# Subset to S. integrifolium and S. perfoliatum (columns 2 and 6-7)
df_long3.1 <- df3[, c(1, 2, 6, 7)] %>%
  pivot_longer(cols = starts_with(c("S.", "H.")), 
               names_to = "Species", values_to = "Percentage")

# Stacked bar plot with different TE types
ggplot(df_long3.1, aes(x = Species, y = Percentage, 
                       fill = reorder(Type, Percentage))) +
  geom_col() +
  theme_minimal(base_size = 15) +
  scale_fill_brewer(palette = "Spectral", name = "TE type") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", size = 0.3),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13, color = "black"),
    axis.text.y = element_text(size = 13, color = "black")
       # legend.position = "top",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 9),
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5, margin = margin(b = 6))) +
  labs( x = NULL, y = "Percentage of the genome (%)") + scale_y_continuous( limits = c(0, 100), breaks = seq(0, 100, by = 10))


ggsave("TE_percentage_silphium.jpeg", width = 4, height = 4, dpi = 200)


##gene count


ggplot(df_long4, aes(x = Species_ordered, y = Count, fill = (Species))) +
  geom_col() +
  theme_minimal(base_size = 15) +
  #facet_wrap(~Class) +
  scale_fill_brewer(palette = "Spectral", name = "Genome") +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey80", size = 0.3),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"),
        legend.position = "none",
        #legend.position = element_text(element_blank()),
        #legend.title = element_text(face = "bold"),
        #legend.text = element_text(size = 9),
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5, margin = margin(b = 6))) +
  labs( x = NULL, y = "Gene count") + scale_y_continuous( limits = c(0, 50000), breaks = seq(0, 50000, by = 10000))

ggsave("gene_count.jpeg", width = 4, height = 4, dpi = 200)
