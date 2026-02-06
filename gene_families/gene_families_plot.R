
library("tidyverse")
library("ggplot2")

df<- read.csv("all_genes.csv")

# Select and reshape data
gene_families <- c("Histone", "Fatty_acid", "Ethylene", "Cytokinin", 
                   "Auxin", "Squalene", "Self_incompatibility", "LRR")

#genome_order <- c("S.integrifolium", "S.perfoliatum", "H.annuus_ha412", "H.annuus_XRQ","S.taccada")
genome_order <- c("S.taccada", "H.annuus_XRQ", "H.annuus_ha412", "S.perfoliatum", "S.integrifolium")
family_order <- c("Histone", "LRR", "Self_incompatibility",  "Fatty_acid", "Squalene", "Auxin", "Cytokinin", "Ethylene")


heatmap_df <- df %>%
  select(Genome, all_of(gene_families)) %>%
  pivot_longer(cols = -Genome, names_to = "Gene_family", values_to = "Count") %>%
  mutate(
    Genome = factor(Genome, levels = genome_order),
    Gene_family = factor(Gene_family, levels = family_order)
  )

ggplot(heatmap_df, aes(x = Gene_family, y = Genome)) +
  # Draw bubbles
  geom_point(
    aes(size = Count, fill = Gene_family),
    shape = 21, color = "black", alpha = 0.9
  ) +
  # Add count labels on top
  geom_text(aes(label = Count), color = "black", size = 5, fontface = "bold") +
  # Scale bubble size range (adjust to taste)
  scale_size_continuous(range = c(8, 23), guide = "none") +
  # Assign one color per gene family, hide legend
  scale_fill_brewer(palette = "Set3", guide = "none") +
  # Titles & labels
  labs(
    title = "Gene family counts",
    x = "Gene family",
    y = "Genome"
  ) +
  # Tighter layout and larger fonts
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # reduce vertical space between rows
    axis.ticks.length = unit(0, "pt"),
    panel.spacing = unit(0.1, "lines") # squeezes extra white space
  )
