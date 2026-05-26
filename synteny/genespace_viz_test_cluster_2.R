################################################################################
# Genome Annotation Classification and Visualization
################################################################################
# Purpose: Classify genome regions based on feature content (genes, TEs, etc);
#          calculate sliding window statistics to visualize genomic landscape;
#          create stacked area plots showing feature distribution across
#          chromosomes
# Input: Genome assembly FASTA, gene GFF3 annotations, repeat element BED files
# Output: Classified feature windows, visualization plots
# Features analyzed: Genes, Helitrons, TIR elements, LTR elements
################################################################################

library(devtools)
# devtools::install_github("jtlovell/GENESPACE", upgrade = F)
library(GENESPACE)

# Data manipulation and handling
library(data.table)     # Fast data frame operations
library(Biostrings)     # Biological sequences (DNA, RNA, proteins)
library(ggplot2)        # Grammar of graphics for visualization

# Genomic range operations
library(IRanges)        # Integer ranges and range operations
library(GenomicRanges)  # Genomic ranges (chromosomes + coordinates)
library(GenomeInfoDb)   # Genome information and metadata
library(rtracklayer)    # Read/write genome track files (GFF, BED, etc)

# ============================================================================
# SECTION 1: DEFINE PATHS AND SETUP
# ============================================================================
# Working directory for GENESPACE analysis
wd <- "/cluster/home/rsouza/genespace/"
if (!dir.exists(wd))
  dir.create(wd)

# Path to S. integrifolium genome assembly
faFile <- file.path(wd, "genomes/Silphium_integrifolium_var_Bad_Astra.mainGenome.fasta")

# Path to gene annotations (GFF3 format)
geneGffFile <- file.path(wd, "genes/sintegrifolium.gene.gff3")

# Path to transposable element annotations (BED format, EDTA pipeline output)
repBedFile <- file.path(wd, "repeats/all.int.fa.mod.EDTA.RM.bed")

# ============================================================================
# SECTION 2: LOAD AND PROCESS GENOME SEQUENCE
# ============================================================================
# Read genome assembly as DNA string set
dnaSS <- Biostrings::readDNAStringSet(faFile)

# Simplify sequence names (remove space-separated annotation)
names(dnaSS) <- sapply(names(dnaSS), function(x) strsplit(x, " ")[[1]][1])

# Filter to major chromosomes/scaffolds (>1 million bp)
dnaSS <- dnaSS[Biostrings::width(dnaSS) > 1e6]

# Extract sequence information for GenomicRanges object
seqInfo <- pull_seqInfo(dnaSS)

# ============================================================================
# SECTION 3: LOAD GENE ANNOTATIONS
# ============================================================================
# Read gene features from GFF3 file
genes <- as.data.frame(rtracklayer::readGFF(
  geneGffFile,
  columns = c("seqid", "start", "end", "type"),   # Required columns
  tags = "Parent",                                  # Extract parent ID tag
  filter = list(type = c("gene"))                  # Keep only gene features
))

# Subset to genes on major chromosomes
genes <- subset(genes, seqid %in% names(dnaSS))

# ============================================================================
# SECTION 4: LOAD TRANSPOSABLE ELEMENT ANNOTATIONS
# ============================================================================
# Read repeat elements from BED file (from EDTA pipeline)
# Columns: chromosome, start, end, element_type
repeats <- fread(
  repBedFile,
  select = c(1:3, 12),                            # Chr, start, end, type
  col.names = c("chr", "start", "end", "type")
)

# Subset to repeats on major chromosomes
repeats <- subset(repeats, chr %in% names(dnaSS))

# ============================================================================
# SECTION 5: ORGANIZE FEATURES INTO BED LIST
# ============================================================================
# Create separate BED dataframes for each feature type
# These will be used for classification and visualization
bedList <- list(
  # Protein-coding genes
  gene = with(subset(genes, type == "gene"), data.frame(
    chr = seqid, start = start, end = end)),
  
  # Transposable element classes
  # Helitrons: Rolling-circle transposons
  Helitron = with(subset(repeats, grepl("Helitron", type)), data.frame(
    chr = chr, start = start, end = end)),
  
  # TIR elements: DNA transposons with Terminal Inverted Repeats
  TIR = with(subset(repeats, grepl("TIR", type)), data.frame(
    chr = chr, start = start, end = end)),
  
  # LTR elements: Retrotransposons with Long Terminal Repeats
  LTR = with(subset(repeats, grepl("LTR", type)), data.frame(
    chr = chr, start = start, end = end))
)

# Alternative: include other repeat types
# otherRepeat = with(subset(repeats, !grepl("Helitron|TIR|LTR", type)), data.frame(
#   chr = chr, start = start, end = end))

# Display structure of feature list
print(lapply(bedList, head))

# Examine repeat element distribution
summary(repeats)

# ============================================================================
# SECTION 6: CLASSIFY GENOME INTO FEATURE CATEGORIES
# ============================================================================
# Use GENESPACE function to classify each genomic position
suppressWarnings(genomeClasses <- classify_genome(
  dnaSS = dnaSS,                      # Genome sequences
  listOfBeds = bedList,               # Feature annotations
  verbose = TRUE                       # Print progress
))

# Display classified regions
print(lapply(genomeClasses, head))

# ============================================================================
# SECTION 7: CALCULATE SLIDING WINDOW STATISTICS
# ============================================================================
# Compute feature proportions in sliding windows across chromosomes
# Window parameters: 5 Mb windows, 1 Mb step size
sws <- slide_genome(
  seqInfo = seqInfo,                  # Chromosome information
  listOfGrs = genomeClasses[c(1, 5, 2, 3, 4)],  # Features to analyze
  windowSize = 5e6,                   # 5 Mb window size
  stepSize = 1e6                      # 1 Mb sliding step
)
head(sws)

# ============================================================================
# SECTION 8: CREATE STACKED AREA VISUALIZATION
# ============================================================================
# Define colors for each feature type
plotCols <- c("#CC6828",              # Gene: tan/brown
              "#FFFFFF",              # Unknown: white
              "#0F4F8B",              # Helitron: dark blue
              "#4C86C6",              # TIR: light blue
              "#AED8E6")              # LTR: pale blue

# Define custom ggplot2 theme for clean genomic visualization
plotTheme <- theme(
  panel.background = element_rect(fill = "black"),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  panel.spacing = unit(0, "cm"),
  strip.background = element_blank()
)

# Order feature types for consistent stacking
sws$id2 <- factor(sws$id, levels = c("gene", "unknown", "Helitron", "TIR", "LTR"))

# Create stacked area plot showing feature distribution across genome
ggplot(sws, aes(x = start/1e6, y = propWind, fill = id2)) +
  geom_area() +                                    # Stacked area plot
  scale_fill_manual(values = plotCols) +          # Custom colors
  facet_grid(. ~ gsub("chr", "", chr), scale = "free", space = "free") +  # Chromosome facets
  scale_x_continuous(expand = c(0.05, 0.05)) +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  plotTheme +
  labs(x = "Chromosomes (scaled by physical size; 5Mb windows, 1Mb steps)",
       y = "Cumulative % of Sequence") +
  theme(axis.text.x = element_text(face = "bold", size = 10))
        axis.text.y = element_text(face="bold",
                                   size=10))

ggsave("plot_new2.png", dpi = 150, width = 15, height = 6)

#p2 <- ggplot(sws,
#            aes(x = start/1e6, y = propWind, fill = id))+
#  geom_area()+
#  scale_fill_manual(values = plotCols) +
#  facet_grid(.~gsub("chr","",chr), scale = "free", space = "free") +
#  scale_x_continuous(expand = c(0.05,0.05))+
#  scale_y_continuous(expand = c(0.005,0.005))+
#  plotTheme +
#  labs(x = "Chromosomes (scaled by physical size; 5Mb windows, 1Mb steps)",
#       y = "Cumulative % of Sequence")

#pdf(file = "/cluster/home/rsouza/genespace/Plot_2.pdf",   # The directory you want to save the file in
#    width = 4, # The width of the plot in inches
#    height = 4)

#gridExtra::grid.arrange(p1, nrow = 1)

#dev.off()