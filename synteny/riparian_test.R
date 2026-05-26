################################################################################
# Genome Synteny and Comparative Genomics Visualization using GENESPACE
################################################################################
# Purpose: Compare genome organization and identify syntenic relationships
#          (conserved gene order) across plant species using GENESPACE pipeline;
#          visualize riparian plots showing synteny relationships between
#          S. integrifolium (reference) and related species
# Input: Complete genome assemblies and gene annotations for 5 plant species
# Output: Synteny block assignments, riparian visualization plots
# Reference: GENESPACE - Genome Synteny Plots Easy Analysis Comparative Exploration
#            (https://github.com/jtlovell/GENESPACE)
################################################################################

# Load GENESPACE library
library(GENESPACE)

# ============================================================================
# SECTION 1: CONFIGURE ANALYSIS PATHS AND TOOLS
# ============================================================================
# NOTE: Modify these paths to match your local system installation

# Working directory where GENESPACE runs and stores intermediate files
wd <- "/cluster/home/rsouza/genespace/riparian_run2"

# Path to MCScanX installation (for synteny block identification)
# MCScanX: Multiple Collinearity Scan toolkit for identifying conserved blocks
path2mcscanx <- "/cluster/home/rsouza/apps/MCScanX-master/"

# Path to OrthoFinder installation (for ortholog identification)
# OrthoFinder: Identifies orthologous genes between species using BLAST
path2orthofinder <- "/cluster/home/rsouza/apps/OrthoFinder/orthofinder"

# Path to DIAMOND (fast BLAST-like alignment for ortholog identification)
path2diamond <- "/cluster/home/rsouza/apps/OrthoFinder/bin/diamond"

# ============================================================================
# SECTION 2: DATA PREPARATION (COMMENTED - ALREADY COMPLETED)
# ============================================================================
# These steps are already completed for this analysis
# Uncomment to repeat if necessary:

# # Create directory for raw genome/annotation files
# dir.create(genomeRepo)

# # Download genome sequences and annotations from NCBI
# rawFiles <- download_exampleData(filepath = genomeRepo)

# # Parse raw FASTA/GFF files to standardized format
# # Aligns headers between sequence and annotation files
# parsedPaths <- parse_annotations(
#   rawGenomeRepo = genomeRepo,
#   genomeDirs = c("carrot", "lettuce", "mimulus", "olives", "sunflower"),
#   genomeIDs = c("carrot", "lettuce", "mimulus", "olives", "sunflower"),
#   presets = "phytozome",         # Use Phytozome genome naming conventions
#   genespaceWd = wd
# )

# ============================================================================
# SECTION 3: INITIALIZE GENESPACE RUN
# ============================================================================
# Set up GENESPACE analysis with tool paths and working directory
# Checks for required files and validates setup
gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx,
  path2orthofinder = path2orthofinder,
  path2diamond = path2diamond
)

# ============================================================================
# SECTION 4: RUN COMPLETE GENESPACE PIPELINE
# ============================================================================
# Execute full analysis pipeline:
#   1. Ortholog identification (OrthoFinder)
#   2. Synteny block identification (MCScanX)
#   3. Block ordering and visualization
out <- run_genespace(gpar)

# ============================================================================
# SECTION 5: CREATE RIPARIAN VISUALIZATION
# ============================================================================
# Define custom ggplot2 theme for riparian plots
ggthemes <- ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"))

# Generate riparian plots showing synteny relationships
# Reference genome: S. integrifolium (Silphium integrifolium)
# Comparison genomes: other Silphium and Helianthus species
# Layout: Chromosomes shown as horizontal bars, colored blocks show syntenic regions
ripDat <- plot_riparian(
  gsParam = out,
  chrFill = "lightgrey",                          # Chromosome background color
  refGenome = "S.integrifolium",                  # Reference species for comparison
  genomeIDs = c("S.taccada", "H.annuus_xrq",    # Species to compare
                 "H.annuus_ha412", "S.perfoliatum", "S.integrifolium"),
  addThemes = ggthemes,                           # Apply custom ggplot theme
  forceRecalcBlocks = FALSE                       # Reuse previously calculated blocks
)
