# Silphium Genome Analysis

Scripts for analysis presented in the manuscript: <i>"Assembly of Silphium interspecific hybrid genomes opens the genus to phylogenomics, ecogenomics, and molecular breeding."</i>

## Overview

This repository contains reproducible analysis pipelines for genome assembly quality control, population genomics, phenotypic data analysis, comparative genomics, and gene family characterization.

## Getting Started

All code for data analysis and figure generation is organized in the folders described below. Source data is included to regenerate results. The analyses were performed on:
- **R version**: 4.4.2 (package dependencies listed in individual scripts)
- **Operating System**: x86_64-redhat-linux-gnu
- **Bash version**: GNU bash 5.1.8(1)-release
- **Job Scheduler**: SLURM v24.11.6 (for computationally intensive analyses)

### Prerequisites
* R (v4.4.2) with required packages specified in each script
* x86_64-redhat-linux-gnu compatible system
* SLURM job scheduler (for HPC cluster submissions)

---

## Directory Structure

### `gea_data/` — Geographic and Climate Data Visualization
**Purpose**: Visualize sample collection locations and environmental context

**Main Script**: `climate_data.R`
- **Workflow**: Integrates geographic coordinates with WorldClim v2.1 climate variables
- **Outputs**:
  - Geographic maps showing Silphium sample collection sites
  - Climate data layers (temperature, precipitation, solar radiation, elevation) for field sites
- **Key Features**:
  - High-resolution climate rasters at 10-minute resolution

---

### `gene_families/` — Gene Family Abundance Analysis
**Purpose**: Compare gene family copy numbers across plant genomes

**Main Script**: `gene_families_plot.R`
- **Workflow**: Analyzes abundance of selected gene families across 5 plant genomes (2 Silphium spp., 2 Helianthus spp., S. taccada)
- **Gene Families Analyzed**:
  - Histone: Chromatin-associated proteins
  - Fatty Acid: Lipid metabolism
  - Ethylene: Hormone signaling
  - Cytokinin: Growth regulators
  - Auxin: Plant development hormone
  - Squalene: Sterol synthesis
  - Self-incompatibility: Reproductive isolation
  - LRR: Disease resistance (Leucine-Rich Repeat)
- **Outputs**: Plots showing gene family copy numbers with counts overlaid

---

### `genome_scans/` — Genome-Wide Association Study (GWAS) and Genome-environment analysis (GEA)
**Purpose**: Identify SNPs significantly associated with quantitative traits

**Main Script**: `gwas_gapit.R`
- **Workflow**: Performs genome-wide association analysis using GAPIT3 package
- **Input Data**:
  - 449K SNPs genotyped across 258 Silphium samples
  - Phenotypes
- **Analysis Methods**:
  - FarmCPU model: Fast, computationally efficient mixed model
  - BLINK model: Penalized regression with false positive control
- **Quality Control**:
  - Biallelic SNPs only
  - Population structure correction via kinship matrix (3 groups)
- **Outputs**: Manhattan plots, Q-Q plots, SNP significance lists with p-values and effect estimates

---

### `heritability_correlations/` — Phenotypic Analysis and G×E Interactions
**Purpose**: Quantify trait heritability and genotype-by-environment interactions

**Scripts**:
- `cor_var_her.R`: Phenotypic correlation and variance component analysis
- `gxe.R`: Genotype-by-environment interaction visualization

**Workflow**:
1. **Phenotypic Data Processing** (`cor_var_her.R`):
   - Load phenotypic data with proper factor/numeric type conversion
   - Subset data by environment (Alabama 2024, Kansas 2024, Alabama 2025, Kansas 2025)
   - Calculate Pearson correlations between traits
   - Perform PCA-based imputation for missing values
   - Statistical significance testing of correlations

2. **G×E Interaction Analysis** (`gxe.R`):
   - Fit linear mixed effects models: `Trait ~ Clone + Rep + Env + Interactions`
   - Extract estimated marginal means (emmeans) for visualization
   - Create interaction plots showing clone-specific responses across environments
   - Overlay population means for reference

**Traits Analyzed**: RFC (Ray Floret Count), SM (Seed Mass), SNC (Seed Number per Capitulum), RD (Receptacle Diameter)

**Outputs**: 
- Correlation matrices and significance tables
- Interaction plots with clone trajectories
- ANOVA tables partitioning variance

---

### `nj_tree/` — Phylogenetic Tree and Population Structure
**Purpose**: Visualize genetic relationships and population substructure

**Main Script**: `nj_tree.R`
- **Workflow**: Constructs neighbor-joining phylogenetic tree from kinship matrix and colors clades by species/subpopulation
- **Input Data**:
  - Identity-by-state (IBS) kinship matrix calculated from SNP data
  - Sample population/species assignments (K3 extended classification)
- **Population Groups Visualized**:
  - S. integrifolium east/west/northeast subpopulations
  - S. integrifolium × S. perfoliatum hybrids
  - S. radula
  - S. perfoliatum
  - Unassigned samples
- **Output**: Circular phylogram with color-coded population groups (dendrogram.png)
- **Interpretation**: Reveals genetic structure
---

### `pca/` — Principal Component Analysis (PCA) Visualization
**Purpose**: Explore population structure and geographic genetic patterns

**Main Script**: `pca_plot.R`
- **Workflow**: Creates 2D and 3D visualizations of genetic principal components
- **Analysis Dimensions**:
  - PC1 vs PC2 colored by species/population groups
  - PC1 vs PC3 colored by geographic location (longitude)
- **Variance Explained**: PC1=8.45%, PC2=2.06%, PC3=1.79%
- **Outputs**:
  - Interactive 3D plot (Plotly format)
  - Static 2D plots colored by species assignment
  - Geographic gradient plots showing east-west population structure
- **Insights**: Identifies geographic structure and admixture patterns in natural populations

---

### `phenotypes/` — Mixed Effects Phenotypic Modeling
**Purpose**: Estimate breeding values and partition phenotypic variance

**Main Script**: `phenos_analysis.R`
- **Workflow**: 
  1. Data cleaning and exploratory visualization
  2. Subsetting by environment and G×E interaction factor creation
  3. Calculation of Best Linear Unbiased Estimates (BLUEs) for each clone
- **Traits Modeled**:
  - Ray floret count
  - Average seed mass
  - Seed number per head
  - Receptacle diameter
- **Outputs**: 
  - Summary statistics per trait
  - Random effect estimates (clone intercepts + BLUPs)
  - Fixed effect estimates (overall means)
  - Full model summaries with variance explained

---

### `snp_stats/` — SNP Statistics and Genetic Diversity
**Purpose**: Analyze SNP patterns, linkage disequilibrium, and heterozygosity

**Scripts**:
- `alleles_ld.R`: Linkage disequilibrium between candidate SNPs
- `heterozygosity1.R`: SNP heterozygosity estimation
- `SNP_density.R`: SNP density plot visualization

**Workflow**:

1. **Linkage Disequilibrium Analysis** (`alleles_ld.R`):
   - Convert genotype strings to numeric encoding
   - Calculate pairwise correlation (r²) between 6 candidate SNPs
   - Output: LD matrix heatmap with correlation coefficients

2. **Heterozygosity Analysis** (`heterozygosity1.R`):
   - Load HapMap format genotype data
   - Convert to genind format (adegenet package)
   - Calculate observed and expected heterozygosity per locus
   - Analyze heterozygosity distribution across loci
   - Output: Summary statistics, distribution plots per population

3. **Manhattan Plot Visualization** (`SNP_density.R`):
   - Show SNP density across chromosomes

---

### `synteny/` — Comparative Genomics and Synteny Analysis
**Purpose**: Compare genome organization across plant species

**Scripts**:
- `riparian_test.R`: GENESPACE pipeline for synteny identification
- `genespace_viz_test_cluster_2.R`: Genome feature classification and visualization

**Workflow**:

1. **Synteny Analysis** (`riparian_test.R`):
   - Run complete GENESPACE pipeline: ortholog identification → synteny block detection → visualization
   - Identify conserved syntenic blocks and genome rearrangements
   - Output: Riparian plots showing synteny relationships

2. **Genome Classification** (`genespace_viz_test_cluster_2.R`):
   - Classify genomic regions by feature type (genes, TIRs, Helitrons, LTRs)
   - Calculate sliding window statistics (5 Mb windows, 1 Mb steps)
   - Visualize feature distribution across chromosomes
   - Output: Stacked area plots showing cumulative % of sequence types per chromosome

---

### `te_families/` — Transposable Element Analysis
**Purpose**: Characterize TE abundance and distribution across genomes

**Main Script**: `plot_te.R`
- **Workflow**: Compares TE content across 5 plant genomes using EDTA pipeline results
- **TE Classes Analyzed**:
  - Class I (Retrotransposons): LTR-RTs, LINEs
  - Class II (DNA Transposons): TIR elements, Helitrons
- **Comparison Metrics**:
  - Total TE sequence counts per species
  - TE base pair content (bp)
  - TE genome percentage (identifies TE-rich genomes)
  - TE type diversity (detailed breakdown for Silphium)
---

### `tess3r/` — Spatially-Informed Population Structure
**Purpose**: Identify population substructure incorporating geographic information

**Main Script**: `tess3r_updated3.R`
- **Workflow**: 
  1. Load numeric genotype data and sample coordinates
  2. Perform TESS3 analysis testing K=1-8 ancestral populations
  3. Use spatial priors to constrain ancestry assignments to geographic patterns
  4. Model selection via cross-validation
  5. Extract Q-matrices (ancestry proportions) for selected K values
- **Analysis Method**: Projected Least Squares with spatial constraints
- **Outputs**:
  - Cross-validation curves for model selection
  - Q-matrices for K=2,3,4,5 (ancestry proportions per individual)
  - Sample location maps with geographic distribution
  - STRUCTURE-like plots (generated from Q-matrices)
- **Interpretation**: Identifies genetically homogeneous subpopulations with geographic coherence

---

## Data Availability

Raw genomic data (reads, assemblies) and phenotypic data are available through [specify repository/database]. Processed data files needed to run these scripts are included in each folder.

## Contact

For questions about these analyses, please contact Renan Souza at rsouza@hudsonalpha.org.
