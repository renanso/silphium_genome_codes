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
**Purpose**: Visualize sample collection locations and environmental context. The main script in this directory is `climate_data.R`, which integrates geographic coordinates with WorldClim v2.1 climate variables. In addition to obtaining the data from WorldClim, the script generates maps showing Silphium sample collection sites and maps with index layers (Aridity and Heat).

---

### `gene_families/` — Gene Family Abundance Analysis
**Purpose**: Compare gene family copy numbers across plant genomes. The main script here is `gene_families_plot.R`, which analyzes abundance of selected gene families across 5 plant genomes (2 Silphium spp., 2 Helianthus spp., S. taccada).

---

### `genome_scans/` — Genome-Wide Association Study (GWAS) and Genome-environment analysis (GEA)
**Purpose**: Identify SNPs significantly associated with capitulum traits and environmental indexes. The script to run this analysis is `gwas_gapit.R` and it uses the R package GAPIT3. The input data are: 449K SNPs genotyped across 258 Silphium samples and the phenotypes: RFC (Ray Floret Count), SM (Seed Mass), SNC (Seed Number per Capitulum), RD (Receptacle Diameter), AI (Aridity Index) and HI (Heat Index). Two models were used: FarmCPU and BLINK.

---

### `heritability_correlations/` — Phenotypic Analysis and G×E Interactions
**Purpose**: Quantify trait heritability and genotype-by-environment interactions. The main scripts are:
`cor_var_her.R`: Phenotypic correlation and variance component analysis.
`gxe.R`: Genotype-by-environment interaction visualization.

The scripts generate correlation matrices, G X E interaction plots and variance components table.

---

### `nj_tree/` — Phylogenetic Tree and Population Structure
**Purpose**: Visualize genetic relationships and population substructure. The core script is `nj_tree.R` and it constructs neighbor-joining phylogenetic tree from a IBS kinship matrix and colors clades by species/subpopulation. 

---

### `pca/` — Principal Component Analysis (PCA) Visualization
**Purpose**: Explore population structure and geographic genetic patterns. The main script is `pca_plot.R`.It generates one PCA plot colored by colored by species assignment and one with geographic gradient showing east-west population structure.

---

### `phenotypes/` — Mixed Effects Phenotypic Modeling
**Purpose**: Data QC, BLUEs calculation and partition phenotypic variance. The main script is `phenos_analysis.R` and it performs data cleaning, exploratory visualization and calculation of Best Linear Unbiased Estimates (BLUEs). 

---

### `snp_stats/` — SNP Statistics and Genetic Diversity
**Purpose**: Analyze SNP patterns, linkage disequilibrium, and heterozygosity. The scripts in this directory are: 
- `alleles_ld.R`: Linkage disequilibrium between candidate SNPs.
- `heterozygosity1.R`: SNP heterozygosity estimation.
- `SNP_density.R`: SNP density plot visualization.

---

### `synteny/` — Comparative Genomics and Synteny Analysis
**Purpose**: Comparative genomics across studied species. The scripts here are built around the Genespace R package:
- `riparian_test.R`: GENESPACE pipeline for synteny identification. 
- `genespace_viz_test_cluster_2.R`: Genomic feature classification and visualization (genes, TIRs, Helitrons, LTRs). 

---

### `te_families/` — Transposable Element Analysis
**Purpose**: Characterize TE abundance and distribution across genomes. The script `plot_te.R` compares TE content across the 5 plant genomes using EDTA pipeline results. 

---

### `tess3r/` — Spatially-Informed Population Structure
**Purpose**: Identify population substructure incorporating geographic information using the Tess3r R package. The script performs TESS3 analysis testing K=1-8 ancestral populations and uses spatial priors to constrain ancestry assignments to geographic patterns. 

---

### `SLURM_scripts/` — SLURM job scripts
**Purpose**: This folder contains several helper scripts to perform QC and analyze sequencing data. These are all shell scripts to run in a HPC.

---

## Data Availability

Genomic data (reads, assemblies) for this research is available in the following NCBI bioprojects:
- PRJNA1417267: Raw sequencing data for genome assembly.
- PRJNA1426741: Sequencing data for the Silphium association panel.
- PRJNA1426744: Sequencing data for the Silphium phylogeny.
- PRJNA1422233: Final genome assembly for Silphium perfoliatum (1JHW).
- PRJNA1422234: Final genome assembly for Silphium integrifolium (Bad Astra).

Genome annotations, chloroplast assembly, annotation and supplementary data are available at Figshare https://figshare.com/s/f75787ae1adf175ad955.

Processed data files needed to run the scripts are included in each folder.

## Contact

For questions about these analyses, please contact Renan Souza at rsouza@hudsonalpha.org.
