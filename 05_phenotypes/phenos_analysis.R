################################################################################
# Phenotypic Data Analysis and Linear Mixed Effects Modeling
################################################################################
# Purpose: Analyze phenotypic variation using mixed effects models
#          to partition genetic and environmental effects; calculate
#          Best Linear Unbiased Estimates (BLUEs) for each genet
# Input: phenotypic_data.csv - Raw phenotypic measurements from field trials
# Output: Summary statistics, mixed effects model results, BLUE predictions
# Traits: max_ray_floret_count, ave_sd_mass_mg, seed_number_per_head,
#         largest_recp_diam_mm
################################################################################

# Load required library for mixed effects models
library("lme4")

# Clear workspace
rm(list = ls())

# ============================================================================
# SECTION 1: DATA INPUT AND INITIAL INSPECTION
# ============================================================================
# Load raw phenotypic data from field trials
data <- read.csv("phenotypic_data.csv", stringsAsFactors = FALSE, fileEncoding = "latin1")

# Examine data structure
head(data)
tail(data)
str(data)

# ============================================================================
# SECTION 2: DATA SUBSETTING AND TYPE CONVERSION
# ============================================================================
# Select relevant columns for analysis (first 5 plus traits of interest)
# Reduces dimensionality while retaining key variables
data2 <- data[, c(1:5, 7, 16, 17, 20:23)]

# Convert first 8 columns to factors (experimental design variables)
# These represent: clone, environment, replication, block, etc.
data2[, 1:8] <- lapply(data2[, 1:8], factor)
str(data2)

# Convert last 4 columns to numeric (quantitative traits)
data2[, 9:12] <- lapply(data2[, 9:12], as.numeric)
str(data2)
summary(data2)

# ============================================================================
# SECTION 3: DATA PERSISTENCE
# ============================================================================
# Save cleaned data as R object (faster loading, preserves structure)
saveRDS(data2, "pheno")

# Load and attach phenotype data for analysis
pheno <- readRDS("pheno")
attach(pheno)
str(pheno)

# ============================================================================
# SECTION 4: EXPLORATORY DATA VISUALIZATION
# ============================================================================
# Boxplots for trait distributions across experimental groups
par(mar = c(1, 1, 1, 1))
par(mfrow = c(2, 2))
for (i in 1:length(pheno[, 9:12])) {
  boxplot(pheno[, 9:12][, i], main = names(pheno[, 9:12][i]), type = "l")
}
dev.off()

# Histograms for trait distributions (assess normality, identify outliers)
par(mar = c(1, 1, 1, 1))
par(mfrow = c(2, 2))
for (i in 1:length(pheno[, 9:12])) {
  hist(pheno[, 9:12][, i], main = names(pheno[, 9:12][i]))
}
dev.off()

# Save clean phenotypes for downstream analysis
saveRDS(pheno, "pheno2")

# ============================================================================
# SECTION 5: RELOAD DATA AND PREPARE FOR ANALYSIS
# ============================================================================
# Clear workspace and reload cleaned data
rm(list = ls())

# Load the clean phenotype dataset
data <- readRDS("pheno2")
attach(data)
str(data)

# ============================================================================
# SECTION 6: INITIAL TRAIT ANALYSIS
# ============================================================================
# Visualize phenotypic variation by clone (main genetic effect)
boxplot(data$max_ray_floret_count ~ clone)
boxplot(data$ave_sd_mass_mg ~ clone)
boxplot(data$seed_number_per_head ~ clone)
boxplot(data$largest_recp_diam_mm ~ clone)

# ============================================================================
# SECTION 7: SUBSETTING BY ENVIRONMENT
# ============================================================================
# Create separate datasets for each environment (location-year combination)
# This allows environment-specific analysis
data0a <- data[data$env == "al24", 1:ncol(data)]  # Alabama 2024
data0b <- data[data$env == "ks24", 1:ncol(data)]  # Kansas 2024
data0c <- data[data$env == "al25", 1:ncol(data)]  # Alabama 2025
data0d <- data[data$env == "ks25", 1:ncol(data)]  # Kansas 2025

# Create genotype-by-environment factor for combined analysis
# This allows modeling of GxE interactions
data$ge <- as.factor(paste0(data$env, data$clone))
str(data)

# ============================================================================
# SECTION 8: DESCRIPTIVE STATISTICS
# ============================================================================
# Generate summary statistics for all traits across genotypes
summary_all <- data.frame(summary(data[9:12], maxsum = 12))
write.csv(summary_all, file = "all_traits_summary.csv")

# ============================================================================
# SECTION 9: LINEAR MIXED EFFECTS MODELS
# ============================================================================
# Fit mixed effects models for combined analysis (all environments)
# Model structure: Trait ~ Fixed effects (clone) + Random effects (env, rep, block, GxE)

# Model 1: Ray Floret Count
# Random effects account for variation due to location, replication, and GxE
fm_test1.1 <- lmer(max_ray_floret_count ~ clone + (1|env/rep) + (1|rep/block) + (1|ge), 
                    data = data)
# Extract BLUEs (conditional means) for each clone
coef(fm_test1.1)[[2]]          # Random effect of clone (Intercept + BLUP)
ranef(fm_test1.1)              # Random effects summary
fix <- fixef(fm_test1.1)        # Fixed effects estimates
summary(fm_test1.1)             # Full model summary
coef(summary(fm_test1.1))       # Fixed effect coefficient table

# Model 2: Seed Mass (mg)
fm_test1.2 <- lmer(ave_sd_mass_mg ~ clone + (1|env/rep) + (1|rep/block) + (1|ge), 
                    data = data)
coef(fm_test1.2)[[2]]
ranef(fm_test1.2)
fix <- fixef(fm_test1.2)
summary(fm_test1.2)
coef(summary(fm_test1.2))

# Model 3: Seed Number per Head
fm_test1.3 <- lmer(seed_number_per_head ~ clone + (1|env/rep) + (1|rep/block) + (1|ge), 
                    data = data)
coef(fm_test1.3)[[2]]
ranef(fm_test1.3)
fix <- fixef(fm_test1.3)
summary(fm_test1.3)
coef(summary(fm_test1.3)) #fixed effect

fm_test1.4<-lmer(largest_recp_diam_mm~ clone + (1|env/rep) + (1|rep/block) + (1|ge), data=data)
coef(fm_test1.4)[[2]] #random effect of clone (Intercept + BLUP)
ranef(fm_test1.4)
fix<-fixef(fm_test1.4)
summary(fm_test1.4)
coef(summary(fm_test1.4)) #fixed effect

##combining all the results
ray_floret_count<- as.data.frame(coef(summary(fm_test1.1)))
ave_sd_mass<- as.data.frame(coef(summary(fm_test1.2)))
seed_number_per_head<- as.data.frame(coef(summary(fm_test1.3)))
largest_recp_diam_mm<- as.data.frame(coef(summary(fm_test1.4)))

write.csv(ray_floret_count, "ray_floret_count_comb.csv")
write.csv(ave_sd_mass, "ave_sd_mass_comb.csv")
write.csv(seed_number_per_head, "seed_number_per_head_comb.csv")
write.csv(largest_recp_diam_mm, "largest_recp_diam_mm_comb.csv")

############ALABAMA 2024###############
fm_test2.1<-lmer(max_ray_floret_count~ clone + (1|rep), data=data0a)
coef(fm_test2.1)[[1]] #random effect of clone (Intercept + BLUP)
ranef(fm_test2.1)
fixef(fm_test2.1)
summary(fm_test2.1)
coef(summary(fm_test2.1)) #fixed effect

fm_test2.2<-lmer(ave_sd_mass_mg~ clone + (1|rep), data=data0a)
coef(fm_test2.2)[[1]] #random effect of clone (Intercept + BLUP)
ranef(fm_test2.2)
fixef(fm_test2.2)
summary(fm_test2.2)
coef(summary(fm_test2.2)) #fixed effect

##combining all the results
ray_floret_count<- as.data.frame(coef(summary(fm_test2.1)))
ave_sd_mass<- as.data.frame(coef(summary(fm_test2.2)))

write.csv(ray_floret_count, "ray_floret_count_al24.csv")
write.csv(ave_sd_mass, "ave_sd_mass_al24.csv")

############KANSAS 2024###############
fm_test3.1<-lmer(max_ray_floret_count~ clone + (1|rep/block), data=data0b)
coef(fm_test3.1)[[1]] #random effect 
ranef(fm_test3.1)
fixef(fm_test3.1)
summary(fm_test3.1)
coef(summary(fm_test3.1)) #fixed effect

fm_test3.2<-lmer(ave_sd_mass_mg~ clone + (1|rep/block), data=data0b)
coef(fm_test3.2)[[1]] #random effect 
ranef(fm_test3.2)
fixef(fm_test3.2)
summary(fm_test3.2)
coef(summary(fm_test3.2)) #fixed effect

##combining all the results
ray_floret_count<- as.data.frame(coef(summary(fm_test3.1)))
ave_sd_mass<- as.data.frame(coef(summary(fm_test3.2)))

write.csv(ray_floret_count, "ray_floret_count_ks24.csv")
write.csv(ave_sd_mass, "ave_sd_mass_ks24.csv")

############ALABAMA 2025###############
colnames(data0c)
traits4<-colnames(data0c[,c(9:12)])
traits4
str(data0c)
dataset4<-list(data0c)
result4<-list()

for (i in traits4) {
  model <- lapply(dataset4, 
                  function(data0c){
                    lmer(reformulate(termlabels="clone + (1|rep/block)", response=i), data0c) })
  vals<-lapply(model,function(x) {coef(summary(x))})
  result4[[length(result4)+1]]<-vals
  print(i)
}

##combining all the results
ray_floret_count<- as.data.frame(result4[[1]][[1]])
ave_sd_mass<- as.data.frame(result4[[2]][[1]])
seed_number_per_head<- as.data.frame(result4[[3]][[1]])
largest_recp_diam_mm<- as.data.frame(result4[[4]][[1]])

write.csv(ray_floret_count, "ray_floret_count_al25.csv")
write.csv(ave_sd_mass, "ave_sd_mass_al25.csv")
write.csv(seed_number_per_head, "seed_number_per_head_al25.csv")
write.csv(largest_recp_diam_mm, "largest_recp_diam_mm_al25.csv")

############KANSAS 2025###############
colnames(data0d)
traits5<-colnames(data0d[,c(9:12)])
traits5
str(data0d)
dataset5<-list(data0d)
result5<-list()

for (i in traits5) {
  model <- lapply(dataset5, 
                  function(data0d){
                    lmer(reformulate(termlabels="clone + (1|rep/block)", response=i), data0d) })
  vals<-lapply(model,function(x) {coef(summary(x))})
  result5[[length(result5)+1]]<-vals
  print(i)
}

##combining all the results
ray_floret_count<- as.data.frame(result5[[1]][[1]])
ave_sd_mass_mg<- as.data.frame(result5[[2]][[1]])
seed_number_per_head<- as.data.frame(result5[[3]][[1]])
largest_recp_diam_mm<- as.data.frame(result5[[4]][[1]])

write.csv(ray_floret_count, "ray_floret_count_ks25.csv")
write.csv(ave_sd_mass_mg, "ave_sd_mass_mg_ks25.csv")
write.csv(seed_number_per_head, "seed_number_per_head_ks25.csv")
write.csv(largest_recp_diam_mm, "largest_recp_diam_mm_ks25.csv")
