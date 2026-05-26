################################################################################
# Gene-by-Environment Interaction Analysis
################################################################################
# Purpose: Analyze genotype-by-environment (GxE) interactions for quantitative
#          traits using linear models; create interaction plots showing trait
#          responses across environments and genotypes
# Input: final_phenos2.csv - Phenotypic data with genotypes and environments
# Output: Model summaries, ANOVA tables, interaction plots for each trait
# Traits: RFC (Ray Floret Count), SM (Seed Mass), SNC (Seed Number),
#         RD (Receptacle Diameter)
################################################################################

# Clear workspace
gc()
rm(list = ls())

# ============================================================================
# LOAD REQUIRED LIBRARIES
# ============================================================================
library("ggplot2")      # Grammar of graphics for publication-quality plots
library("emmeans")      # Estimated marginal means and interaction plots
library("dplyr")        # Data manipulation (group_by, summarise, mutate)
library("ggpubr")       # Publication-ready ggplot2 theme extensions

# ============================================================================
# SECTION 1: DATA INPUT AND STRUCTURE
# ============================================================================
# Load phenotypic data with proper encoding
data <- read.csv("final_phenos2.csv", stringsAsFactors = FALSE, 
                 fileEncoding = "latin1", na.strings = NA)

# Examine data structure
head(data)
tail(data)
str(data)

# Convert categorical variables (columns 1-21) to factors
data[, 1:21] <- lapply(data[, 1:21], factor)

# Convert quantitative traits (columns 22-27) to numeric
data[, 22:27] <- lapply(data[, 22:27], as.numeric)

# Verify conversions
str(data)
summary(data)

# Attach data frame for direct column access
attach(data)

# ============================================================================
# SECTION 2: FIT LINEAR MODELS WITH GENOTYPE-BY-ENVIRONMENT INTERACTIONS
# ============================================================================
# Model formula: Trait ~ clone + rep + env + clone:rep + clone:env
# This partitions variance into main effects and interaction terms

# Ray Floret Count model
fm1 <- lm(RFC ~ clone * rep * env, data = data)
summary(fm1)
anova(fm1)

# Seed Mass model
fm2 <- lm(SM ~ clone * rep * env, data = data)
summary(fm2)
anova(fm2)

# Seed Number per Capitulum model
fm3 <- lm(SNC ~ clone * rep * env, data = data)
summary(fm3)
anova(fm3)

# Receptacle Diameter model
fm4 <- lm(RD ~ clone * rep * env, data = data)
summary(fm4)
anova(fm4)

# ============================================================================
# SECTION 3: CALCULATE TRAIT MEANS BY ENVIRONMENT
# ============================================================================
# Calculate mean phenotypic values across all genotypes for each environment
# These means are overlaid on genotype-specific trends for visual reference

# Ray Floret Count means by environment
env_means_rfc <- data %>%
  group_by(env) %>%
  summarise(mean_value = mean(RFC, na.rm = TRUE)) %>%
  # Create numeric environment index for plotting
  mutate(env_num = as.numeric(factor(env)),
         trait = "RFC")

# Seed Mass means by environment
env_means_sm <- data %>%
  group_by(env) %>%
  summarise(mean_value = mean(SM, na.rm = TRUE)) %>%
  mutate(env_num = as.numeric(factor(env)),
         trait = "SM")

# Seed Number per Capitulum means
# Note: Only 2 environments have SNC data (al25 and ks25)
env_means_snc <- data %>%
  group_by(env) %>%
  summarise(mean_value = mean(SNC, na.rm = TRUE)) %>%
  mutate(env_num = as.numeric(factor(env)),
         trait = "SNC")
# Subset and renumber environments
env_means_snc <- env_means_snc[c(2, 4), ]
env_means_snc$env_num <- c(1, 2)

# Receptacle Diameter means (2 environments: al25, ks25)
env_means_rd <- data %>%
  group_by(env) %>%
  summarise(mean_value = mean(RD, na.rm = TRUE)) %>%
  mutate(env_num = as.numeric(factor(env)),
         trait = "RD")
# Subset and renumber environments
env_means_rd <- env_means_rd[c(2, 4), ]
env_means_rd$env_num <- c(1, 2)

# ============================================================================
# SECTION 4: CREATE GENOTYPE-BY-ENVIRONMENT INTERACTION PLOTS
# ============================================================================
# Each plot shows trait values for each genotype across environments
# Light gray lines = individual genotypes
# Black line = population mean (for visual reference)

# Ray Floret Count interaction plot
plot_rfc <- emmip(fm1, clone ~ env) +
  theme(legend.position = 'none') +
  ylab("Ray Floret Count") +
  xlab(NULL) +
  theme(text = element_text(size = 18)) +
  # All genotype lines in light gray
  scale_colour_manual(values = rep("lightgray", 300)) +
  # Overlay population mean as solid black line
  geom_line(data = env_means_rfc, aes(x = env_num, y = mean_value),
            colour = "black", linewidth = 1.2, linetype = "solid", inherit.aes = FALSE)

# Seed Mass interaction plot
plot_sm <- emmip(fm2, clone ~ env) +
  theme(legend.position = 'none') +
  ylab("Seed Mass (mg)") +
  xlab(NULL) +
  theme(text = element_text(size = 18)) +
  scale_colour_manual(values = rep("lightgray", 300)) +
  geom_line(data = env_means_sm, aes(x = env_num, y = mean_value),
            colour = "black", linewidth = 1.2, linetype = "solid", inherit.aes = FALSE)

# Seed Number per Capitulum interaction plot
plot_snc <- emmip(fm3, clone ~ env) +
  theme(legend.position = 'none') +
  ylab("Seed Number per Capitulum") +
  xlab(NULL) +
  theme(text = element_text(size = 18)) +
  scale_colour_manual(values = rep("lightgray", 300)) +
  geom_line(data = env_means_snc, aes(x = env_num, y = mean_value),
            colour = "black", linewidth = 1.2, linetype = "solid", inherit.aes = FALSE)

# Receptacle Diameter interaction plot
plot_rd <- emmip(fm4, clone ~ env) +
  theme(legend.position = 'none') +
  ylab("Receptacle Diameter (mm)") +
  xlab(NULL) +
  theme(text = element_text(size = 18)) +
  scale_colour_manual(values = rep("lightgray", 300)) +
  geom_line(data = env_means_rd, aes(x = env_num, y = mean_value),
            colour = "black", linewidth = 1.2, linetype = "solid", inherit.aes = FALSE)

# ============================================================================
# SECTION 5: DISPLAY PLOTS
# ============================================================================
# Print Ray Floret Count plot to console/graphics window
plot_rfc
plot_sm
plot_snc
plot_rd

tiff("g_x_e_plots.jpeg", width = 30, height = 30, res = 100, units = "cm")
ggarrange(plot_rfc,plot_sm,plot_snc,plot_rd, ncol = 2, nrow = 2)
dev.off()
