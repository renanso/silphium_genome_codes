################################################################################
# Phenotypic Correlations and Heritability Analysis
################################################################################
# Purpose: Calculate phenotypic correlations, test significance, and create
#          correlation plots for multiple traits across different environments
# Input: final_phenos2.csv - Phenotypic data with factors and traits
# Output: Correlation matrices, p-value matrices, and visualization plots
# Traits analyzed: RFC (Ray Floret Count), SM (Seed Mass), SNC (Seed Number 
#                  per Capitulum), RD (Receptacle Diameter)
################################################################################

# Clear workspace
gc()
rm(list = ls())

# ============================================================================
# LOAD REQUIRED LIBRARIES
# ============================================================================
library("Hmisc")           # Statistical tools including correlation functions
library("missMDA")         # Imputation for missing multivariate data
library("corrplot")        # Visualization of correlation matrices
library("VCA")             # Variance component analysis

# ============================================================================
# SECTION 1: DATA INPUT AND STRUCTURE INSPECTION
# ============================================================================
# Load phenotypic data with proper encoding for special characters
# Includes factors (clone, environment, rep, block) and quantitative traits
data <- read.csv("final_phenos2.csv", stringsAsFactors = FALSE, 
                 fileEncoding = "latin1", na.strings = NA)

# Examine data structure
head(data)
tail(data)
str(data)

# ============================================================================
# SECTION 2: DATA TYPE CONVERSION
# ============================================================================
# Columns 1-21: Categorical factors (clone, environment, replication, block, etc.)
# Convert to factors for proper analysis
data[, 1:21] <- lapply(data[, 1:21], factor)

# Columns 22-27: Quantitative traits (RFC, SM, SNC, RD, AI, HI)
# Convert to numeric for statistical analysis
data[, 22:27] <- lapply(data[, 22:27], as.numeric)

# Verify conversions
str(data)
summary(data)

# Attach data frame for direct column access
attach(data)

# ============================================================================
# SECTION 3: EXPLORATORY DATA VISUALIZATION
# ============================================================================
# Boxplots for trait distributions by environment
par(mar = c(1, 1, 1, 1))
par(mfrow = c(2, 3))
for (i in 1:length(data[, 22:27])) {
  boxplot(data[, 22:27][, i], main = names(data[, 22:27][i]), type = "l")
}
dev.off()

# Histograms for trait distributions
par(mar = c(1, 1, 1, 1))
par(mfrow = c(2, 3))
for (i in 1:length(data[, 22:27])) {
  hist(data[, 22:27][, i], main = names(data[, 22:27][i]))
}
dev.off()

# ============================================================================
# SECTION 4: DATA SUBSETTING BY ENVIRONMENT
# ============================================================================
# Create separate datasets for each environment (location x year combination)
# Select relevant columns: factors + key quantitative traits
data0a <- data[data$env == "al24", c(1:23, 26, 27)]  # Alabama 2024
data0b <- data[data$env == "ks24", c(1:23, 26, 27)]  # Kansas 2024
data0c <- data[data$env == "al25", 1:ncol(data)]     # Alabama 2025
data0d <- data[data$env == "ks25", 1:ncol(data)]     # Kansas 2025

# ============================================================================
# SECTION 5: CALCULATE PHENOTYPIC CORRELATIONS BY ENVIRONMENT
# ============================================================================
# Use Pearson correlation on quantitative traits
# rcorr() returns both correlation coefficients and p-values

# Alabama 2024: 4 traits (columns 22:25)
correlations_al24 <- rcorr(as.matrix(data0a[, c(22:25)]), type = "pearson")
write.csv(round(correlations_al24$r, 2), file = "phenotypic_correlations_al24.csv")
write.csv(correlations_al24$P, file = "phenotypic_corr_significance_al24.csv")

# Kansas 2024: 4 traits
correlations_ks24 <- rcorr(as.matrix(data0b[, c(22:25)]), type = "pearson")
write.csv(round(correlations_ks24$r, 2), file = "phenotypic_correlations_ks24.csv")
write.csv(correlations_ks24$P, file = "phenotypic_corr_significance_ks24.csv")

# Alabama 2025: 6 traits (columns 22:27)
correlations_al25 <- rcorr(as.matrix(data0c[, c(22:27)]), type = "pearson")
write.csv(round(correlations_al25$r, 2), file = "phenotypic_correlations_al25.csv")
write.csv(correlations_al25$P, file = "phenotypic_corr_significance_al25.csv")

# Kansas 2025: 6 traits
correlations_ks25 <- rcorr(as.matrix(data0d[, c(22:27)]), type = "pearson")
write.csv(round(correlations_ks25$r, 2), file = "phenotypic_correlations_ks25.csv")
write.csv(correlations_ks25$P, file = "phenotypic_corr_significance_ks25.csv")

# ============================================================================
# SECTION 6: PRINCIPAL COMPONENT IMPUTATION FOR MISSING VALUES
# ============================================================================
# Use PCA-based imputation to handle missing data while preserving structure
# Impute using 2 principal components
res.comp1 <- imputePCA(data0a[, c(22:25)], ncp = 2)
res.comp2 <- imputePCA(data0b[, c(22:25)], ncp = 2)
res.comp3 <- imputePCA(data0c[, c(22:27)], ncp = 2)
res.comp4 <- imputePCA(data0d[, c(22:27)], ncp = 2)

# Extract complete imputed datasets
cast_data1 <- data.frame(res.comp1$completeObs)
cast_data2 <- data.frame(res.comp2$completeObs)
cast_data3 <- data.frame(res.comp3$completeObs)
cast_data4 <- data.frame(res.comp4$completeObs)

# ============================================================================
# SECTION 7: CALCULATE CORRELATIONS FROM IMPUTED DATA
# ============================================================================
# Calculate correlation matrices for visualization
M1 <- cor(cast_data1)  # Alabama 2024
M2 <- cor(cast_data2)  # Kansas 2024
M3 <- cor(cast_data3)  # Alabama 2025
M4 <- cor(cast_data4)  # Kansas 2025

# ============================================================================
# SECTION 8: CALCULATE SIGNIFICANCE MATRIX FOR CORRELATIONS
# ============================================================================
# Custom function to calculate p-values for all pairwise correlations
# This allows visualization of statistical significance on correlation plots
cor.mtest <- function(mat, ...) {
  # Convert to matrix if needed
  mat <- as.matrix(mat)
  n <- ncol(mat)
  
  # Initialize p-value matrix
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  
  # Calculate p-value for each pairwise correlation
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  
  # Assign dimension names
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# Define color palette for visualization (red = negative, blue = positive)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

# ============================================================================
# SECTION 9: SAVE CORRELATION VISUALIZATION
# ============================================================================
# Create high-resolution image of correlation matrix
tiff("phenotypic_correlations.png", width = 25, height = 25, res = 200, units = "cm")

par(mfrow = c(2, 2))
p.mat1 <- cor.mtest(cast_data1)
corrplot(M1, method="color", col=col(200), 
         type="lower", order="original", 
         addCoef.col = "black",number.cex = 1.5, # Add coefficient of correlation
         tl.col="black", tl.srt=45, tl.cex = 1.5, #Text label color and rotation
         title = "AL24", mar = c(0, 0, 3, 0),
         cex.main = 2,
         # Combine with significance
         p.mat = p.mat1, sig.level = 0.001, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE) 
p.mat2 <- cor.mtest(cast_data2)
corrplot(M2, method="color", col=col(200), 
         type="lower", order="original", 
         addCoef.col = "black",number.cex = 1.5, # Add coefficient of correlation
         tl.col="black", tl.srt=45, tl.cex = 1.5, #Text label color and rotation
         title = "KS24", mar = c(0, 0, 3, 0),
         cex.main = 2,
         # Combine with significance
         p.mat = p.mat2, sig.level = 0.001, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE) 

p.mat3 <- cor.mtest(cast_data3)
corrplot(M3, method="color", col=col(200), 
         type="lower", order="original", 
         addCoef.col = "black",number.cex = 1.5, # Add coefficient of correlation
         tl.col="black", tl.srt=45, tl.cex = 1.5, #Text label color and rotation
         title = "AL25", mar = c(0, 0, 3, 0),
         cex.main = 2,
         # Combine with significance
         p.mat = p.mat3, sig.level = 0.001, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE) 

p.mat4 <- cor.mtest(cast_data4)
corrplot(M4, method="color", col=col(200), 
         type="lower", order="original", 
         addCoef.col = "black",number.cex = 1.5, # Add coefficient of correlation
         tl.col="black", tl.srt=45, tl.cex = 1.5, #Text label color and rotation
         title = "KS25", mar = c(0, 0, 3, 0),
         cex.main = 2,
         # Combine with significance
         p.mat = p.mat4, sig.level = 0.001, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE) 

dev.off()

############################################
#######BLUE correlations#######
############################################

#dataset filter to separate by env
#data0e<-rbind(data0c,data0d)

#including a G x E interaction factor for the combined analysis
data$ge<-as.factor(paste0(data$env,data$clone))
#data0e$ge<-as.factor(paste0(data0e$env,data0e$clone))
#str(data0e)

library(lme4)
fm_test1.1<-lmer(RFC~ clone + (1|env/rep) + (1|rep/block) + (1|ge), data=data)
rfc_blue<-data.frame(fixef(fm_test1.1))
rfc_blue$rfc <- rfc_blue$fixef.fm_test1.1. + rfc_blue[1,1]

fm_test1.2<-lmer(SM~ clone + (1|env/rep) + (1|rep/block) + (1|ge), data=data)
sm_blue<-data.frame(fixef(fm_test1.2))
sm_blue$sm <- sm_blue$fixef.fm_test1.2. + sm_blue[1,1]

fm_test1.3<-lmer(SNC~ clone + (1|env/rep) + (1|rep/block) + (1|ge), data=data)
snc_blue<-data.frame(fixef(fm_test1.3))
snc_blue$snc <- snc_blue$fixef.fm_test1.3. + snc_blue[1,1]

fm_test1.4<-lmer(RD~ clone + (1|env/rep) + (1|rep/block) + (1|ge), data=data)
rd_blue<-data.frame(fixef(fm_test1.4))
rd_blue$rd <- rd_blue$fixef.fm_test1.4. + rd_blue[1,1]

rfc_blue$name<- row.names(rfc_blue)
sm_blue$name<- row.names(sm_blue)
snc_blue$name<- row.names(snc_blue)
rd_blue$name<- row.names(rd_blue)

rfc_blue$name <- gsub("clone", "", rfc_blue$name)
sm_blue$name <- gsub("clone", "", sm_blue$name)
snc_blue$name <- gsub("clone", "", snc_blue$name)
rd_blue$name <- gsub("clone", "", rd_blue$name)

table1 <- merge(rfc_blue,sm_blue,  by = "name", all = TRUE)
table2 <- merge(snc_blue,rd_blue,  by = "name", all = TRUE)

table3 <- merge(table1,table2,  by = "name", all = TRUE)

table4<- table3[-1,c(1,3,5,7,9)]
colnames(table4)<-c("clone","RFC","SM","SNC","RD")
him<- read.csv("./env_index/hit_og.csv", header = T)
ai<- read.csv("./env_index/ai_summer_mean.csv", header = T)


table5 <- merge(table4,ai,  by = "clone", all = TRUE)
table6 <- merge(table5,him,  by = "clone", all = TRUE)

table7<- table6 [complete.cases(table6$AI), ]

table8<-table7[-c(61,94,146:148,158,201),]

## Impute NAS
res.comp <- imputePCA(table8[,c(2:7)],ncp=2)
comp_data<- data.frame(res.comp$completeObs)

##plot

M<-cor(comp_data)

corrplot(M, type="lower", order="hclust", tl.col="black", tl.srt=45)

##add significance
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(comp_data)
head(p.mat[, 1:5])

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
tiff("blue_correlations.png", width = 20, height = 20, res = 100, units = "cm")
corrplot(M, method="color", col=col(200), 
         type="lower", order="original", 
         addCoef.col = "black",number.cex = 1.5, # Add coefficient of correlation
         tl.col="black", tl.srt=45, tl.cex = 1.5,
         #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.001, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)
dev.off()


################################################
####variance component analysis and heritability
################################################
str(data)
#trait variance components
rfc<-lmer(RFC~ 1 + (1|clone) + (1|env/rep) + (1|ge), data=data)
asm<-lmer(SM~ 1 + (1|clone) + (1|env/rep) + (1|ge), data=data)
snc<-lmer(SNC~ 1 + (1|clone) + (1|env/rep) + (1|ge), data=data)
rd<-lmer(RD~ 1 + (1|clone) + (1|env/rep) + (1|ge), data=data)


vcs<-list(VarCorr(rfc),
          VarCorr(asm),VarCorr(snc),VarCorr(rd))

vas<-data.frame(print(vcs,comp=c("Variance")))
vas<-cbind(vas[,c(1,4,9,14,19)])

colnames(vas)<-c("component", "ray_floret_count", "average_seed_mass", "seed_number_per_capitulum", 
                 "receptacle_diameter") 

##Entry-mean based heritability equation: h2=v_RIL/[v_RIL + (v_RIL*E/E) + (v_e/RE)] 

#heritability rfc (4 environments, 2 reps each)
he_rfc<-vas[2,2]/(vas[2,2]+(vas[1,2]/4 + vas[5,2]/8))
he_rfc

#heritability asm (4 environments, 2 reps each)
he_asm<-vas[2,3]/(vas[2,3]+(vas[1,3]/4 + vas[5,3]/8))
he_asm

#heritability snc (2 environments, 2 reps each)
he_snc<-vas[2,4]/(vas[2,4]+(vas[1,4]/2 + vas[5,4]/4))
he_snc

#heritability rd (2 environments, 2 reps each)
he_rd<-vas[2,5]/(vas[2,5]+(vas[1,5]/2 + vas[5,5]/4))
he_rd

#heritability table
her<-c("heritability",he_rfc,he_asm,he_snc,he_rd)

her_table<-rbind(vas,her)

write.csv(her_table, file = "heritability_table.csv")


################################################
####variance component significance
################################################

#replication nested within the environment
# testing factors significance
#fit <- anovaVCA(protein~ env + name + env:rep + ge, data5)

fit5 <- anovaMM(RFC ~ clone + env + env/rep + ge, data)
fit6 <- anovaMM(SM ~ clone + env + env/rep + ge, data)
fit8 <- anovaMM(SNC ~ clone + env + env/rep + ge, data)
fit9 <- anovaMM(RD ~ clone + env + env/rep + ge, data)

fit5
fit6
fit8
fit9


############################
####Important QTL variances and effects
###single QTL model
#############################
###Variance Components: treating the QTL as a random effect

m1<-lmer(SNC~ 1 + (1|Chr01_695173749), data=data)
#m2<-lmer(RFC~ 1 + (1|Chr01_695173749), data=data)
#m3<-lmer(HI~ 1 + (1|Chr01_347199779), data=data)
#m4<-lmer(AI~ 1 + (1|Chr02_892024560), data=data)
#m5<-lmer(RFC~ 1 + (1|Chr04_222292225), data=data)
#m6<-lmer(HI~ 1 + (1|Chr03_78530134), data=data)
#m8<-lmer(AI~ 1 + (1|Chr06_274305900), data=data)
#m9<-lmer(RD~ 1 + (1|Chr01_1583902655), data=data)
#m10<-lmer(RFC~ 1 + (1|Chr01_1583902655), data=data)
#m11<-lmer(SNC~ 1 + (1|Chr01_1583902655), data=data)

vc <- as.data.frame(VarCorr(m1))
total_var <- sum(vc$vcov, na.rm = TRUE)
vc$pct_variance <- round(100 * vc$vcov / total_var, 4)
vc ##the variance for the markers is too high with these models.

#effects
ranef(m1)
coef(m1)
summary(m1)
coef(summary(m1)) #intercept

##A more precise way to determine the QTL effect is obtaining the R2 from the
## regression with the different allele combinations.
m1b<-lm(SNC~ Chr01_695173749, data=data)
r_squared_value <- summary(m1b)$adj.r.squared
r_squared_value

###################
##The R2 can also be calculated by the ratio of ss_pred / total_ss
##variance calculation with the 08_run_lmm_merged_traits.R script
##this calculation was based on the SS from ANOVA.######

m1.1<-lm(SNC~ as.factor(Chr01_695173749), data=data)
s <- summary(m1.1)
an <- anova(m1.1)
coef_df <- as.data.frame(coef(s))
coef_df$term <- rownames(coef_df)
rownames(coef_df) <- NULL

# Calculate variance explained from ANOVA (sum of squares)
an_df <- as.data.frame(an)
an_df$term <- rownames(an_df)
rownames(an_df) <- NULL

# find the predictor row (contains 'geno')
pred_row <- grep("Chr01_695173749", an_df$term, ignore.case = TRUE)
resid_row <- grep("Residuals", an_df$term, ignore.case = TRUE)
length(pred_row) 
length(resid_row)

ss_pred <- an_df[pred_row, "Sum Sq"]
ss_resid <- an_df[resid_row, "Sum Sq"]
total_ss <- ss_pred + ss_resid

pct <- round(100 * ss_pred / total_ss, 6)

variance_df <- data.frame(term = an_df$term[c(pred_row, resid_row)],
                          df = an_df$Df[c(pred_row, resid_row)],
                          sumsq = an_df$`Sum Sq`[c(pred_row, resid_row)],
                          pct_variance = c(pct, round(100 * ss_resid / total_ss, 6)))
 
#### effect calculation
###Mean value of the genotype

eff_stats <- aggregate(SNC ~ Chr01_695173749, data = data, FUN = function(x) c(n = length(x), mean = mean(x), sd = sd(x)))
eff_df <- do.call(data.frame, eff_stats)
colnames(eff_df) <- c("geno", "n", "mean", "sd")
eff_df$se <- with(eff_df, ifelse(n > 0, sd / sqrt(n), NA)) 


