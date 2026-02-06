gc()
rm(list = ls());ls ()

library("Hmisc")
library("missMDA") # impute NAs
library("corrplot") # plot correlations
library("VCA") # variance component significance
#loading the phenotypes file
data <- read.csv("final_phenos2.csv", stringsAsFactors=FALSE, fileEncoding="latin1", na.strings = NA)
head(data)
tail(data)
str(data)

# adjusting the factors and variables
data[,1:21] <- lapply(data[,1:21], factor)
data[,22:27] <- lapply(data[,22:27], as.numeric)
str(data)
summary(data)
attach(data)

# data visualization
par(mar=c(1,1,1,1))
par(mfrow=c(2,3))
for (i in 1:length(data[,22:27])) {
  boxplot(data[,22:27][,i], main=names(data[,22:27][i]), type="l")}
dev.off()

#quantitative
par(mar=c(1,1,1,1))
par(mfrow=c(2,3))
for (i in 1:length(data[,22:27])) {
  hist(data[,22:27][,i], main=names(data[,22:27][i]))}
dev.off()

#data subset

data0a<- data[data$env=="al24",c(1:23,26,27)]
data0b<- data[data$env=="ks24",c(1:23,26,27)]
data0c<- data[data$env=="al25",1:ncol(data)]
data0d<- data[data$env=="ks25",1:ncol(data)]

###correlations

##ALL ENV
correlations_al24<- rcorr(as.matrix(data0a[,c(22:25)]), type= "pearson")
write.csv(round(correlations_al24$r,2), file = "phenotypic_correlations_al24.csv")
write.csv(correlations_al24$P, file = "phenotypic_corr_significance_al24.csv")

correlations_ks24<- rcorr(as.matrix(data0b[,c(22:25)]), type= "pearson")
write.csv(round(correlations_ks24$r,2), file = "phenotypic_correlations_ks24.csv")
write.csv(correlations_ks24$P, file = "phenotypic_corr_significance_ks24.csv")

correlations_al25<- rcorr(as.matrix(data0c[,c(22:27)]), type= "pearson")
write.csv(round(correlations_al25$r,2), file = "phenotypic_correlations_al25.csv")
write.csv(correlations_al25$P, file = "phenotypic_corr_significance_al25.csv")

correlations_ks25<- rcorr(as.matrix(data0d[,c(22:27)]), type= "pearson")
write.csv(round(correlations_ks25$r,2), file = "phenotypic_correlations_ks25.csv")
write.csv(correlations_ks25$P, file = "phenotypic_corr_significance_ks25.csv")

##plots

res.comp1 <- imputePCA(data0a[,c(22:25)],ncp=2)
res.comp2 <- imputePCA(data0b[,c(22:25)],ncp=2)
res.comp3 <- imputePCA(data0c[,c(22:27)],ncp=2)
res.comp4 <- imputePCA(data0d[,c(22:27)],ncp=2)

cast_data1<- data.frame(res.comp1$completeObs)
cast_data2<- data.frame(res.comp2$completeObs)
cast_data3<- data.frame(res.comp3$completeObs)
cast_data4<- data.frame(res.comp4$completeObs)

M1<-cor(cast_data1)
M2<-cor(cast_data2)
M3<-cor(cast_data3)
M4<-cor(cast_data4)

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

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

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

m1<-lmer(SNC~ 1 + (1|Chr01_695173749), data=data)
m2<-lmer(RFC~ 1 + (1|Chr01_695173749), data=data)
#m3<-lmer(HI~ 1 + (1|Chr01_347199779), data=data)
#m4<-lmer(AI~ 1 + (1|Chr02_892024560), data=data)
m5<-lmer(RFC~ 1 + (1|Chr04_222292225), data=data)
#m6<-lmer(HI~ 1 + (1|Chr03_78530134), data=data)
#m8<-lmer(AI~ 1 + (1|Chr06_274305900), data=data)
m9<-lmer(RD~ 1 + (1|Chr01_1583902655), data=data)
m10<-lmer(RFC~ 1 + (1|Chr01_1583902655), data=data)
m11<-lmer(SNC~ 1 + (1|Chr01_1583902655), data=data)


vc <- as.data.frame(VarCorr(m1))
total_var <- sum(vc$vcov, na.rm = TRUE)
vc$pct_variance <- round(100 * vc$vcov / total_var, 4)
vc


