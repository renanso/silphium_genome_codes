rm(list = ls())
#devtools::install_github("jiabowang/GAPIT3",force=TRUE)
#options(download.file.method="libcurl", url.method="libcurl")
#source("http://www.zzlab.net/GAPIT/GAPIT.library.R")
#source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
library(GAPIT)
library(BiocParallel)
bplapply(1:10, print, BPPARAM = MulticoreParam (workers = 8))
###fast reading dataframe
library(data.table)

## genotyping information
myG<-fread(file="011_silphium_smiss1_miss03_maf01_449K_258n_updated_2n_std.Ihapmap",  sep = "\t") # reading only the first column of the file

###File QC###
#removing 1st column
myG<-myG[2:nrow(myG),]

str(myG)

#selecting SNPs with only two alleles
myG<-subset(myG, nchar(as.character(alleles)) <= 3) 

##subset in two parts to compare with the phenotyping information

# general information
myG1<-myG[,1:11]

# actual snp calls
myG2<-myG[,12:ncol(myG)]
myG2[1:20,1:30]

## phenotyping information
myY<- fread("rfc.txt", sep = "\t")
myY2<-data.frame(myY[myY$Taxa %in% colnames(myG2),])


## selecting genotyping information of lines with phenotyping information
myG3<-myG2[, colnames(myG2) %in% myY2$Taxa, with=FALSE]

## reorder
col_order <- c(colnames(myG3))
myY3<-myY2[match(col_order, myY2$Taxa),]

#checking if the names match
all(myY3$Taxa == colnames(myG3))
myG4<- data.frame(myG1,myG3)
str(myY3)
str(myG4)
myG4[1,]<-names(myG4)
names(myG4) <- NULL

#test line
geno_test<- myG4[1:50000, 1:ncol(myG4)]

#GWAS with 2 models methods
myGAPIT <- GAPIT(
 Y=myY3,
  G=myG4,
  PCA.total=3,
  model=c("FarmCPU", "Blink"),
  NJtree.group = 3,
  Multiple_analysis=TRUE)
