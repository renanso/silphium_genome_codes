library("adegenet")
library("hierfstat")
library("pegas")
library("data.table")

##replace test.hapmap with the actual hapmap
data1 <- read.table("test.hapmap", header = TRUE)

data1$snp<- paste0(data1$chr,"_",data1$pos)
rownames(data1)<- data1$snp
n<-ncol(data1)-1

data2<- t(data1[,3:n])
tree_id<- row.names(data2)
data3<-cbind(tree_id,data2)

info <- read.csv("info_nj_groups2.csv", header = TRUE)
info2<- info[,c(1,2,3,2)]
colnames(info2)<- c("tree_id", "country","state","country_state")

data4<- merge(data3,info2,by = "tree_id")
ncols<- ncol(data4)-3

new_cols<- c(1,(ncols+1),(ncols+2),(ncols+3),2:ncols)

data5<- data4[,new_cols]

#############

#str(data5)

locus <- data5[, -c(1, 2, 3, 4)]    

ind <- as.character(data5$tree_id) # labels of the individuals
population <- as.character(data5$state) # labels of the populations
Mydata1 <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep = "")
Mydata1

nAll(Mydata1) # Number of alleles per locus
saveRDS(Mydata1,"Mydata1")

div <- summary(Mydata1)
saveRDS(div,"div1")
#mean observed heterozygosity
mean(div$Hobs)

names(div)
jpeg(file="observed_het1.jpeg", width = 640, height = 640, units = "px",  quality = 100)
plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus")
dev.off()

jpeg(file="het-distribution1.jpeg", width = 640, height = 640, units = "px",  quality = 100)
hist(div$Hobs)
dev.off()









