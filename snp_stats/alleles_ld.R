rm(list=ls())
library(SNPRelate)

# Read the alleles.csv file
data <- read.csv("alleles.csv", header=TRUE, stringsAsFactors=FALSE)

# Function to convert genotype string to numeric (0,1,2)
convert_geno <- function(geno, alleles) {
  a1 <- substr(geno, 1, 1)
  a2 <- substr(geno, 2, 2)
  ref <- alleles[1]
  alt <- alleles[2]
  if (a1 == ref && a2 == ref) return(0)
  if (a1 == alt && a2 == alt) return(2)
  return(1)
}

# Get alleles for each locus
alleles_list <- lapply(data, function(x) sort(unique(unlist(strsplit(x, "")))))

# Create genotype matrix (rows: samples, columns: SNPs)
geno_matrix <- matrix(0, nrow = nrow(data), ncol = ncol(data))
for (i in 1:ncol(data)) {
  alleles <- alleles_list[[i]]
  geno_matrix[, i] <- sapply(data[, i], convert_geno, alleles = alleles)
}

# Create GDS file (transpose for SNPs as rows)
snpgdsCreateGeno("alleles.gds",
  genmat = t(geno_matrix),
  sample.id = paste0("sample", 1:nrow(data)),
  snp.id = colnames(data),
  snp.position = 1:ncol(data),
  snp.chromosome = rep(1, ncol(data)),
  snp.allele = sapply(alleles_list, function(a) paste(a, collapse = "/"))
)

# Open the GDS file
genofile <- snpgdsOpen("alleles.gds")

# Calculate LD matrix
ld <- snpgdsLDMat(genofile, slide = -1, method = "r")

# Extract r2 matrix
r2_matrix <- ld$LD^2
rownames(r2_matrix) <- colnames(data)
colnames(r2_matrix) <- colnames(data)

r2_matrix

row.names(r2_matrix)<- c("Chr01_695173749\nRFC/SNC","Chr01_1583902655\nRFC/RD/SNC","Chr03_721440854\nRFC/RD",
                         "Chr04_222292225\nRFC", "Chr06_274305900\nAI","Chr03_78530134\nHI")

colnames(r2_matrix)<- c("Chr01_695173749\n  RFC/SNC","Chr01_1583902655\n  RFC/RD/SNC","Chr03_721440854\n  RFC/RD",
                         "Chr04_222292225\n RFC", "Chr06_274305900\n  AI","Chr03_78530134\n  HI")

# plot
#image(t(r2_matrix), col=terrain.colors(16))

library(corrplot)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(r2_matrix, method="color", col=col(200), 
         type="lower", order="original", 
         addCoef.col = "black",number.cex = 0.8, # Add coefficient of correlation
         tl.col="black", tl.srt=45, tl.cex = 0.8, #Text label color and rotation
         # Combine with significance
         #p.mat = p.mat, sig.level = 0.001, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

# Save the r2 matrix to a CSV file
write.csv(r2_matrix, "ld_r2_matrix.csv")


