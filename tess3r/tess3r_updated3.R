rm(list=ls())

library("Rcpp")
library("tess3r")
library("maps")
library("BiocParallel")
library("data.table")
library("ggplot2")
library("ggrepel")

bplapply(1:10, print, BPPARAM = MulticoreParam (workers = 8))

## genotyping information
##replace the test.txt file for the numeric hapmap
data2<-fread(file="test.txt",  sep = "\t")
colnames(data2)[1]<- "clone"

#data2<- read.table("RSouza_silphium.Ihapmap.150gps.numeric2.hmp.txt", header = T)
#colnames(data2)[1]<- "clone"

#coord2<-read.table("coordinates2.txt", header = T)
coord2<-read.table("coordinates_integrifolium.txt", header = T)

#merge to match genotypes with coordinates
n1<-ncol(data2) + 1
n2<- ncol(data2) + 2

coord3<- (merge(data2, coord2, by = 'clone'))[,n2:n1]
coord3<-as.matrix(coord3)

write.table(coord3, file = "coordinates_matrix.txt", sep = "\t", row.names = F, col.names = T)
data3<- (merge(data2, coord2, by = 'clone'))[,1:ncol(data2)]

##plot samples
png(filename="plot.png", width = 10000, height = 10000, res=300, units = "px")
#labels <- coord2$clone   # or df$SampleName, for example
#plot(coord2[,c(3,2)], pch = 19, cex = .5, 
#     xlab = "Longitude (°E)", ylab = "Latitude (°N)")
#map('state', add = T, interior =T)
# Add non-overlapping labels
#pointLabel(coord2[,3], coord2[,2], labels = labels, cex = 0.7)
#text(coord2[,3], coord2[,2], labels = labels,
#     pos = 4,     # position relative to the point: 1=below, 2=left, 3=above, 4=right
#     cex = 0.7,   # text size scaling
#     col = "blue")

# Convert coordinates to a data frame
df <- as.data.frame(coord3)
df$label <- labels

ggplot(coord2, aes(x = long, y = lat)) +
  geom_point(size = 1.5) +
  borders("state") +
  geom_text_repel(aes(label = clone), size = 3, color = "blue",
                  max.overlaps = Inf) +
  labs(x = "Longitude (°E)", y = "Latitude (°N)") +
  coord_quickmap(
    xlim = c(-100, -85),
    ylim = c(30, 45),
    expand = FALSE
  ) +
  theme_minimal(base_size = 14)
dev.off()

##preparing genotyping dataset for tess3r
rownames(data3)<- data3$clone
data3$clone<-NULL

tess3.obj <- tess3(X = data3, coord = coord3, K = 1:8, 
                   method = "projected.ls", ploidy = 2, openMP.core.num = 8) 

# Save an object to a file
saveRDS(tess3.obj, file = "tess3.rds")
# Restore the object
#readRDS(file = "tess3.rds")

png(filename="plot1.png", width = 480, height = 480, units = "px")
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")
dev.off()


# retrieve tess3 Q matrix for clusters 
q.matrix2 <- qmatrix(tess3.obj, K =  2)
q.matrix3 <- qmatrix(tess3.obj, K =  3)
q.matrix4 <- qmatrix(tess3.obj, K =  4)
q.matrix5 <- qmatrix(tess3.obj, K =  5)

saveRDS(q.matrix2, file = "q.matrix2.rds")
saveRDS(q.matrix3, file = "q.matrix3.rds")
saveRDS(q.matrix4, file = "q.matrix4.rds")
saveRDS(q.matrix5, file = "q.matrix5.rds")

# STRUCTURE-like barplot for the Q-matrix 
##CHANGE THE PALETTE

my.colors <- c("tomato", "orange", "lightblue", "wheat","olivedrab")
my.palette <- CreatePalette(my.colors, 5)

#library(RColorBrewer)
#my.colors <- brewer.pal(n = 5, name = "Dark2")
#my.palette <- CreatePalette(my.colors)

png(filename="plot2.png", width = 1800, height = 1000, units = "px")
par(mfrow = c(2, 2))
barplot(q.matrix2, border = NA, space = 0, sort.by.Q =TRUE,
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "K=2", 
        col.palette = my.palette) -> bp
axis(1, at = 1:nrow(q.matrix2), labels = bp$order, las = 3, cex.axis = 1) 

barplot(q.matrix3, border = NA, space = 0, sort.by.Q =TRUE,
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "K=3", 
        col.palette = my.palette) -> bp
axis(1, at = 1:nrow(q.matrix3), labels = bp$order, las = 3, cex.axis = 1)

barplot(q.matrix4, border = NA, space = 0, sort.by.Q =TRUE,
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "K=4", 
        col.palette = my.palette) -> bp
axis(1, at = 1:nrow(q.matrix4), labels = bp$order, las = 3, cex.axis = 1)

barplot(q.matrix5, border = NA, space = 0, sort.by.Q =TRUE,
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "K=5", 
        col.palette = my.palette) -> bp
axis(1, at = 1:nrow(q.matrix5), labels = bp$order, las = 3, cex.axis = 1)
dev.off()

###
#png(filename="plot3.png", width = 1000, height = 1300, units = "px")
#par(mfrow = c(2, 2))

#plot(q.matrix2, coord3, method = "map.max", interpol = FieldsKrigModel(10),  
#     main = "Ancestry model K=2",
#     xlab = "Longitude", ylab = "Latitude", 
#     resolution = c(300,300), cex = .5,
#     col.palette = my.palette)
#map('state', add = T, interior =T)
#map('lakes', add=TRUE, fill=TRUE, col='white', boundary='black')

#plot(q.matrix3, coord3, method = "map.max", interpol = FieldsKrigModel(10),  
#     main = "Ancestry model K=3",
#     xlab = "Longitude", ylab = "Latitude", 
#     resolution = c(300,300), cex = .5,
#     col.palette = my.palette)
#map('state', add = T, interior =T)
#map('lakes', add=TRUE, fill=TRUE, col='white', boundary='black')

#plot(q.matrix4, coord3, method = "map.max", interpol = FieldsKrigModel(10),  
#     main = "Ancestry model K=4",
#     xlab = "Longitude", ylab = "Latitude", 
#     resolution = c(300,300), cex = .5,
#     col.palette = my.palette)
#map('state', add = T, interior =T)
#map('lakes', add=TRUE, fill=TRUE, col='white', boundary='black')

#plot(q.matrix5, coord3, method = "map.max", interpol = FieldsKrigModel(10),  
#     main = "Ancestry model K=5",
#     xlab = "Longitude", ylab = "Latitude", 
#     resolution = c(300,300), cex = .5,
#     col.palette = my.palette)
#map('state', add = T, interior =T)
#map('lakes', add=TRUE, fill=TRUE, col='white', boundary='black')

#dev.off()

# retrieve tess3 results for different K
p.values2 <- pvalue(tess3.obj, K = 2)
p.values3 <- pvalue(tess3.obj, K = 3)
p.values4 <- pvalue(tess3.obj, K = 4)
p.values5 <- pvalue(tess3.obj, K = 5)

saveRDS(p.values2, file = "p.values2.rds")
saveRDS(p.values3, file = "p.values3.rds")
saveRDS(p.values4, file = "p.values4.rds")
saveRDS(p.values5, file = "p.values5.rds")

png(filename="plot4.png", width = 1000, height = 700, units = "px")
par(mfrow = c(2, 2))
hist(p.values2, col = "lightblue", main = "K=2") 
hist(p.values3, col = "lightblue", main = "K=3") 
hist(p.values4, col = "lightblue", main = "K=4") 
hist(p.values5, col = "lightblue", main = "K=5") 
dev.off()

##Fst calculations
# Benjamini-Hochberg algorithm
fdr.level = 1e-4
L = length(p.values2)

w2 = which(sort(p.values2) < fdr.level * (1:L)/L)
w3 = which(sort(p.values3) < fdr.level * (1:L)/L)
w4 = which(sort(p.values4) < fdr.level * (1:L)/L)
w5 = which(sort(p.values5) < fdr.level * (1:L)/L)

candidates2 = order(p.values2)[w2]
candidates3 = order(p.values3)[w3]
candidates4 = order(p.values4)[w4]
candidates5 = order(p.values5)[w5]

length(candidates2)
length(candidates3)
length(candidates4)
length(candidates5)

# Manhattan plot 
png(filename="plot5.png", width = 1000, height = 1300, units = "px")
par(mfrow = c(2, 2))

plot(p.values2, main = "Fst plot K=2", 
     xlab = "Locus id", ylab = "-log10(P-values)",
     cex = .3, col = "grey")
points(candidates2, -log10(p.values2)[candidates2], 
     pch = 19, cex = .2, col = "blue")

plot(p.values3, main = "Fst plot K=3", 
     xlab = "Locus id", ylab = "-log10(P-values)",
     cex = .3, col = "grey")
points(candidates3, -log10(p.values3)[candidates3], 
       pch = 19, cex = .2, col = "blue")

plot(p.values4, main = "Fst plot K=4", 
     xlab = "Locus id", ylab = "-log10(P-values)",
     cex = .3, col = "grey")
points(candidates4, -log10(p.values4)[candidates4], 
       pch = 19, cex = .2, col = "blue")

plot(p.values5, main = "Fst plot K=5", 
     xlab = "Locus id", ylab = "-log10(P-values)",
     cex = .3, col = "grey")
points(candidates5, -log10(p.values5)[candidates5], 
       pch = 19, cex = .2, col = "blue")
dev.off()