rm(list = ls())

library("ggtree")
library("ggplot2")
library("dendextend")

theKin <- read.table("kinship_centered_ibs.txt", header =F, row.names = 1)

clones<-row.names(theKin)
#write.csv(clones, file = "clones.csv")

#theKin <- read.csv("GAPIT.Genotype.Kin_Zhang.csv", header =F, row.names = 1)
distance.matrix <- dist(theKin, upper = TRUE)
hc <- hclust(distance.matrix) #, method = kinship.cluster)
#hc <- ape::as.phylo(hc)
plot(hc, type = "fan")

ggtree(hc, layout='circular')
info<- read.csv("info_nj_groups.csv", header =T)
#info$State<-as.factor(info$State)
#levels(info$State)
#tiff("dendrogram.png", width = 30, height = 25, res = 150, units = "cm")
ggtree(hc,layout='circular') + geom_tiplab(size=1.5, aes(angle=angle))
#dev.off()
#group by species
#hc2<-as.dendrogram(hc)
##plot(hc2, type = "fan")
#rect.dendrogram(hc2, k=3, lty = 5, lwd = 0, x=1, col=rgb(0.1, 0.2, 0.4, 0.1) ) 

str(info)

group<- list(
  S.integrifolium_E=info[info$K3_ext_species == "s_integrifolium_east",][,1],
  S.integrifolium_W=info[info$K3_ext_species == "s_integrifolium_west",][,1],
  S.integrifolium_NE=info[info$K3_ext_species == "s_integrifolium_northeast",][,1],
  S.integrifolium_x_perfoliatum=info[info$K3_ext_species == "s_integrifolium_x_perfoliatum",][,1],
  S.radula=info[info$K3_ext_species == "s_radula",][,1],
  S.perfoliatum=info[info$K3_ext_species == "s_perfoliatum",][,1],
  S.integrifolium_NA=info[info$K3_ext_species == "s_integrifolium_NA",][,1])

#group<- list(
#  EAST=info$Taxa[c(1:62)],
#  WEST=info$Taxa[c(63:142)],
#  N_A=info$Taxa[c(143:234)])


tiff("dendrogram.png", width = 30, height = 25, res = 150, units = "cm")

p_iris <- ggtree(hc, layout = 'fan', size=1.1) + geom_tiplab(size=1.5, aes(angle=angle))

#p2<- flip(p_iris, 124, 123)

p3<-groupOTU(p_iris, group, 'group') + aes(color=group) +
        theme(legend.position="right", legend.text = element_text(size=12), 
              legend.title = element_text(size=12),
              legend.key.size = unit(1, 'cm')) +
       scale_color_manual(values=c("#EC2D01","#D3D3D3","#E7B800","#6BAED6", "#9B5FE0", "#42F7C0", "#000000"))

p3

dev.off()