############################
#PCA PLOT
############################
library("ggplot2")
library("plotly")
library("scales")

dat<-read.csv("pca_data.csv", header =T)

plot_ly(x=dat$PC1, y=dat$PC2, z=dat$PC3, type="scatter3d", mode="markers", color=dat$K3_ext_species)

p1<-ggplot(dat,aes(x=PC1,y=PC2,color=K3_ext_species, label = row.names(dat))) +
    geom_point(size = 4) +
    theme( plot.title = element_text(size=18, hjust = 0.5), 
           panel.background = element_rect(fill = "white", colour = "grey50"),
           axis.text = element_text(face="bold", size=14),
           axis.title = element_text(face="bold", size=14),
           legend.text = element_text( size=12),
           legend.title = element_text( size=12 )) +
  scale_color_manual(
    values = c(
      "s_integrifolium_east"  = "#EC2D01",
      "s_integrifolium_west" = "#6BAED6",
      "s_integrifolium_northeast"  = "#E7B800",
      "s_integrifolium_x_perfoliatum" = "#9B5FE0",
      "s_radula" = "#000000",
      "s_perfoliatum"  = "#3BB143",
      "s_integrifolium_NA" = "#dcdcdb"
      ),
    name = "Groups"   # legend title
  )
p2<- p1 + xlab("PC1 (8.45%)") + ylab("PC2 (2.06%)")
p2
ggsave("pca_plot_k3_groups_pc1-2.png", plot=p2, height=5, width=8, dpi=300)

##longitude plot PC1 vs PC2
p3<-ggplot(dat,aes(x=PC1,y=PC2)) +
  geom_point(aes(color = long), size = 4) +
  geom_point(
    data = subset(dat, is.na(long)),
    shape = 21, fill = "white", colour = "black", size = 4, stroke = 0.3
  ) +
  theme( plot.title = element_text(size=18, hjust = 0.5), 
         panel.background = element_rect(fill = "white", colour = "grey50"),
         axis.text = element_text(face="bold", size=14),
         axis.title = element_text(face="bold", size=14),
         legend.text = element_text( size=12),
         legend.title = element_text( size=12 )) +
  scale_colour_viridis_c(begin = 0, end = 1, direction = -1, na.value = NA, option = "inferno", name = "Longitude") + 
  xlab("PC1 (8.45%)") + ylab("PC2 (2.06%)")
p3 
ggsave("pca_plot_long_pc1-2.png", plot=p3, height=5, width=6, dpi=300)

dev.off()

##longitude plot PC1 vs PC3
p3<-ggplot(dat,aes(x=PC1,y=PC3)) +
  geom_point(aes(color = long), size = 4) +
  geom_point(
    data = subset(dat, is.na(long)),
    shape = 21, fill = "white", colour = "black", size = 4, stroke = 0.3
  ) +
  theme( plot.title = element_text(size=18, hjust = 0.5), 
         panel.background = element_rect(fill = "white", colour = "grey50"),
         axis.text = element_text(face="bold", size=14),
         axis.title = element_text(face="bold", size=14),
         legend.text = element_text(face="bold", size=12),
         legend.title = element_text(face="bold", size=12 )) +
  scale_colour_viridis_c(begin = 0, end = 1, direction = -1, na.value = NA, option = "inferno") + 
  xlab("PC1 (8.45%)") + ylab("PC3 (1.79%)")
p3 
ggsave("pca_plot_long_pc1-3_new.png", plot=p3, height=5, width=6, dpi=300)

dev.off()


