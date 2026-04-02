library(devtools)
#devtools::install_github("jtlovell/GENESPACE", upgrade = F)
library(GENESPACE)

library(data.table)
library(Biostrings)
library(ggplot2)
library(IRanges)
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)

wd <- "/cluster/home/rsouza/genespace/"
if(!dir.exists(wd))
  dir.create(wd)

faFile <- file.path(wd, "genomes/Silphium_integrifolium_var_Bad_Astra.mainGenome.fasta")
geneGffFile <- file.path(wd, "genes/sintegrifolium.gene.gff3")
repBedFile <- file.path(wd, "repeats/all.int.fa.mod.EDTA.RM.bed")

dnaSS <- Biostrings::readDNAStringSet(faFile)
names(dnaSS) <- sapply(names(dnaSS), function(x) strsplit(x, " ")[[1]][1])
dnaSS <- dnaSS[Biostrings::width(dnaSS) > 1e6]
seqInfo <- pull_seqInfo(dnaSS)

#x<-as.data.frame(seqInfo)
#ADAPTATION FOR RUNNING
#seqInfo@is_circular<- c(F,F,F,F,F,F,F)
#seqInfo@genome<-c("A","B","C","D","E","F","G")

genes <- as.data.frame(rtracklayer::readGFF(
  geneGffFile,
  columns = c("seqid", "start", "end", "type"), 
  tags = "Parent", 
  filter = list(type = c("gene"))))
genes <- subset(genes, seqid %in% names(dnaSS))

library(data.table)
repeats <- fread(
  repBedFile, 
  select = c(1:3, 12),
  col.names = c("chr", "start", "end", "type"))
repeats <- subset(repeats, chr %in% names(dnaSS))

bedList <- list(
  gene = with(subset(genes, type == "gene"), data.frame(
    chr = seqid, start = start, end = end)),
  #exon = with(subset(genes, type == "exon"), data.frame(
  #  chr = seqid, start = start, end = end)),
  #intron = with(subset(genes, type == "intron"), data.frame(
  #  chr = seqid, start = start, end = end)),
  Helitron = with(subset(repeats, grepl("Helitron", type)), data.frame(
    chr = chr, start = start, end = end)),
  TIR = with(subset(repeats, grepl("TIR", type)), data.frame(
    chr = chr, start = start, end = end)),
  LTR = with(subset(repeats, grepl("LTR", type)), data.frame(
    chr = chr, start = start, end = end)))
  #otherRepeat = with(subset(repeats, !grepl("Helitron|TIR|LTR", type)), data.frame(
   # chr = chr, start = start, end = end)),
  
print(lapply(bedList, head))

summary(repeats)

suppressWarnings(genomeClasses <- classify_genome(
  dnaSS = dnaSS, listOfBeds = bedList, verbose = T))

print(lapply(genomeClasses, head))

sws <- slide_genome(
  seqInfo = seqInfo,
  listOfGrs = genomeClasses[c(1,5,2,3,4)],
  windowSize = 5e6,
  stepSize = 1e6)
head(sws)

library(ggplot2)
plotCols <- c("#CC6828", "#FFFFFF", "#0F4F8B", "#4C86C6", "#AED8E6")
plotTheme <- theme(
  panel.background = element_rect(fill = "black"),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  panel.spacing = unit(0, "cm"),
  strip.background = element_blank())

sws$id2<- factor(sws$id, levels = c("gene", "unknown","Helitron","TIR","LTR"))

ggplot(sws, aes(x = start/1e6, y = propWind, fill = id2))+
  geom_area()+
  scale_fill_manual(values = plotCols) +
  facet_grid(.~ gsub("chr","",chr), scale = "free", space = "free") +
  scale_x_continuous(expand = c(0.05,0.05))+
  scale_y_continuous(expand = c(0.005,0.005))+
  plotTheme +
  labs(x = "Chromosomes (scaled by physical size; 5Mb windows, 1Mb steps)",
       y = "Cumulative % of Sequence") +
theme(axis.text.x = element_text(face="bold",
                                   size=10),
        axis.text.y = element_text(face="bold",
                                   size=10))

ggsave("plot_new2.png", dpi = 150, width = 15, height = 6)

#p2 <- ggplot(sws,
#            aes(x = start/1e6, y = propWind, fill = id))+
#  geom_area()+
#  scale_fill_manual(values = plotCols) +
#  facet_grid(.~gsub("chr","",chr), scale = "free", space = "free") +
#  scale_x_continuous(expand = c(0.05,0.05))+
#  scale_y_continuous(expand = c(0.005,0.005))+
#  plotTheme +
#  labs(x = "Chromosomes (scaled by physical size; 5Mb windows, 1Mb steps)",
#       y = "Cumulative % of Sequence")

#pdf(file = "/cluster/home/rsouza/genespace/Plot_2.pdf",   # The directory you want to save the file in
#    width = 4, # The width of the plot in inches
#    height = 4)

#gridExtra::grid.arrange(p1, nrow = 1)

#dev.off()