library(GENESPACE)

###############################################
# -- change paths to those valid on your system
#genomeRepo <- "/cluster/home/rsouza/genespace/riparian_run2/genomeRepo"
wd <- "/cluster/home/rsouza/genespace/riparian_run2"
path2mcscanx <- "/cluster/home/rsouza/apps/MCScanX-master/"
path2orthofinder <-"/cluster/home/rsouza/apps/OrthoFinder/orthofinder"
path2diamond <- "/cluster/home/rsouza/apps/OrthoFinder/bin/diamond"
###############################################

# -- download raw data from NCBI for human and chicken genomes
#dir.create(genomeRepo)
#rawFiles <- download_exampleData(filepath = genomeRepo)

# -- parse the annotations to fastas with headers that match a gene bed file

#parsedPaths <- parse_annotations(
#  rawGenomeRepo = genomeRepo,
#  genomeDirs = c("carrot", "lettuce", "mimulus", "olives", "sunflower"),
#  genomeIDs = c("carrot", "lettuce", "mimulus", "olives", "sunflower"),
#  presets = "phytozome",
#  genespaceWd = wd)

# -- initalize the run and QC the inputs
gpar <- init_genespace(
  wd = wd, 
  path2mcscanx = path2mcscanx,
  path2orthofinder = path2orthofinder,
  path2diamond = path2diamond)

# -- accomplish the run
out <- run_genespace(gpar)

##
ggthemes<-ggplot2::theme(panel.background = ggplot2::element_rect(fill="white"))

ripDat <- plot_riparian(
  gsParam = out, 
  chrFill = "lightgrey",
  refGenome = "S.integrifolium", 
  genomeIDs = c("S.taccada","H.annuus_xrq", "H.annuus_ha412","S.perfoliatum","S.integrifolium"),
  addThemes = ggthemes,
  forceRecalcBlocks = FALSE)
