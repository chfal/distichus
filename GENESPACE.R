## Activate R studio from the terminal after you've activated orthofinder conda env

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("Biostrings", "rtracklayer"))

#a

#if (!requireNamespace("devtools", quietly = TRUE))
#  install.packages("devtools")
#devtools::install_github("jtlovell/GENESPACE")

library(GENESPACE)
library(ggplot2)

###############################################
# -- change paths to those valid on your system
genomeRepo <- "/Users/jmt410/Documents/Projects/Anolis/genespace/genomeRepo"
wd <- "/Users/jmt410/Documents/Projects/Anolis/genespace"
path2mcscanx <- "/Users/jmt410/Documents/Projects/Anolis/genespace/MCScanX-master/"
###############################################

# -- download raw data from NCBI for human and chicken genomes
dir.create(genomeRepo)
rawFiles <- download_exampleData(filepath = genomeRepo)

# -- parse the annotations to fastas with headers that match a gene bed file
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = c("human", "chicken"),
  genomeIDs = c("human", "chicken"),
  presets = "ncbi",
  genespaceWd = wd)

# -- initalize the run and QC the inputs
gpar <- init_genespace(
  wd = wd, 
  path2mcscanx = path2mcscanx)

# -- accomplish the run
out <- run_genespace(gpar)

# -- Specifying a set of regions of interest
#Here, we’ll do the query on a single bed file-like data.frame object containing three regions of interest:
  
roi <- data.frame(
    genome = c("human", "human", "chicken"),
    chr = c("6", "X", "Z"),
    start = c(0, 0, 0),    # c(3.8e6, 0, 0) #Original tutorial values
    end = c(Inf, Inf, Inf))  # c(4.5e6, Inf, 1e6) # Orginial tutorial values 

# -- Pulling hits against a reference for a specific region
# For each line in the regions of interest bed parameters, we extract all hits that map to it. 
# We can either pull ALL hits or just those in synteny. 
# Since synteny is so strong in this set, its reasonable to look only at syntenic hits

qreturn <- query_hits(gsParam = out, bed = roi, synOnly = TRUE)

# -- Making dotplots from the hits
# We offer gghits as a method to make a dotplot from hits stored in memory as a data.table. 
# We can do this here for one of the regions of interest: just chicken and human X chr:
  
xvx <- subset(qreturn[["human, X: 0-Inf"]], genome2 == "chicken")
gghits(xvx, useOrder = F)

# -- Querying pan-genes
# A set of pan-genes are in the same orthogroup and can be placed into a synteny position against a reference genome. 
# We can extract these as above using query_pangenes:
  
test <- query_pangenes(
  gsParam = out, bed = roi)

# -- riparian plots for regions of interest
# See the riparian plotting guide for more details, 
# but we can add a color column to the bed and only look at these regions, 
# either with the full synteny map as a background:
  
roibed <- roi[,c("genome", "chr")]
roibed$color <- c("gold", "purple", "green")
ripd <- plot_riparian(
  gsParam = out, 
  useRegions = FALSE, 
  highlightBed = roibed)

# -- just looking at the focal regions of interest:
  
ripd <- plot_riparian(
    gsParam = out,  
    useRegions = FALSE,
    highlightBed = roibed, 
    backgroundColor = NULL)

#### Actually trying this on Anolis data 

###############################################
# -- change paths to those valid on your system
genomeRepo <- "/Users/jmt410/Documents/Projects/Anolis/genespace/genomeRepo"
wd <- "/Users/jmt410/Documents/Projects/Anolis/genespace"
path2mcscanx <- "/Users/jmt410/Documents/Projects/Anolis/genespace/MCScanX-master/"
###############################################

# -- download raw data from NCBI for human and chicken genomes
#dir.create(genomeRepo)
#rawFiles <- download_exampleData(filepath = genomeRepo)

list.files(genomeRepo)

# JT - I've already created the .bed and .fa peptide files so this step wasnt needed
# -- parse the annotations to fastas with headers that match a gene bed file
#parsedPaths <- parse_annotations(
#  rawGenomeRepo = genomeRepo,
#  genomeDirs = c("sagrei"),
#  genomeIDs = c("sagrei"),
#  presets = "ncbi",
#  genespaceWd = wd)

# -- initalize the run and QC the inputs
gpar <- init_genespace(
  wd = wd, 
  path2mcscanx = path2mcscanx)

# -- accomplish the run
out <- run_genespace(gpar)

# -- Reload genespace parameter object 
load('/Users/jmt410/Documents/Projects/Anolis/genespace/results/gsParams.rda', verbose = TRUE)

# -- Specifying a set of regions of interest
#Here, we’ll do the query on a single bed file-like data.frame object containing three regions of interest:

roi <- data.frame(
  genome = c("distichus", "distichus", "carolinensis","carolinensis","carolinensis","carolinensis","carolinensis","carolinensis","carolinensis","carolinensis","carolinensis","carolinensis","carolinensis","carolinensis","carolinensis","carolinensis","carolinensis","carolinensis","carolinensis","carolinensis","carolinensis","carolinensis","carolinensis","carolinensis","carolinensis","carolinensis","carolinensis"),
  chr = c("8", "12", "chrUn0010","chrUn0273","chrUn0184","chrUn0259","chrUn0407","chrUn0336","chrUn0494","LGb","chrUn0172","chrUn0090","chrUn0146","chrUn0225","chrUn0231","chrUn0090","chrUn0756","chrUn0135","chrUn0193","chrUn0071","chrUn0426","chrUn0315","chrUn0210","chrUn0499","chrUn0457","chrUn0420","chrUn0664"),
  start = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),    # c(3.8e6, 0, 0) #Original tutorial values
  end = c(Inf, Inf, Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf))  # c(4.5e6, Inf, 1e6) # Orginial tutorial values 


# -- Pulling hits against a reference for a specific region
# For each line in the regions of interest bed parameters, we extract all hits that map to it. 
# We can either pull ALL hits or just those in synteny. 
# Since synteny is so strong in this set, its reasonable to look only at syntenic hits

qreturn <- query_hits(gsParam = out, bed = roi, synOnly = FALSE)
#qreturn <- query_hits(gsParam = gsParam, bed = roi, synOnly = TRUE) # JT I used this for plotting without rerunning GENESPACE

# -- Making dotplots from the hits
# We offer gghits as a method to make a dotplot from hits stored in memory as a data.table. 
# We can do this here for one of the regions of interest: just chicken and human X chr:
# EXAMPLE: 
#xvx <- subset(qreturn[["human, X: 0-Inf"]], genome2 == "chicken")
#gghits(xvx, useOrder = F)

# -- Querying pan-genes
# A set of pan-genes are in the same orthogroup and can be placed into a synteny position against a reference genome. 
# We can extract these as above using query_pangenes:

test <- query_pangenes(
  gsParam = out, bed = roi)

# -- riparian plots for regions of interest
# See the riparian plotting guide for more details, 
# but we can add a color column to the bed and only look at these regions, 
# either with the full synteny map as a background:

roibed <- roi[,c("genome", "chr")]

#Specify colours for each chromosome - mostly used for plotting specific chromosomes only
roibed$color <- c("#B66DFF", "#FFD700", "#009292", "#006DDB", "#F13030", "#d6a184", "#CF1259", "#83E8BA", "#F4E76E", "#8D3B72", "#B5E2FA")

# Colours for distichus 
roibed$color <- c("#EE2A7B", "#A4509F","#F34791", "#D8266E", "#F06C9F", "#C32060", 
                  "#FA548B", "#BA1A53", "#E13D86", "#F872A8", "#C91054", "#F55C94",  
                  "#E84682", "#D41E67", "#F96A9C", "#BF1C56", "#F03F85", "#B366B2", 
                  "#8E3D8C", "#BC7AC0", "#7D3580", "#9B4893", "#C38AC8", "#732B6B", 
                  "#D09FD4", "#863F89", "#AA5BA7")
         

## This is the all against all plot
ripd <- plot_riparian(
  gsParam = gsParam,
  refGenome = "distichus",
  useRegions = T,
  minChrLen2plot = 100
 )

# -- just looking at the focal regions of interest:
# The details for the extra parameters to customize plot can be found here: 
# https://htmlpreview.github.io/?https://github.com/jtlovell/tutorials/blob/main/riparianGuide.html

ripd <- plot_riparian(
  gsParam = gsParam,
  #genomeIDs = gsParam$genomeIDs,
  genomeIDs = c("carolinensis", "distichus", "sagrei"), # This sets the order of the genomes from the bottom, up
  refGenome = "sagrei", # Sets the reference genome - I use this as an achor to construct the figure - not in a biologically meaningful way
  customRefChrOrder = c("11","7"), # Set order for reference chromosomes
  forceRecalcBlocks = TRUE,
  syntenyWeight = 1,
  minChrLen2plot = 0,  # adjusts default chromosome size to 0 to plot small chr's
  braidAlpha = 0.8,  
  chrLabFontSize = 8,
  chrExpand = 0.1,
  chrBorderLwd = 0.5,
  useRegions = TRUE,
  highlightBed = roibed,
  backgroundColor = NULL,
  chrBorderCol = "grey",
  addThemes = theme_minimal(   # From ggplot 
    base_line_size = 0), 
  )
