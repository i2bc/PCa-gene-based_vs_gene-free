################
# Author : Thi Ngoc Ha Nguyen
# Date   : 4/06/2020
# Email  : thi-ngoc-ha.nguyen@i2bc.paris-saclay.fr
################

rm(list=ls())
set.seed(12678)

################
# Load libraries
################

library(data.table)
################
# Load list of functions 
################
source("useful_functions.R")
#######################################################################
## scripts
#######################################################################

# Directory to store result
dir.storeGene <- paste0("Result_artical/gene_level")
dir.create(file.path(dir.storeGene), showWarnings = FALSE, recursive = TRUE)

dir.storeCtg <- paste0("Result_artical/kmer_level")
dir.create(file.path(dir.storeCtg), showWarnings = FALSE, recursive = TRUE)

# File store data LR vs HR
gene_discover_path = paste0( dir.storeGene, "/data-sig-gene-tcga.tsv")
gene_valid_path = paste0( dir.storeGene, "/data-sig-gene-icgc.tsv")

contig_discover_path = paste0( dir.storeCtg, "/data-sig-contig-tcga.tsv")
contig_valid_path = paste0( dir.storeCtg, "/data-sig-contig-icgc.tsv")

# Load data
dataGeneDis <- as.data.frame(fread(gene_discover_path, sep="\t", header = TRUE))
dataGeneValid <- as.data.frame(fread(gene_valid_path, sep="\t", header = TRUE))

dataContigDis <- as.data.frame(fread(contig_discover_path, sep="\t", header = TRUE))
dataContigValid <- as.data.frame(fread(contig_valid_path, sep="\t", header = TRUE))

###########################################

# Box plot for signature in TCGA
boxPlotGeneDis <- boxPlot(dataPlot = dataGeneDis, modeLevel = "gene" , dir.name = dir.storeGene, status = "tcga")
png(filename = paste0(dir.storeGene, "/box-plot-sig-gene-tcga.png"), width = 3200, height = 1600, units = "px", res = 300)
plot(boxPlotGeneDis)
dev.off()

boxPlotContigDis <- boxPlot(dataPlot = dataContigDis, modeLevel = "contig" , dir.name = dir.storeCtg, status = "tcga")
png(filename = paste0(dir.storeCtg, "/box-plot-sig-contig-tcga.png"), width = 3200, height = 1600, units = "px", res = 300)
plot(boxPlotContigDis)
dev.off()

# Box plot for signature in ICGC
boxPlotGeneValid <- boxPlot(dataPlot = dataGeneValid, modeLevel = "gene" , dir.name = dir.storeGene, status = "icgc")
png(filename = paste0(dir.storeGene, "/box-plot-sig-gene-icgc.png"), width = 3200, height = 1600, units = "px", res = 300)
plot(boxPlotGeneValid)
dev.off()

boxPlotContigValid <- boxPlot(dataPlot = dataContigValid, modeLevel = "contig" , dir.name = dir.storeCtg, status = "icgc")
png(filename = paste0(dir.storeCtg, "/box-plot-sig-contig-icgc.png"), width = 3200, height = 1600, units = "px", res = 300)
plot(boxPlotContigValid)
dev.off()
