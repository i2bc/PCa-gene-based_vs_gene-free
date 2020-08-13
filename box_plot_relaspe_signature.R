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
dir.storeGene <- paste0("Result_article/Relapse/gene_based/")
dir.storeCtg <- paste0("Result_article/Relapse/gene_free/")

# File store data 
gene_discover_path = paste0( dir.storeGene, "data-sig-gene-tcga.tsv")
gene_valid_icgc = paste0( dir.storeGene, "data-sig-gene-icgc.tsv")
gene_valid_stelloo = paste0( dir.storeGene, "data-sig-gene-stelloo.tsv")

contig_discover_path = paste0( dir.storeCtg, "data-sig-contig-tcga.tsv")
contig_valid_icgc = paste0( dir.storeCtg, "data-sig-contig-icgc.tsv")
contig_valid_stelloo = paste0( dir.storeCtg, "data-sig-contig-stelloo.tsv")

# Load data
dataGeneDis <- as.data.frame(fread(gene_discover_path, sep="\t", header = TRUE))
dataGeneIcgc <- as.data.frame(fread(gene_valid_icgc, sep="\t", header = TRUE))
dataGeneStelloo <- as.data.frame(fread(gene_valid_stelloo, sep="\t", header = TRUE))

dataContigDis <- as.data.frame(fread(contig_discover_path, sep="\t", header = TRUE))
dataContigIcgc <- as.data.frame(fread(contig_valid_icgc, sep="\t", header = TRUE))
dataContigStelloo <- as.data.frame(fread(contig_valid_stelloo, sep="\t", header = TRUE))

###########################################

# Box plot for signature in TCGA
boxPlotGeneDis <- boxPlot(dataPlot = dataGeneDis, modeLevel = "gene" , dir.name = dir.storeGene, status = "tcga")
png(filename = paste0(dir.storeGene, "box-plot-sig-gene-tcga.png"), width = 3200, height = 1600, units = "px", res = 300)
plot(boxPlotGeneDis)
dev.off()

boxPlotContigDis <- boxPlot(dataPlot = dataContigDis, modeLevel = "contig" , dir.name = dir.storeCtg, status = "tcga")
png(filename = paste0(dir.storeCtg, "box-plot-sig-contig-tcga.png"), width = 3200, height = 1600, units = "px", res = 300)
plot(boxPlotContigDis)
dev.off()

# Box plot for signature in ICGC
boxPlotGeneIcgc <- boxPlot(dataPlot = dataGeneIcgc, modeLevel = "gene" , dir.name = dir.storeGene, status = "icgc")
png(filename = paste0(dir.storeGene, "box-plot-sig-gene-icgc.png"), width = 3200, height = 1600, units = "px", res = 300)
plot(boxPlotGeneIcgc)
dev.off()

boxPlotContigIcgc <- boxPlot(dataPlot = dataContigIcgc, modeLevel = "contig" , dir.name = dir.storeCtg, status = "icgc")
png(filename = paste0(dir.storeCtg, "box-plot-sig-contig-icgc.png"), width = 3200, height = 1600, units = "px", res = 300)
plot(boxPlotContigIcgc)
dev.off()

# Box plot for signature in STELLOO
boxPlotGeneStelloo <- boxPlot(dataPlot = dataGeneStelloo, modeLevel = "gene" , dir.name = dir.storeGene, status = "stelloo")
png(filename = paste0(dir.storeGene, "box-plot-sig-gene-stelloo.png"), width = 3200, height = 1600, units = "px", res = 300)
plot(boxPlotGeneStelloo)
dev.off()

boxPlotContigStelloo <- boxPlot(dataPlot = dataContigStelloo, modeLevel = "contig" , dir.name = dir.storeCtg, status = "stelloo")
png(filename = paste0(dir.storeCtg, "box-plot-sig-contig-stelloo.png"), width = 3200, height = 1600, units = "px", res = 300)
plot(boxPlotContigStelloo)
dev.off()
