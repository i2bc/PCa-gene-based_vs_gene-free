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
library(ggplot2)
################
#######################################################################
## scripts
#######################################################################

# Directory to store result
dir.storeGene <- "Result_article/Risk/gene_level"
dir.storeCtg <- "Result_article/Risk/kmer_level"

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
boxPlot <- function(dataPlot, modeLevel, dir.name, status){
  # Plot box plot 
  #
  # Args:
  #   dataPlot: Data frame containes expression of selected probes
  #   modeLevel : String indicates gene-level or kmer-level
  #   dir.name : Path storing
  #   status : String indicates with database will be use (tcga or icgc)
  # Results:
  #   Figure for boxplot
  
  #Plot boxplot for the count of signature
  white_background<-theme(axis.line = element_line(colour = "black"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank())
  
  if (modeLevel == "contig"){
    for (i in 1: length(names(dataPlot))){
      if (names(dataPlot)[i] == "COL1A1"){
        names(dataPlot)[i] = "COL1A1.a"
      }
      if(names(dataPlot)[i] == "COL1A1.v20"){
        names(dataPlot)[i] = "COL1A1.b"
      }
      if(names(dataPlot)[i] =="COL1A1.v5"){
        names(dataPlot)[i] = "COL1A1.c"
      }
      if(names(dataPlot)[i]=="ASPN.v24"){
        names(dataPlot)[i] = "ASPN.a"
      }
      if(names(dataPlot)[i]=="ASPN"){
        names(dataPlot)[i] = "ASPN.b"
      }
    }
  }
  #Unlog the dataPlot into CPM/CPB expression
  unlog <- exp(dataPlot[-c(1,2)])
  dataNormalized <- cbind(dataPlot[c(1,2)], unlog)
  
  nameProbe <- names(dataNormalized[3:length(names(dataNormalized))])
  my.data <- melt(dataNormalized,measure.vars=c(names(dataNormalized[3:length(names(dataNormalized))])))
  my.data$value <- as.numeric(my.data$value)
  if(modeLevel == "contig"){xlabel = "Contigs"}
  if(modeLevel == "gene"){xlabel = "Genes"}
  
  boxplots_frame<-data.frame()
  
  for(i in 1: length(nameProbe)){
    
    HRvalues <- my.data[which(my.data[,"condition"]=="HR" & my.data$variable == nameProbe[i]),]$value
    
    meanHR <- mean(HRvalues)
    
    LRvalues <- my.data[which(my.data[,"condition"]=="LR" & my.data$variable == nameProbe[i]),]$value
    
    meanLR<- mean(LRvalues)
    
    logFC<- log(meanHR/ meanLR)
    
    wilcoxon_pvalue<-as.numeric(format(wilcox.test(LRvalues,HRvalues)$p.value,scientific=T,digits=4,quote=F))
    
    boxplots_frame<-rbind(boxplots_frame,data.frame(probes=nameProbe[i],mean_LR = round(meanLR, digits=2),mean_HR = round(meanHR, digits=2),logFC=round(logFC, digits=2),wilcoxon_pvalue=wilcoxon_pvalue)) 
    
  }
  
  # The order of contig based on logFC decreasing
  ctgLogFC <- nameProbe[order(abs(boxplots_frame$logFC), decreasing = TRUE)]
  
  frameLogFC <- boxplots_frame[match(ctgLogFC,boxplots_frame$probes),]
  # Save
  
  if (status == "tcga"){
    name.file <- paste0(dir.name, "/logFC-sig-", modeLevel, "-tcga.tsv")
  }else {
    name.file <- paste0(dir.name, "/logFC-sig-", modeLevel, "-icgc.tsv")
  }
  
  write.table(frameLogFC, file = name.file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Reoder the dataframe
  dataPlot$condition = factor(dataPlot$condition, levels = c("LR", "HR"))
  data.plot <- melt(dataPlot,measure.vars= ctgLogFC)
  
  boxPlot<-ggplot(data.plot, aes(x=variable,y=value,fill=condition))+ theme_bw() + geom_boxplot(outlier.size =0.5)+
    ylab("log(normalized expression)")+xlab(xlabel)+ggtitle(paste("Distribution of probes expression\n","nb probes = ",
                                                                  length(unique(my.data$variable)),"\nnb samples = ",length(unique(my.data$Sample)),
                                                                  " (LR = ",length(unique(my.data[which(my.data$condition=="LR"),]$Sample))," ; HR = ",length(unique(my.data[which(my.data$condition=="HR"),]$Sample)),")",sep=""))+
    theme(axis.text.x=element_text(color="black",size=12,angle=45,hjust=1),axis.text.y = element_text(color="black",size=14),plot.title = element_text(hjust = 0.5,size=16),
          axis.title=element_text(size=16,face="bold"),legend.text=element_text(size=14),legend.title=element_text(size=16))+
    scale_fill_manual(values=c("LR"="cornflowerblue","HR"="red"),name="conditions") + white_background
}
##################################
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
