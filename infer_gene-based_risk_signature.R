################
# Author : Thi Ngoc Ha Nguyen
# Date   : 13/07/2020
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

# Discovery data using top 500 NB5 genes
dir.discovery <- "Data_discovery/Risk/"

topGene <- paste0(dir.discovery, "top500_TCGA_gene_norm.nb5.tsv")
sampleTCGA <- paste0(dir.discovery, "sample_conditions.tsv")

# Validation data
dir.validation <- "Data_validation/Risk/"		   
sampleICGC <- paste0(dir.validation, "sample_conditions.tsv")
geneICGC <- paste0(dir.validation,"gene-counts-ICGC148-LRHR.norm.tsv")

NUM_RUNS=100

# Directory to store result
dir.store <- paste0("Result_article/Risk/gene_based")
dir.create(file.path(dir.store), showWarnings = FALSE, recursive = TRUE)

#################################
pipeline <- function(topProbesPath, samplesConditionDisPath, dataValidPath, samplesConditionValidPath, numruns){
  
  # loading top 500 genes in TCGA discovery dataset
  countTopProbe <- as.data.frame(fread(topProbesPath, sep="\t", header = TRUE))
  
  countTopProbe <- data.frame(row.names = countTopProbe$feature, countTopProbe[,-c(1,2)], check.names=F)

  ########## processing of conditions ##############
  samplesConditionDis<-as.data.frame(fread(samplesConditionDisPath,sep="\t", header= FALSE ,check.names=F))
  names(samplesConditionDis)<-c("Sample","condition")
  
  #Transpose
  countTopProbe <- as.data.frame(t(countTopProbe))
  
  countTopProbe$Sample <- row.names(countTopProbe)
  
  # Mapping count table and condition in TCGA
  dataTCGA<-merge(countTopProbe,samplesConditionDis, 
                  by.x="Sample",
                  by.y="Sample",
                  all.x=TRUE,
                  all.y=FALSE)
  dataTCGA$condition <- factor(dataTCGA$condition)
  
  ########## processing of feature selection ##############
  
  # Finding the signature in TCGA discovery dataset
  dataModel <-dataTCGA[,-which(names(dataTCGA) %in% c("condition","Sample"))]
  dataModel$condition <-factor(dataTCGA$condition)
  
  # Set threshold for Lasso logistic regression
  th_lasso = 0.5
  
  sigTCGA <- extractSignatureStb(dataModel, thres = th_lasso)

  # From gene symbol to gene name
  signature.idsym <- fromGeneIDtakeGenName(sigTCGA)
  signature.save <- signature.idsym$SYMBOL
  
  print(paste("List", length(signature.save),"genes in signature:"))
  print(signature.save)
  
  name.file <- paste0(dir.store, "/sig-gene.tsv")
  
  write.table(signature.idsym, file = name.file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  # dataframe of signature in TCGA dataset
  dataSigTCGA <- dataTCGA[,c("condition",sigTCGA)]
  names(dataSigTCGA) <- c("condition", signature.save)

  ########## finding signature in validation set ##############
  # dataframe of signature in ICGC dataset
  dataValid <- as.data.frame(fread(dataValidPath, sep="\t", header = TRUE))
  dataValid <- data.frame(row.names = dataValid$feature, dataValid[,-1], check.names=F)
  dataSignatureValid <- dataValid[sigTCGA,]
  
  #create sample condition for signature
  samplesConditionValid<-as.data.frame(fread(samplesConditionValidPath,sep="\t", header= FALSE ,check.names=F))
  names(samplesConditionValid)<-c("Sample","condition")
  
  #Transpose
  dataSignatureValid <- as.data.frame(t(dataSignatureValid))
  
  dataSignatureValid$Sample <- row.names(dataSignatureValid)
  
  # Mapping count table and condition in ICGC
  dataICGC<-merge(dataSignatureValid,samplesConditionValid, 
                  by.x="Sample",
                  by.y="Sample",
                  all.x=TRUE,
                  all.y=FALSE)
  dataICGC$condition <- factor(dataICGC$condition)
  
  # dataframe of signature in ICGC dataset
  dataSigICGC <- dataICGC[,c("condition",sigTCGA)]  
  names(dataSigICGC) <- c("condition", signature.save)
  
  resSign <- takeDataReturnAUC(frameTrain = dataSigTCGA, frameTest = dataSigICGC, absentContig = c(), status = "risk")

  finalTab <- rbind(resSign[[3]], resSign[[1]], resSign[[2]])
  rownames(finalTab) <- c("AUC signature gene in ICGC", 
			 "Mean AUC cross-valid (TCGA)", "SD AUC cross-valid (TCGA)")
  print(finalTab)
  
  name.file <- paste0(dir.store, "/auc-sig-gene.tsv")
  write.table(finalTab, file = name.file, quote = FALSE, row.name = TRUE, col.names = FALSE)   
  
  # Store important variable 
  name.file <- paste0(dir.store, "/import-sig-gene.tsv")
  importVar <- as.data.frame(resSign[[4]]$importance)
  importVar$Gene <- rownames(resSign[[4]]$importance)
  names(importVar) <- c("Overall", "Gene")
  dataStore <- importVar[round(order(importVar$Overall, decreasing = T),4),]
  print("Contribution of each gene in signature:")
  print(dataStore[, c("Gene","Overall")])
  write.table(dataStore[, c("Gene", "Overall")], file = name.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

  name.file <- paste0(dir.store, "/data-sig-gene-tcga.tsv")
  dataSaveTCGA <- dataTCGA[,c("Sample","condition",sigTCGA)]
  colnames(dataSaveTCGA) <- c("Sample", "condition", signature.save)
  write.table(dataSaveTCGA, file = name.file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  name.file <- paste0(dir.store, "/data-sig-gene-icgc.tsv")
  dataSaveICGC <- dataICGC[,c("Sample","condition",sigTCGA)]
  colnames(dataSaveICGC) <- c("Sample", "condition", signature.save)
  write.table(dataSaveICGC, file = name.file, sep = "\t", quote = FALSE, row.names = FALSE)
  
}
#################################
resGene <- pipeline(topProbesPath = topGene, samplesConditionDisPath = sampleTCGA, dataValidPath = geneICGC,
                    samplesConditionValidPath = sampleICGC, numruns = NUM_RUNS)
