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

# Discovery data using top 500 NB5 gene
topGene <-"Data_discovery/top2k_TCGA_gene_norm.nb5.out"
sampleTCGA <-"Data_discovery/sample_conditions.tsv"

# Validation data
sampleICGC <- "Data_validation/sample_conditions.tsv"
geneICGC <- "Data_validation/gene-counts-ICGC148-LRHR.norm.tsv"

NUM_RUNS=100
nKeep = 500

# Directory to store result
dir.store <- paste0("Result_infer_signature/gene_level")
dir.create(file.path(dir.store), showWarnings = FALSE, recursive = TRUE)

#################################
pipeline <- function(topProbesPath, samplesConditionDisPath, dataValidPath, samplesConditionValidPath, numruns){
  
  # loading top 500 genes/contigs in TCGA discovery dataset based on approach level
  countTopProbe <- as.data.frame(fread(topProbesPath, sep="\t", header = TRUE))
  
  countTopProbe <- data.frame(row.names = countTopProbe$feature, countTopProbe[,-c(1,2)], check.names=F)
    
  countTopProbe <- countTopProbe[1:nKeep,]
  
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
  
  resSign <- takeDataReturnAUC_ICGC(frameTrain = dataSigTCGA, frameTest = dataSigICGC, absentContig = c())
  
  print(paste("Performance of gene signature in TCGA:", resSign[[1]], "+/-", resSign[[2]]))
  print(paste("Performance of gene signature in ICGC:", resSign[[3]]))
  
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
