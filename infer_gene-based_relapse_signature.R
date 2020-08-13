################
# Author : Thi Ngoc Ha Nguyen
# Date   : 05/07/2020
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
dir.discovery <- "Data_discovery/Relapse/"

topGene <- paste0(dir.discovery, "top500_TCGA_gene_norm.nb5.out")
sampleTCGA <- paste0(dir.discovery, "sample_conditions.tsv")

# Validation ICGC data
dir.validation <- "Data_validation/Relapse/"
		     
sampleICGC <- paste0(dir.validation, "sample_cond_icgc.tsv")
geneICGC <- paste0(dir.validation, "gene-norm-ICGC-Relapse-2-2.tsv")

# Validation Stello data
sampleStelloo <- paste0(dir.validation, "sample_cond_stelloo.tsv")
geneStelloo <- paste0(dir.validation, "gene-filter-norm-Stello-Relapse-3-5.tsv")

NUM_RUNS=100

# Directory to store result
dir.store <- paste0("Result_article/Relapse/gene_based")
dir.create(file.path(dir.store), showWarnings = FALSE, recursive = TRUE)
############################################################################################
evaluate_gene <- function( sigGeneDis, dataSigDis, dataValidPath, samplesConditionValidPath){
  
  # dataframe of signature in validation dataset
  dataValid <- as.data.frame(fread(dataValidPath, sep="\t", header = TRUE))
  dataValid <- data.frame(row.names = dataValid$feature, dataValid[,-1], check.names=F)
  dataSignatureValid <- dataValid[sigGeneDis,]
  
  ########## processing of conditions ##############
  #create sample condition for signature
  samplesConditionValid<-as.data.frame(fread(samplesConditionValidPath,sep="\t", header= FALSE ,check.names=F))
  names(samplesConditionValid)<-c("Sample","condition")
  
  for (i in 1:dim(samplesConditionValid)[1]){
    if(samplesConditionValid$condition[i] == "no"){
      samplesConditionValid$condition[i] = "No"
    }
    if(samplesConditionValid$condition[i] == "yes"){
      samplesConditionValid$condition[i] = "Yes"
    }
  }
  
  #Transpose
  dataSignatureValid <- as.data.frame(t(dataSignatureValid))
  
  dataSignatureValid$Sample <- row.names(dataSignatureValid)
  
  # Mapping count table and condition in validation set
  dataValid<-merge(dataSignatureValid,samplesConditionValid, 
                   by.x="Sample",
                   by.y="Sample",
                   all.x=TRUE,
                   all.y=FALSE)
  dataValid$condition <- factor(dataValid$condition)
  
  # dataframe of signature in validation set
  dataSigValid <- dataValid[,c("condition",sigGeneDis)]  
  
  # From gene symbol to gene name
  signature.idsym <- fromGeneIDtakeGenName(sigGeneDis)
  signature.save <- signature.idsym$SYMBOL
  dataSaveValid <- dataValid[,c("Sample", "condition",sigGeneDis)]
  names(dataSaveValid) <- c("Sample", "condition",signature.save)
  
  names(dataSigValid) <- c("condition",signature.save)
  
  resSign <- takeDataReturnAUC(frameTrain = dataSigDis, frameTest = dataSigValid, absentContig = c(), status = "relapse")  
  
  return(list(resSign, dataSaveValid))
}

#################
pipeline <- function(topProbesPath, samplesConditionDisPath, numruns){
  
  # loading top  genes/contigs in TCGA discovery dataset based on approach level
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
  
  print(paste("List", length(signature.save), "genes in signature that found in TCGA:"))
  print(signature.save)
  
  name.file <- paste0(dir.store, "/sig-gene-tcga.txt")
  
  write.table(signature.idsym, file = name.file, quote = FALSE, row.names = FALSE, col.names = FALSE)

  # dataframe of signature in TCGA dataset
  dataSigTCGA <- dataTCGA[,c("condition", sigTCGA)]
  names(dataSigTCGA) <- c("condition", signature.save)

  name.file <- paste0(dir.store, "/data-sig-gene-tcga.tsv")
  dataSaveTCGA <- dataTCGA[,c("Sample","condition",sigTCGA)]
  colnames(dataSaveTCGA) <- c("Sample", "condition", signature.save)
  write.table(dataSaveTCGA, file = name.file, sep = "\t", quote = FALSE, row.names = FALSE)

  ########## finding signature in validation set ##############

  #Valid in ICGC
  resSigICGC <- evaluate_gene(dataSigDis = dataSigTCGA, sigGeneDis = sigTCGA, 
                              dataValidPath = geneICGC, samplesConditionValidPath = sampleICGC) 

  # Store important variable 
  name.file <- paste0(dir.store, "/import-sig-gene-icgc.tsv")
  importVar <- as.data.frame(resSigICGC[[1]][[4]]$importance)
  importVar$Gene <- rownames(resSigICGC[[1]][[4]]$importance)
  names(importVar) <- c("Overall", "Gene")
  dataStore <- importVar[round(order(importVar$Overall, decreasing = T),4),]
  print(dataStore[, c("Gene","Overall")])
  write.table(dataStore[, c("Gene", "Overall")], file = name.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

  #Store data frame valid
  name.file <- paste0(dir.store, "/data-sig-gene-icgc.tsv")
  write.table(resSigICGC[[2]], file = name.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  #Valid in Stelloo
  resSigStelloo <- evaluate_gene(dataSigDis = dataSigTCGA, sigGeneDis = sigTCGA, 
                                 dataValidPath = geneStelloo, samplesConditionValidPath = sampleStelloo)
  # Store important variable 
  name.file <- paste0(dir.store, "/import-sig-gene-stelloo.tsv")
  importVar <- as.data.frame(resSigStelloo[[1]][[4]]$importance)
  importVar$Gene <- rownames(resSigStelloo[[1]][[4]]$importance)
  names(importVar) <- c("Overall", "Gene")
  dataStore <- importVar[round(order(importVar$Overall, decreasing = T),4),]
  print(dataStore[, c("Gene","Overall")])
  write.table(dataStore[, c("Gene", "Overall")], file = name.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  #Store data frame valid
  name.file <- paste0(dir.store, "/data-sig-gene-stelloo.tsv")
  write.table(resSigStelloo[[2]], file = name.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  finalTab <- rbind(resSigICGC[[1]][[3]], resSigStelloo[[1]][[3]], resSigICGC[[1]][[1]], resSigICGC[[1]][[2]])
  rownames(finalTab) <- c("AUC signature gene in ICGC", "AUC signature gene in Stelloo", 
						  "Mean AUC cross-valid (TCGA)", "SD AUC cross-valid (TCGA)")
  print(finalTab)
  
  name.file <- paste0(dir.store, "/auc-sig-gene.tsv")
  write.table(finalTab, file = name.file, quote = FALSE, row.name = TRUE, col.names = FALSE)  
}
#################################
#Run the script to infer
resGene <- pipeline(topProbesPath = topGene, samplesConditionDisPath = sampleTCGA, numruns = NUM_RUNS)
