################
# Author : Thi Ngoc Ha Nguyen
# Date   : 04/10/2020
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

topGene <- paste0(dir.discovery, "top500_TCGA_gene_norm.nb5.tsv")
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
dir.store <- paste0("Result_infer_signature/Relapse/gene_based")
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
  
  resSignROC <- takeDataReturnAUC(frameTrain = dataSigDis, frameTest = dataSigValid, absentContig = c(), status = "relapse")  
  
  resSignPR <- takeDataReturnPR(frameTrain = dataSigDis, frameTest = dataSigValid, absentContig = c(), status = "relapse")  
  
  return(list(resSignROC, dataSaveValid, resSignPR))
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
  
  #Store data frame valid
  name.file <- paste0(dir.store, "/data-sig-gene-icgc.tsv")
  write.table(resSigICGC[[2]], file = name.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

  finalTab <- rbind(resSigICGC[[1]][[1]], resSigICGC[[1]][[2]], resSigICGC[[1]][[3]], resSigICGC[[1]][[4]],
                    resSigICGC[[1]][[5]], resSigICGC[[1]][[6]], resSigICGC[[1]][[7]], resSigICGC[[1]][[8]],
                    resSigICGC[[1]][[9]], resSigICGC[[1]][[10]], resSigICGC[[1]][[11]], resSigICGC[[1]][[12]])
  rownames(finalTab) <- c( 
    "Mean ROC_AUC cv (TCGA)", "SD ROC_AUC cv (TCGA)", 
    "Mean ROC_AUC down cv (TCGA)", "SD ROC_AUC down cv (TCGA)",
    "Mean ROCAUC up cv (TCGA)", "SD ROC_AUC up cv (TCGA)",
    "ROC_AUC sig gene in ICGC", "PrAUC sig gene in ICGC",
    "ROC_AUC down sig gene in ICGC", "PrAUC down sig gene in ICGC",
    "ROC_AUC up sig gene in ICGC", "PrAUC up sig gene in ICGC")
  print(finalTab)

  name.file <- paste0(dir.store, "/Roc_auc-sig-gene-icgc.tsv")
  write.table(finalTab, file = name.file, quote = FALSE, row.name = TRUE, col.names = FALSE) 
  
  finalTab <- rbind(resSigICGC[[3]][[1]], resSigICGC[[3]][[2]], resSigICGC[[3]][[3]], resSigICGC[[3]][[4]],
                    resSigICGC[[3]][[5]], resSigICGC[[3]][[6]], resSigICGC[[3]][[7]], resSigICGC[[3]][[8]],
                    resSigICGC[[3]][[9]], resSigICGC[[3]][[10]], resSigICGC[[3]][[11]], resSigICGC[[3]][[12]])
  rownames(finalTab) <- c( 
    "Mean PR_AUC cv (TCGA)", "SD PR_AUC cv (TCGA)", 
    "Mean PR_AUC down cv (TCGA)", "SD PR_AUC down cv (TCGA)",
    "Mean PR_AUC up cv (TCGA)", "SD PR_AUC up cv (TCGA)",
    "ROC_AUC sig gene in ICGC", "PrAUC sig gene in ICGC",
    "ROC_AUC down sig gene in ICGC", "PrAUC down sig gene in ICGC",
    "ROC_AUC up sig gene in ICGC", "PrAUC up sig gene in ICGC")
  print(finalTab)  
  
  name.file <- paste0(dir.store, "/Pr_auc-sig-gene-icgc.tsv")
  write.table(finalTab, file = name.file, quote = FALSE, row.name = TRUE, col.names = FALSE)   
    
  #Valid in Stelloo
  resSigStelloo <- evaluate_gene(dataSigDis = dataSigTCGA, sigGeneDis = sigTCGA, 
                                 dataValidPath = geneStelloo, samplesConditionValidPath = sampleStelloo)
 
  #Store data frame valid
  name.file <- paste0(dir.store, "/data-sig-gene-stelloo.tsv")
  write.table(resSigStelloo[[2]], file = name.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  finalTab <- rbind(resSigStelloo[[1]][[1]], resSigStelloo[[1]][[2]], resSigStelloo[[1]][[3]], resSigStelloo[[1]][[4]],
                    resSigStelloo[[1]][[5]], resSigStelloo[[1]][[6]], resSigStelloo[[1]][[7]], resSigStelloo[[1]][[8]],
                    resSigStelloo[[1]][[9]], resSigStelloo[[1]][[10]], resSigStelloo[[1]][[11]], resSigStelloo[[1]][[12]])
  rownames(finalTab) <- c( 
    "Mean ROC_AUC cv (TCGA)", "SD ROC_AUC cv (TCGA)", 
    "Mean ROC_AUC down cv (TCGA)", "SD ROC_AUC down cv (TCGA)",
    "Mean ROCAUC up cv (TCGA)", "SD ROC_AUC up cv (TCGA)",
    "ROC_AUC sig gene in Stelloo", "PrAUC sig gene in Stelloo",
    "ROC_AUC down sig gene in Stelloo", "PrAUC down sig gene in Stelloo",
    "ROC_AUC up sig gene in Stelloo", "PrAUC up sig gene in Stelloo")
  print(finalTab)
  
  name.file <- paste0(dir.store, "/Roc_auc-sig-gene-stelloo.tsv")
  write.table(finalTab, file = name.file, quote = FALSE, row.name = TRUE, col.names = FALSE) 
  
  finalTab <- rbind(resSigStelloo[[3]][[1]], resSigStelloo[[3]][[2]], resSigStelloo[[3]][[3]], resSigStelloo[[3]][[4]],
                    resSigStelloo[[3]][[5]], resSigStelloo[[3]][[6]], resSigStelloo[[3]][[7]], resSigStelloo[[3]][[8]],
                    resSigStelloo[[3]][[9]], resSigStelloo[[3]][[10]], resSigStelloo[[3]][[11]], resSigStelloo[[3]][[12]])
  rownames(finalTab) <- c( 
    "Mean PR_AUC cv (TCGA)", "SD PR_AUC cv (TCGA)", 
    "Mean PR_AUC down cv (TCGA)", "SD PR_AUC down cv (TCGA)",
    "Mean PR_AUC up cv (TCGA)", "SD PR_AUC up cv (TCGA)",
    "ROC_AUC sig gene in Stelloo", "PrAUC sig contig in Stelloo",
    "ROC_AUC down sig gene in Stelloo", "PrAUC down sig gene in Stelloo",
    "ROC_AUC up sig gene in Stelloo", "PrAUC up sig gene in Stelloo")
  print(finalTab)  
  
  name.file <- paste0(dir.store, "/Pr_auc-sig-gene-stelloo.tsv")
  write.table(finalTab, file = name.file, quote = FALSE, row.name = TRUE, col.names = FALSE)   
  
}
#################################
#Run the script to infer
resGene <- pipeline(topProbesPath = topGene, samplesConditionDisPath = sampleTCGA, numruns = NUM_RUNS)
