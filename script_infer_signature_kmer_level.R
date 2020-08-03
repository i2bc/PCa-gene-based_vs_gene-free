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

library(knitr)
library(data.table)
library(RColorBrewer)
library(edgeR)
################
# Load list of functions 
################
source("useful_functions.R")
#######################################################################
## scripts
#######################################################################

# Discovery data kmer level
topContig <- "Data_discovery/top2k_contig_merge_norm.nb5.out"
sampleTCGA <-"Data_discovery/sample_conditions.tsv"

# Validation data
sampleICGC <- "Data_validation/sample_conditions.tsv"
kmerICGC <- "Data_validation/raw-counts.tsv.gz"
totalKmers <-"Data_validation/sum_counts.tsv"


NUM_RUNS=100
nKeep = 500

# Directory to store result
dir.store <- paste0("Result_article/kmer_level")
dir.create(file.path(dir.store), showWarnings = FALSE, recursive = TRUE)
############################################################################################
normalizeContig <-function(validCountPath, libSizePath){
  
  # Normalization based on logCPM
  countKmerICGC <- as.data.frame(fread(validCountPath, sep="\t", header = TRUE))
  
  rownames(countKmerICGC) <- countKmerICGC$tag
  
  countKmerICGC <- countKmerICGC[,-1]
  
  # Normalization based on logCPM
  inforKmerICGC <- as.data.frame(fread(libSizePath, sep="\t", header = FALSE))
  names(inforKmerICGC) <- c("Sample", "Total_kmers")
  rownames(inforKmerICGC) <- inforKmerICGC$Sample
  
  #Sort libSize base on sample name
  inforKmerICGC <- inforKmerICGC[colnames(countKmerICGC),]
  libSizeICGC <- inforKmerICGC$Total_kmers
  
  # compute logBPM
  logCpmICGC <- log(countKmerICGC/expandAsMatrix(libSizeICGC/1e+09, dim(countKmerICGC))+1)
  
  return (logCpmICGC)
}
#################################
pipeline <- function(topProbesPath, samplesConditionDisPath, dataValidPath, samplesConditionValidPath, numruns){
  
  # loading top 500 contigs in TCGA discovery dataset
  countTopProbe <- as.data.frame(fread(topProbesPath, sep="\t", header = TRUE))
    
  countTopProbe <- data.frame(row.names = countTopProbe$feature, countTopProbe[,-c(1:4)], check.names=F)    
    
  countTopProbe <- countTopProbe[1:nKeep,]

  ########## processing of conditions ##############
  samplesConditionDis<-as.data.frame(fread(samplesConditionDisPath,sep="\t", header= FALSE ,check.names=F))
  names(samplesConditionDis)<-c("Sample","condition")
  
  #Transpose
  countTopProbe <- as.data.frame(t(countTopProbe))
  
  #Change the double value in right format
  for (i in 1:dim(countTopProbe)[2]){
    countTopProbe[,i] <- as.double(as.character(countTopProbe[,i]))
  }

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
  
  # Mapping contig to genenome
  # Create fa file
  insertStr <- paste0(">ctg_",seq(1:length(sigTCGA)))
  preFA <- as.vector(t(cbind(insertStr, matrix(sigTCGA, length(sigTCGA), byrow=T))))
  
  name.file <- paste0(dir.store, "/sig-contig-tcga.fa")
  write.table(preFA, file = name.file, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Table assig ctg -> contig
  assignCtg <- as.data.frame(cbind(paste0("ctg_",seq(1:length(sigTCGA))), sigTCGA))
  names(assignCtg) <- c("Alias", "contig")
  # Call blast command
  cmdRun <- cmdRun <- paste0("bash -c \"blastn -db /store/EQUIPES/SSFA/Index/Gencode/gencode.v34.transcripts.fa -query ",
                   dir.store, "/sig-contig-tcga.fa -evalue 1e-3 -max_target_seqs 1 -outfmt 6 -max_hsps 1 -word_size 10 ",
                   " > ", dir.store,"/sig-contig-tcga-geneName.tsv\"")
  system(cmdRun)    
  
  #Open file geneName
  contigName <- as.data.frame(fread(input = paste0(dir.store,"/sig-contig-tcga-geneName.tsv"), sep = "\t", header = FALSE))
  contigName <- contigName[,c(1,2)]
  names(contigName) <- c("Alias", "Name")
  
  for (i in 1:dim(contigName)[1]){
    
    contigName$Name[[i]] = unlist(strsplit(contigName$Name[[i]],"\\|"))[6]
    print(contigName$Name[[i]])
  }
  
  # Mapping 
  mapAlliasName <- merge(assignCtg, contigName,
                         by.x = "Alias",
                         by.y = "Alias",
                         all.x = TRUE,
                         all.y = FALSE)
  uniqueGene <- c()
  dupGene <- c()
  for (i in 1:dim(mapAlliasName)[1]){
    if(is.na(mapAlliasName$Name[i])){
      
      uniqueGene[i] = mapAlliasName$Alias[i]
      dupGene[i] = mapAlliasName$Alias[i]
      
    }else{
      if(mapAlliasName$Name[i] %in% uniqueGene){
        dupGene[i] = paste0(mapAlliasName$Name[i], ".v",i)
      }else{
        dupGene[i] = mapAlliasName$Name[i]
      }
      uniqueGene[i] = mapAlliasName$Name[i]
    }
  }
  mapAlliasName$aliasGene = dupGene
  mapAlliasName$uniqueGene = uniqueGene
  contigMapp <- mapAlliasName[, c("contig", "aliasGene")]
  contigMappSave <- mapAlliasName[, c("contig", "uniqueGene")]
  
  #Save to file
  name.file <- paste0(dir.store, "/sig-contig-tcga.tsv")
  write.table(mapAlliasName[, c("contig")], file = name.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)   
  
  name.file <- paste0(dir.store, "/sig-contig-genemap-tcga.tsv")
  write.table(contigMappSave, file = name.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)   

  # dataframe of signature in TCGA dataset
  dataSigTCGA <- dataTCGA[,c("condition",sigTCGA)]
  
  #Order the geneName base on the contig sequence
  # Transpose contigMapp with colname is contig
  contigName <- as.data.frame(t(contigMapp[,-c(1)]))
  colnames(contigName) <- contigMapp$contig
  contigName = as.data.frame(contigName[, sigTCGA])

  ctgNameDis = unname(unlist(contigName[1,]))

  # Change the name of dataSigTCGA to mapping genenome
  names(dataSigTCGA) <- c("condition", ctgNameDis)  

  name.file <- paste0(dir.store, "/data-sig-contig-tcga.tsv")
  dataSaveTCGA <- dataTCGA[,c("Sample","condition",sigTCGA)]
  colnames(dataSaveTCGA) <- c("Sample", "condition", ctgNameDis)
  write.table(dataSaveTCGA, file = name.file, sep = "\t", quote = FALSE, row.names = FALSE)
  ########## finding signature in validation set ##############
    
  # Call the sript infer
  cmdRun <- paste0("bash -c \"/store/EQUIPES/SSFA/MEMBERS/thi-ngoc-ha.nguyen/Haoliang_Find_Kmer/kmer-filter/kmerFilter -n -k 31 -f ",
                   dir.store, "/sig-contig-tcga.tsv <(zcat ", dataValidPath ,") > ", dir.store,"/sig-contig-valid-ICGC.tsv\"")
  system(cmdRun)
  
  #Total kmers in this signature
  totalKmerContig <- fromContig2Kmers(sigDiscovery = sigTCGA)
  
  # loading list of kmer found in ICGC dataset based on method
  dataValidPath <- paste0(dir.store,"/sig-contig-valid-ICGC.tsv")
  
  countKmerValid <- normalizeContig(validCountPath = dataValidPath, libSizePath = totalKmers )
  
  # Classify kmer belong to contig in ICGC dataset
  kmerValid <- rownames(countKmerValid)
  listClf <- list()
  for(i in 1:length(totalKmerContig)){
    
    tmp <- c()
    for(kmer in kmerValid){
      
      if(kmer%in%totalKmerContig[[i]])
        tmp <- c(tmp, kmer)
      
    }
    listClf[[i]]<- tmp
  }
  names(listClf) <- sigTCGA[1:length(listClf)]
  
  #Extract contigs present in validation dataset
  listClf[sapply(listClf, is.null)] <- NULL
  
  # Calculate the count for each contig - based on mean count
  dataValidRepre <- data.frame(matrix(ncol = dim(countKmerValid)[2], nrow = 0))
  
  for(i in 1:length(listClf)){
    
    frame <- countKmerValid[c(listClf[[i]]),]
    repreCtg <- apply(data.matrix(frame),2, FUN = median)
    dataValidRepre <- rbind(dataValidRepre, repreCtg)
  }
  rownames(dataValidRepre) <- names(listClf)
  colnames(dataValidRepre) <- names(countKmerValid)
  
  strInfo <-paste("Contig signature include:", length(sigTCGA),"contigs that are created", length(unlist(totalKmerContig)),"kmer(s)\n",
                  "There are: ",length(kmerValid),"kmers found in validation and belong to",length(listClf),"contigs")
  
  print(strInfo)
  name.file <- paste0(dir.store, "/info-icgc.txt")
  write.csv(strInfo, name.file)
  
  probesID <- rownames(dataValidRepre)
  dataSignatureValid <- dataValidRepre

  ########## processing of conditions ##############
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
  dataSigICGC <- dataICGC[,c("condition", probesID)]
  #print(probesID)
  ctgValidName <- c()
  # Mapping contig to genome
  for ( i in 1: length(probesID)){
    
    ctgValidName[i] = unlist(unname(contigName[probesID[i]]))
    print(ctgValidName[i])
  }
  
  names(dataSigICGC) <- c("condition", ctgValidName)
 
  name.file <- paste0(dir.store, "/data-sig-contig-icgc.tsv")
  dataSaveICGC <- dataICGC[,c("Sample","condition", probesID)]
  colnames(dataSaveICGC) <- c("Sample", "condition", ctgValidName)
  write.table(dataSaveICGC, file = name.file, sep = "\t", quote = FALSE, row.names = FALSE)
  ########## signature performance ##############
  #Combine data in TCGA and ICGC
  absentCtg <- setdiff(names(dataSigTCGA),names(dataSigICGC))
  print(paste("List", length(names(dataSigTCGA))-1,"contigs in signature that found in TCGA."))
  print(names(dataSigTCGA)[-1])
  
  print(paste("There are", length(unique(uniqueGene)), "unique genes in signature in TCGA."))
  
  print("List contig signature are found in ICGC")
  print(names(dataSigICGC)[-1])
  print("Contig Absentce:")
  print(absentCtg)
  if (length(absentCtg) != 0){
    
    for (contig in absentCtg){
      dataSigICGC$contig <- 0
      names(dataSigICGC)[length(names(dataSigICGC))] <- contig
    }
  }
  dataSigICGC <- dataSigICGC[,names(dataSigTCGA)]   
  
  resSign <- takeDataReturnAUC(frameTrain = dataSigTCGA, frameTest = dataSigICGC, absentContig = absentCtg, status = "risk")
  
  finalTab <- rbind(resSign[[3]], resSign[[1]], resSign[[2]])
  rownames(finalTab) <- c("AUC signature contig in ICGC", 
			"Mean AUC cross-valid (TCGA)", "SD AUC cross-valid (TCGA)")
  print(finalTab)
  
  name.file <- paste0(dir.store, "/auc-sig-cotig.tsv")
  write.table(finalTab, file = name.file, quote = FALSE, row.name = TRUE, col.names = FALSE)  
  
  # Store important variable 
  name.file <- paste0(dir.store, "/import-sig-contig.tsv")
  importVar <- as.data.frame(resSign[[4]]$importance)
  importVar$Gene <- rownames(resSign[[4]]$importance)
  names(importVar) <- c("Overall", "Contig_map_to_gene")
  dataStore <- importVar[round(order(importVar$Overall, decreasing = T),4),]

  print("Contribution of each contig in signature:")
  print(dataStore[, c("Contig_map_to_gene","Overall")])
  write.table(dataStore[, c("Contig_map_to_gene", "Overall")], file = name.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
}
#################################
#Run the script to infer

resContig <- pipeline(topProbesPath = topContig, samplesConditionDisPath = sampleTCGA, dataValidPath = kmerICGC,
                    samplesConditionValidPath = sampleICGC, numruns = NUM_RUNS)

#################################
