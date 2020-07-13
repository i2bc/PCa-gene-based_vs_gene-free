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

library(knitr)
library(data.table)
library(RColorBrewer)
################
# Load list of functions 
################
source("useful_functions.R")
#######################################################################
## scripts
#######################################################################

# Discovery data using top 500 NB5 gene and contig

topGene <-"Data_discovery/top2k_TCGA_gene_norm.nb5.out"
topContig <- "Data_discovery/top2k_contig_norm_merge.nb5.out"

sampleTCGA <-"Data_discovery/sample_conditions.tsv"

# Validation data
sampleICGC <- "Data_validation/sample_conditions.tsv"
geneICGC <- "Data_validation/gene-counts-ICGC148-LRHR.norm.tsv"
kmerICGC <- "Data_validation/raw-counts.tsv.gz"
totalKmers <-"Data_validation/sum_counts.tsv"


NUM_RUNS=100
nKeep = 500

# Directory to store result
dir.store <- paste0("Result_infer_signature/",nKeep)
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
pipeline <- function(topProbesPath, samplesConditionDisPath, dataValidPath, samplesConditionValidPath, numruns, appLevel){
  
  # loading top 500 genes/contigs in TCGA discovery dataset based on approach level
  countTopProbe <- as.data.frame(fread(topProbesPath, sep="\t", header = TRUE))
  if(appLevel == "gene"){
    
    countTopProbe <- data.frame(row.names = countTopProbe$feature, countTopProbe[,-c(1,2)], check.names=F)
    
  }
  if(appLevel == "contig"){
    
    countTopProbe <- data.frame(row.names = countTopProbe$feature, countTopProbe[,-c(1:4)], check.names=F)    
  
  }
  countTopProbe <- countTopProbe[1:nKeep,]
  top1k <- c(rownames(countTopProbe))
  ########## processing of conditions ##############
  samplesConditionDis<-as.data.frame(fread(samplesConditionDisPath,sep="\t", header= FALSE ,check.names=F))
  names(samplesConditionDis)<-c("Sample","condition")
  
  #Transpose
  countTopProbe <- as.data.frame(t(countTopProbe))
  
  if(appLevel == "contig"){
    #Change the double value in right format
    for (i in 1:dim(countTopProbe)[2]){
      countTopProbe[,i] <- as.double(as.character(countTopProbe[,i]))
    }
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
  
  
  if(appLevel == "contig"){
    
    print(sigTCGA)
    
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
    cmdRun <- paste0("\"blastn -db /store/EQUIPES/SSFA/Index/Gencode/gencode.v34.transcripts.fa -query ",
                     dir.store, "/sig-contig-tcga.fa -evalue 1e-3 -max_target_seqs 1 -outfmt 6 -max_hsps 1 -word_size 10 ",
                     " > ", dir.store,"/sig-",appLevel,"-tcgc-geneName.tsv\"")
    system(cmdRun)    
    
    #Open file geneName
    contigName <- as.data.frame(fread(input = paste0(dir.store,"/sig-",appLevel,"-tcgc-geneName.tsv"), sep = "\t", header = FALSE))
    contigName <- contigName[,c(1,2)]
    names(contigName) <- c("Alias", "Name")

    for (i in 1:dim(contigName)[1]){
      
      contigName$Name[[i]] = unlist(strsplit(contigName$Name[[i]],"\\|"))[6]
    }
    
    # Mapping 
    mapAlliasName <- merge(assignCtg, contigName,
                           by.x = "Alias",
                           by.y = "Alias",
                           all.x = TRUE,
                           all.y = FALSE)
    mapAlliasName$Alias <- levels(droplevels(mapAlliasName$Alias))
    for (i in 1:dim(mapAlliasName)[1]){
      if(is.na(mapAlliasName$Name[i])){
        
        mapAlliasName$Name[i] <- mapAlliasName$Alias[i]

      }

    }
    #Save to file
    name.file <- paste0(dir.store, "/contig-alias-geneName.tsv")
    write.table(mapAlliasName, file = name.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)   
    
    # Save contig signature
    name.file <- paste0(dir.store, "/sig-contig-tcga.txt")
    write.table(sigTCGA, file = name.file, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  if(appLevel == "gene"){
    
    # From gene symbol to gene name
    signature.idsym <- fromGeneIDtakeGenName(sigTCGA)
    signature.save <- signature.idsym$SYMBOL
    
    print(signature.save)
    
    name.file <- paste0(dir.store, "/sig-gene.tsv")
    
    write.table(signature.idsym, file = name.file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    
  }
  
  
  # dataframe of signature in TCGA dataset
  dataSigTCGA <- dataTCGA[,c("condition",sigTCGA)]
  
  ########## finding signature in validation set ##############
  if(appLevel == "contig"){
    
    # Call the sript infer
    cmdRun <- paste0("bash -c \"/store/EQUIPES/SSFA/MEMBERS/thi-ngoc-ha.nguyen/Haoliang_Find_Kmer/kmer-filter/kmerFilter -n -k 31 -f ",
                     dir.store, "/sig-contig-tcga.tsv <(zcat ", dataValidPath ,") > ", dir.store,"/sig-",appLevel,"-valid-ICGC.tsv\"")
    system(cmdRun)
    
    #Total kmers in this signature
    totalKmerContig <- fromContig2Kmers(sigDiscovery = sigTCGA)
    
    # loading list of kmer found in ICGC dataset based on method
    dataValidPath <- paste0(dir.store,"/sig-",appLevel,"-valid-ICGC.tsv")
    
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
    nameProbe <- "contig"
    
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
    
    ########## signature performance ##############
    #Combine data in TCGA and ICGC
    absentCtg <- setdiff(names(dataSigTCGA),names(dataSigICGC))
    if (length(absentCtg) != 0){
      
      for (contig in absentCtg){
        dataSigICGC$contig <- 0
        names(dataSigICGC)[length(names(dataSigICGC))] <- contig
      }
    }
    dataSigICGC <- dataSigICGC[,names(dataSigTCGA)]   
    
    resSign <- takeDataReturnAUC_ICGC(frameTrain = dataSigTCGA, frameTest = dataSigICGC, absentContig = absentCtg)
    
    name.file <- paste0(dir.store, "/coef-sig-contig.tsv")
    
    dataStore <- resSign[[3]][-c(1)]
    
    #Mapp contig to gene nome
    names(dataStore) <- fromGeneIDtakeGenName(names(dataStore))$SYMBOL
    dataStore <- round(dataStore[order(dataStore, decreasing = TRUE)],4)
  }
  if(appLevel == "gene"){
    
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
    
    resSign <- takeDataReturnAUC_ICGC(frameTrain = dataSigTCGA, frameTest = dataSigICGC, absentContig = c())
    
    name.file <- paste0(dir.store, "/coef-sig-gene.tsv")
    dataStore <- resSign[[3]][-c(1)]
    names(dataStore) <- fromGeneIDtakeGenName(names(dataStore))$SYMBOL
    dataStore <- round(dataStore[order(dataStore, decreasing = TRUE)],4)
  }
  
  print(paste("Performance of", appLevel,"in TCGA:", resSign[[1]]))
  print(paste("Performance of", appLevel,"in ICGC:", resSign[[2]]))
  
  # Store coefficent 
  write.table(dataStore, file = name.file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
  
  resSignCross <- takeDataReturnAUC_cross(frame = dataSigTCGA, NUM_RUNS = NUM_RUNS)
  print(paste("Performance of", appLevel,"in TCGA (cross validation):", resSignCross[[1]], "+/-", resSignCross[[2]]))
  
  if(appLevel == "gene"){
    name.file <- paste0(dir.store, "/data-sig-gene-tcga.tsv")
    dataSaveTCGA <- dataTCGA[,c("Sample","condition",sigTCGA)]
    colnames(dataSaveTCGA) <- c("Sample", "condition", signature.save)
    write.table(dataSaveTCGA, file = name.file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    name.file <- paste0(dir.store, "/data-sig-gene-icgc.tsv")
    dataSaveICGC <- dataICGC[,c("Sample","condition",sigTCGA)]
    colnames(dataSaveICGC) <- c("Sample", "condition", signature.save)
    write.table(dataSaveICGC, file = name.file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    return (list(signature.save, resSignCross, dataSaveTCGA, dataSaveICGC, resSign[[2]], resSign[[3]]))
  }
  
  if(appLevel == "contig"){
    name.file <- paste0(dir.store, "/data-sig-contig-tcga.tsv")
    dataSaveTCGA <- dataTCGA[,c("Sample","condition",sigTCGA)]
    write.table(dataSaveTCGA, file = name.file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    name.file <- paste0(dir.store, "/data-sig-contig-icgc.tsv")
    dataSaveICGC <- dataICGC[,c("Sample","condition",probesID)]
    write.table(dataSaveICGC, file = name.file, sep = "\t", quote = FALSE, row.names = FALSE)
    return (list(sigTCGA, resSignCross, dataSaveTCGA, dataSaveICGC, probesID, resSign[[2]]))
  }

}
#################################
#Run the script to infer
resGene <- pipeline(topProbesPath = topGene, samplesConditionDisPath = sampleTCGA, dataValidPath = geneICGC,
                       samplesConditionValidPath = sampleICGC, numruns = NUM_RUNS, appLevel = "gene")
#resContig <- pipeline(topProbesPath = topContig, samplesConditionDisPath = sampleTCGA, dataValidPath = kmerICGC,
#                    samplesConditionValidPath = sampleICGC, numruns = NUM_RUNS, appLevel = "contig")

#################################
#Box plot
boxPlotGeneTCGA <- boxPlot(dataPlot = resGene[[3]], signature = resGene[[1]], modeLevel = "gene", dir.name = dir.store, status = "tcga")
png(filename = paste0(dir.store, "/box-plot-sig-gene-TCGA.png"), width = 3200, height = 1600, units = "px", res = 300)
plot(boxPlotGeneTCGA)
dev.off()

boxPlotGeneICGC <- boxPlot(dataPlot = resGene[[4]], signature = resGene[[1]], modeLevel = "gene", dir.name = dir.store, status = "icgc")
png(filename = paste0(dir.store, "/box-plot-sig-gene-ICGC.png"), width = 3200, height = 1600, units = "px", res = 300)
plot(boxPlotGeneICGC)
dev.off()
################################
# Process for coefficent
coefGene <- as.data.frame(fread(input = paste0(dir.store, "/coef-sig-gene.tsv"), sep = "\t", header = FALSE))
names(coefGene) <- c("ID", "coef")
coefGene$coef <- round(coefGene$coef,4)
coefGene$ID <- gsub("\\..*","",coefGene$ID)
# GeneNam- gene ID
geneName <- as.data.frame(fread(input = paste0(dir.store, "/sig-gene-tcga.tsv"), sep = "\t", header = FALSE))
names(geneName) <- c("Name", "ID")                          
#Mapping
geneMerge <- merge(geneName, coefGene,
                  by.x = "ID",
                  by.y = "ID",
                  all.x = TRUE,
                  all.y = TRUE)
geneMerge <- geneMerge[-c(1)]
geneMerge <- geneMerge[order(geneMerge$coef, decreasing = TRUE),]
name.file <- paste0(dir.store, "/coef-gene-order.tsv")
write.table(geneMerge, file = name.file, sep = "\t", quote = FALSE, row.names = FALSE)