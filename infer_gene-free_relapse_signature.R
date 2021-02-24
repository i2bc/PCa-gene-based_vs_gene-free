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

# Discovery data using top 500 contigs
dir.discovery <- "Data_discovery/Relapse/"

topContig <- paste0(dir.discovery, "top500_contig_merge_norm.nb5.tsv")
sampleTCGA <- paste0(dir.discovery, "sample_conditions.tsv")

# Validation ICGC data
dir.validation <- "Data_validation/Relapse/"

sampleICGC <- paste0(dir.validation, "sample_cond_icgc.tsv")
kmerICGC <- paste0(dir.validation, "raw-counts.tsv.gz")
totalKmersICGC <- paste0(dir.validation, "sum_counts_icgc.tsv")

# Validation Stello data
sampleStelloo <- paste0(dir.validation, "sample_cond_stelloo.tsv")
kmerStelloo <- paste0(dir.validation, "raw-counts-Stello-Relapse-3-3.tsv.gz")
totalKmersStelloo <- paste0(dir.validation, "sum_counts_stelloo.tsv")

# Number of time for sampling
NUM_RUNS=100

# Directory of kmerFilter tool 
dir.kmerFilter <- "/store/EQUIPES/SSFA/MEMBERS/thi-ngoc-ha.nguyen/Haoliang_Find_Kmer/kmer-filter/kmerFilter"

# Directory of BLAST database
dir.blastDB <- "/store/EQUIPES/SSFA/Index/Gencode/gencode.v34.transcripts.fa"

# Directory to store result
dir.store <- paste0("Result_infer_signature/Relapse/gene_free")
dir.create(file.path(dir.store), showWarnings = FALSE, recursive = TRUE)
############################################################################################
evaluate_contig <- function(ctgMappDis, dataSigDis, sigContigDis, sigContigDisPath, dataValidPath, samplesConditionValidPath, modeValid){
 
 # Call the sript infer
 cmdRun <- paste0("bash -c \"", dir.kmerFilter, " -n -k 31 -f ",
                  sigContigDisPath, " <(zcat ", dataValidPath ,") > ", dir.store,"/sig-contig-valid-", modeValid, ".tsv\"")
 system(cmdRun)
 
 #Total kmers in this signature
 totalKmerContig <- fromContig2Kmers(sigDiscovery = sigContigDis)
 
 # loading list of kmer found in ICGC dataset based on method
 dataValidPath <- paste0(dir.store,"/sig-contig-valid-", modeValid, ".tsv")
 
 if(modeValid == "icgc"){
   countKmerValid <- normalizeContig(validCountPath = dataValidPath, libSizePath = totalKmersICGC )
 }else{
   countKmerValid <- normalizeContig(validCountPath = dataValidPath, libSizePath = totalKmersStelloo )
 }
 
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
 names(listClf) <- sigContigDis[1:length(listClf)]
 
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
 
 strInfo <-paste("Contig signature include:", length(sigContigDis),"contigs that are created", length(unlist(totalKmerContig)),"kmer(s)\n",
                 "There are: ",length(kmerValid),"kmers found in validation and belong to",length(listClf),"contigs")
 
 print(strInfo)
 name.file <- paste0(dir.store, "/info-", modeValid, ".txt")
 write.csv(strInfo, name.file)
 
 probesID <- rownames(dataValidRepre)
 dataSignatureValid <- dataValidRepre
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
 
 # dataframe of signature in ICGC dataset
 dataSigValid <- dataValid[,c("condition", probesID)]
 dataSave <- dataValid[,c("Sample", "condition", probesID)]
 
 print(paste("List of contig from TCGA found in", modeValid, "validation:"))
 print(probesID)
 ctgValidName <- c()
 # Mapping contig to genome
 for ( i in 1: length(probesID)){
   
   ctgValidName[i] = ctgMappDis[ctgMappDis$contig == probesID[i],]$UniqueName
   print(ctgValidName[i])
 }
 
 names(dataSigValid) <- c("condition", ctgValidName)
 names(dataSave) <- c("Sample","condition", ctgValidName)
 ########## signature performance ##############
 #Combine data in TCGA and validation set
 absentCtg <- setdiff(names(dataSigDis),names(dataSigValid))
 print(paste(length(absentCtg), "contigs are not found in", modeValid, "valid set:"))
 print(absentCtg)
 if (length(absentCtg) != 0){
   
   for (contig in absentCtg){
     dataSigValid$contig <- 0
     names(dataSigValid)[length(names(dataSigValid))] <- contig
   }
 }
 dataSigValid <- dataSigValid[,names(dataSigDis)]   
 
 
 resSignROC <- takeDataReturnAUC(frameTrain = dataSigDis, frameTest = dataSigValid, absentContig = absentCtg, status = "relapse")
 
 resSignPR <- takeDataReturnPR(frameTrain = dataSigDis, frameTest = dataSigValid, absentContig = absentCtg, status = "relapse")
 
 return (list(resSignROC, dataSave, resSignPR))
}

################################
pipeline <- function(topProbesPath, samplesConditionDisPath, numruns){
 
 # loading top  genes/contigs in TCGA discovery dataset based on approach level
 countTopProbe <- as.data.frame(fread(topProbesPath, sep="\t", header = TRUE))
 
 countTopProbe <- data.frame(row.names = countTopProbe$feature, countTopProbe[,-c(1:4)], check.names=F)    
 
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
 
 print("List of contig in the signature:")
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
 cmdRun <- paste0("bash -c \"blastn -db ", dir.blastDB, " -query ",
                  dir.store, "/sig-contig-tcga.fa -evalue 1e-3 -max_target_seqs 1 -outfmt 6 -max_hsps 1 -word_size 10 ",
                  #   '|awk -F"\t"', "'{print $1" , '"\t"', "$2}'", '|awk -F"|"', " '{print $1",' "\t" ', "$6}'",'|awk -F"\t"'," '{print $1",' "\t"', " $3}'",
                  " > ", dir.store,"/sig-contig-tcga-geneName.tsv\"")
 system(cmdRun)      
 
 #Open file geneName
 contigName <- as.data.frame(fread(input = paste0(dir.store,"/sig-contig-tcga-geneName.tsv"), sep = "\t", header = FALSE))
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
 uniqueGene <- c()
 for (i in 1:dim(mapAlliasName)[1]){
   if(is.na(mapAlliasName$Name[i])){
     
     uniqueGene[i] = mapAlliasName$Alias[i]
     
   }else{
     if(mapAlliasName$Name[i] %in% uniqueGene){
       uniqueGene[i] = paste0(mapAlliasName$Name[i], ".v",i)
     }else{
       uniqueGene[i] = mapAlliasName$Name[i]
     }
   }
 }
 
 mapAlliasName$UniqueName = uniqueGene
 contigMapp <- mapAlliasName[, c("contig", "UniqueName")]
 
 #Save to file
 name.file <- paste0(dir.store, "/sig-contig-tcga.tsv")
 write.table(mapAlliasName[, c("contig")], file = name.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)   
 
 name.file <- paste0(dir.store, "/sig-contig-genemap-tcga.tsv")
 write.table(contigMapp, file = name.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)       
 
 # dataframe of signature in TCGA dataset
 dataSigTCGA <- dataTCGA[,c("condition", sigTCGA)]
 
 #Order the geneName base on the contig sequence
 # Transpose contigMapp with colname is contig
 contigName <- as.data.frame(t(contigMapp[,-c(1)]))
 colnames(contigName) <- contigMapp$contig
 contigName = as.data.frame(contigName[, sigTCGA])
 
 ctgNameDis = unname(unlist(contigName[1,]))
 
 print("List of contig mapping to genome: ")
 print(ctgNameDis)
 
 name.file <- paste0(dir.store, "/data-sig-contig-tcga.tsv")
 dataSaveTCGA <- dataTCGA[,c("Sample","condition",sigTCGA)]
 colnames(dataSaveTCGA) <- c("Sample", "condition", ctgNameDis)
 write.table(dataSaveTCGA, file = name.file, sep = "\t", quote = FALSE, row.names = FALSE)
 names(dataSigTCGA) <- c("condition", ctgNameDis)
 ########## finding signature in validation set ##############
 
 #Valid in ICGC
 resSigICGC <- evaluate_contig(ctgMappDis = contigMapp, dataSigDis = dataSigTCGA, sigContigDis = sigTCGA, sigContigDisPath = paste0(dir.store, "/sig-contig-tcga.tsv"), 
                               dataValidPath = kmerICGC, samplesConditionValidPath = sampleICGC, modeValid = "icgc")
 

 #Store data frame valid
 name.file <- paste0(dir.store, "/data-sig-contig-icgc.tsv")
 write.table(resSigICGC[[2]], file = name.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

 finalTab <- rbind(resSigICGC[[1]][[1]], resSigICGC[[1]][[2]], resSigICGC[[1]][[3]], resSigICGC[[1]][[4]],
                   resSigICGC[[1]][[5]], resSigICGC[[1]][[6]], resSigICGC[[1]][[7]], resSigICGC[[1]][[8]],
                   resSigICGC[[1]][[9]], resSigICGC[[1]][[10]], resSigICGC[[1]][[11]], resSigICGC[[1]][[12]])
 rownames(finalTab) <- c( 
   "Mean ROC_AUC cv (TCGA)", "SD ROC_AUC cv (TCGA)", 
   "Mean ROC_AUC down cv (TCGA)", "SD ROC_AUC down cv (TCGA)",
   "Mean ROCAUC up cv (TCGA)", "SD ROC_AUC up cv (TCGA)",
   "ROC_AUC sig contig in ICGC", "PrAUC sig contig in ICGC",
   "ROC_AUC down sig contig in ICGC", "PrAUC down sig contig in ICGC",
   "ROC_AUC up sig contig in ICGC", "PrAUC up sig contig in ICGC")
 print(finalTab)
 
 name.file <- paste0(dir.store, "/Roc_auc-sig-contig-icgc.tsv")
 write.table(finalTab, file = name.file, quote = FALSE, row.name = TRUE, col.names = FALSE) 
 
 finalTab <- rbind(resSigICGC[[3]][[1]], resSigICGC[[3]][[2]], resSigICGC[[3]][[3]], resSigICGC[[3]][[4]],
                   resSigICGC[[3]][[5]], resSigICGC[[3]][[6]], resSigICGC[[3]][[7]], resSigICGC[[3]][[8]],
                   resSigICGC[[3]][[9]], resSigICGC[[3]][[10]], resSigICGC[[3]][[11]], resSigICGC[[3]][[12]])
 rownames(finalTab) <- c( 
   "Mean PR_AUC cv (TCGA)", "SD PR_AUC cv (TCGA)", 
   "Mean PR_AUC down cv (TCGA)", "SD PR_AUC down cv (TCGA)",
   "Mean PR_AUC up cv (TCGA)", "SD PR_AUC up cv (TCGA)",
   "ROC_AUC sig contig in ICGC", "PrAUC sig contig in ICGC",
   "ROC_AUC down sig contig in ICGC", "PrAUC down sig contig in ICGC",
   "ROC_AUC up sig contig in ICGC", "PrAUC up sig contig in ICGC")
 print(finalTab)  
 
 name.file <- paste0(dir.store, "/Pr_auc-sig-contig-icgc.tsv")
 write.table(finalTab, file = name.file, quote = FALSE, row.name = TRUE, col.names = FALSE)   
 
 
  
 # Valid in Stelloo
 resSigStelloo <- evaluate_contig(ctgMappDis = contigMapp, dataSigDis = dataSigTCGA,sigContigDis = sigTCGA, sigContigDisPath = paste0(dir.store, "/sig-contig-tcga.tsv"), 
                                  dataValidPath = kmerStelloo, samplesConditionValidPath = sampleStelloo, modeValid = "stello")
 
 

 #Store data frame valid
 name.file <- paste0(dir.store, "/data-sig-contig-stelloo.tsv")
 write.table(resSigStelloo[[2]], file = name.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
 
 finalTab <- rbind(resSigStelloo[[1]][[1]], resSigStelloo[[1]][[2]], resSigStelloo[[1]][[3]], resSigStelloo[[1]][[4]],
                   resSigStelloo[[1]][[5]], resSigStelloo[[1]][[6]], resSigStelloo[[1]][[7]], resSigStelloo[[1]][[8]],
                   resSigStelloo[[1]][[9]], resSigStelloo[[1]][[10]], resSigStelloo[[1]][[11]], resSigStelloo[[1]][[12]])
 rownames(finalTab) <- c( 
   "Mean ROC_AUC cv (TCGA)", "SD ROC_AUC cv (TCGA)", 
   "Mean ROC_AUC down cv (TCGA)", "SD ROC_AUC down cv (TCGA)",
   "Mean ROCAUC up cv (TCGA)", "SD ROC_AUC up cv (TCGA)",
   "ROC_AUC sig contig in Stelloo", "PrAUC sig contig in Stelloo",
   "ROC_AUC down sig contig in Stelloo", "PrAUC down sig contig in Stelloo",
   "ROC_AUC up sig contig in Stelloo", "PrAUC up sig contig in Stelloo")
 print(finalTab)
 
 name.file <- paste0(dir.store, "/Roc_auc-sig-contig-stelloo.tsv")
 write.table(finalTab, file = name.file, quote = FALSE, row.name = TRUE, col.names = FALSE) 
 
 finalTab <- rbind(resSigStelloo[[3]][[1]], resSigStelloo[[3]][[2]], resSigStelloo[[3]][[3]], resSigStelloo[[3]][[4]],
                   resSigStelloo[[3]][[5]], resSigStelloo[[3]][[6]], resSigStelloo[[3]][[7]], resSigStelloo[[3]][[8]],
                   resSigStelloo[[3]][[9]], resSigStelloo[[3]][[10]], resSigStelloo[[3]][[11]], resSigStelloo[[3]][[12]])
 rownames(finalTab) <- c( 
   "Mean PR_AUC cv (TCGA)", "SD PR_AUC cv (TCGA)", 
   "Mean PR_AUC down cv (TCGA)", "SD PR_AUC down cv (TCGA)",
   "Mean PR_AUC up cv (TCGA)", "SD PR_AUC up cv (TCGA)",
   "ROC_AUC sig contig in Stelloo", "PrAUC sig contig in Stelloo",
   "ROC_AUC down sig contig in Stelloo", "PrAUC down sig contig in Stelloo",
   "ROC_AUC up sig contig in Stelloo", "PrAUC up sig contig in Stelloo")
 print(finalTab)  
 
 name.file <- paste0(dir.store, "/Pr_auc-sig-contig-stelloo.tsv")
 write.table(finalTab, file = name.file, quote = FALSE, row.name = TRUE, col.names = FALSE)  
}
#################################
#Run the script to infer
resContig <- pipeline(topProbesPath = topContig, samplesConditionDisPath = sampleTCGA, numruns = NUM_RUNS)
