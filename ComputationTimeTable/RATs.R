
library(gtools)
library(rats)
library(data.table)
library(tidyr)
#infReps <- "GibbsThin100"
infReps <- "Boot"
#dataset <- "SQCC"
dataset <- "GEUV1"
DRIMSeqFiltering <- TRUE

TwentySamplesTotalAnalysis <- TRUE



source("/Users/Scott/Documents/Dissertation/code/CompFunctions.R")

if(dataset=="SQCC"){
   def_wd <- "/Users/Scott/Documents/Dissertation Data/SQCCDataReproduceOldResBeforeCommonCode/"
   save_dir <- "/Users/Scott/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/SQCC/"
   npartstouse <- 10
   
   wd2 <- "/Users/Scott/Documents/Dissertation/res/SQCCReproduceOldResBeforeCommonCode"
   
}else if(dataset=="GEUV1"){
  def_wd <- "/Users/Scott/Documents/Dissertation Data/GEUV1Data/"
  if(TwentySamplesTotalAnalysis==TRUE){
    save_dir <- "/Users/Scott/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/GEUV1TwentySamples/"
  }else{
    save_dir <- "/Users/Scott/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/GEUV1/"
  }
  npartstouse <- 100
  
  wd2 <- "/Users/Scott/Documents/Dissertation/res/GEUV1Data"
}

setwd(def_wd)

func_loc <- "/Users/Scott/Documents/Dissertation/code/CompFunctions.R"
source(func_loc)
if(infReps == "Boot" & DRIMSeqFiltering==TRUE){
  RATs_InfRep_save_dir <- paste0(def_wd, "RATsInfRepsBootDRIMSeqFiltering/")
}else if(infReps == "Boot" & DRIMSeqFiltering==FALSE){
  RATs_InfRep_save_dir <- paste0(def_wd, "RATsInfRepsBoot/")
}else if(infReps=="GibbsThin100"){
  RATs_InfRep_save_dir <- paste0(def_wd, "RATsInfReps/")
}

 
 # #Only load the ilrMeansCovs to get the list of genes used to match the Compositional Analysis
 # if(DRIMSeqFiltering==TRUE){
 #   ilrMeansCovs <- loadRData(paste0(def_wd,"ilrMeansCovsBoot/ilrMeansCovsNoOtherGroupsFilteredPart1.RData"))
 #   if(npartstouse > 1){
 #     for(i in 2:npartstouse){
 #       ilrMeansCovs <- c(ilrMeansCovs, loadRData(paste0(def_wd,"ilrMeansCovsBoot/ilrMeansCovsNoOtherGroupsFilteredPart", i, ".RData")))
 #     }
 #   }
 # }else if(DRIMSeqFiltering==FALSE){
 #   ilrMeansCovs <- loadRData(paste0(def_wd,"ilrMeansCovsBoot/ilrMeansCovsOtherGroupsPart1.RData"))
 #   if(npartstouse > 1){
 #     for(i in 2:npartstouse){
 #       ilrMeansCovs <- c(ilrMeansCovs, loadRData(paste0(def_wd,"ilrMeansCovsBoot/ilrMeansCovsOtherGroupsPart", i, ".RData")))
 #     }
 #   }
 # }
 #genestouse <- names(ilrMeansCovs)[1:500]
 
 if(DRIMSeqFiltering==TRUE){
   FilteringMethod <- "DRIMSeq"

   load("cntDatasetscntsScaledTPMNoOtherGroupsFiltered.RData")
   
   cntDatasets <- cntDatasetsFiltered
   
 }else if(DRIMSeqFiltering==FALSE){
   FilteringMethod <- "OtherGroups"
   load("cntDatasetscntsScaledTPM.RData")
 }

 #genestouse <- names(cntDatasets)[1:100]
 genestouse <- names(cntDatasets)
 cntDatasetsToUse <- cntDatasets[genestouse]
 #load("cntGenecntsScaledTPM.RData")
 cntGene <- loadRData("cntGenecntsScaledTPMFiltered.RData")

# if(infReps=="GibbsThin100"){
#   sub_dir1 <- "SalmonReproduceResBeforeCommonCodeThin100/"
# }else if(infReps=="Gibbs"){
#   sub_dir1 <- "SalmonReproduceResBeforeCommonCode/"
# }else if(infReps=="Boot"){
#   sub_dir1 <- "SalmonReproduceResBeforeCommonCodeBootSamps/"
# }


# for(i in 1:10){
#   curr_samp <- paste0("Sample", i)
#   StartT <- proc.time()
#   curr <- reshapeinfRepDataforRATs(curr_samp, cntDatasets = cntDatasets[1:100], GibbsSamps_dir = GibbsSamps_dir)
#   assign(paste0("RATsInfRepData", curr_samp, "CompTime"), proc.time() - StartT)
#   
#   assign(paste0("RATsInfRepData", curr_samp), curr)
#   
#   save(list = c(paste0("RATsInfRepData", curr_samp), paste0("RATsInfRepData", curr_samp, "CompTime")), file = paste0("RATsInfRepData", curr_samp, ".RData"))
#   
# }


# if(!file.exists("InfRepsForRATs.RData")){
#   sT <- proc.time()
#   InfRepsListCondA <- lapply(paste0("Sample", 1:5), reshapeinfRepDataforRATs, cntDatasets = cntDatasets, GibbsSamps_dir = GibbsSamps_dir)
#   names(InfRepsListCondA) <- paste0("Sample", 1:5)
#   InfRepsListCondB <- lapply(paste0("Sample", 6:10), reshapeinfRepDataforRATs, cntDatasets = cntDatasets, GibbsSamps_dir = GibbsSamps_dir)
#   names(InfRepsListCondB) <- paste0("Sample", 6:10)
#   infRepDatTime <- proc.time() - sT
#   save(InfRepsListCondA, InfRepsListCondB, infRepDatTime, file = "InfRepsForRATs.RData")
# }else{
#   load("InfRepsForRATs.RData")
# }
 
 if(dataset=="SQCC"){
   key_to_use <- key
   sampsCondA <- subset(key_to_use$Identifier, key$Condition=="A")
   sampsCondB <- subset(key_to_use$Identifier, key$Condition=="B")
   
   cntGeneSub <- cntGene
 }else if(dataset=="GEUV1"){
   if(TwentySamplesTotalAnalysis==TRUE){
     load("/Users/Scott/Documents/Dissertation Data/GEUV1Data/FilesForPowerAnalysis1.RData")
     
     s_c1 <- as.numeric(rownames(sub_key)[sub_key$Condition=="CEU"][1:10])
     s_c2 <- as.numeric(rownames(sub_key)[sub_key$Condition=="GBR"][1:10])
     
     #Order the sample numbers from smallest to largest
     s_ordered <- sort(c(s_c1, s_c2))
     
     #sampstouse <- c(paste0("Sample", s_c1), paste0("Sample", s_c2))
     sampstouse <- paste0("Sample", s_ordered)
     
     
     sampsCondA <- paste0("Sample", sort(s_c1))
     sampsCondB <- paste0("Sample", sort(s_c2))
     
     key_to_use <- subset(sub_key, sub_key$Identifier %in% sampstouse)
     #key_to_use <- subset(key, key$Identifier %in% sampstouse)
     
     #Remove the unneeded factor levels
     #key_to_use$Condition <- factor(key_to_use$Condition)
     #nsamp <- length(sampstouse)
     
     #Group <- key_to_use$Condition
   }else{
     key_to_use <- key
     sampstouse <- key$Identifier
   }
   
   cntGeneSub <- cntGene[,c("gene_id", "tx_id", paste0(sampstouse, "Cnt"))]
   
 }
InfRepsListCondA <- vector(mode = "list", length = length(sampsCondA))

for(k in 1:length(sampsCondA)){
  curr_samp <- sampsCondA[k]
  InfRepsListCondA[[k]] <- loadRData(paste0(RATs_InfRep_save_dir, "RATsInfRepData", curr_samp, ".RData"))
  if(DRIMSeqFiltering==TRUE){
    InfRepsListCondA[[k]] <- subset(InfRepsListCondA[[k]], InfRepsListCondA[[k]]$feature_id %in% cntGeneSub$tx_id)
  }
  
}

names(InfRepsListCondA) <- sampsCondA



InfRepsListCondB <- vector(mode = "list", length = length(sampsCondB))

for(k in 1:length(sampsCondB)){
  curr_samp <- sampsCondB[k]
  InfRepsListCondB[[k]] <- loadRData(paste0(RATs_InfRep_save_dir, "RATsInfRepData", curr_samp, ".RData"))
  if(DRIMSeqFiltering==TRUE){
    InfRepsListCondB[[k]] <- subset(InfRepsListCondB[[k]], InfRepsListCondB[[k]]$feature_id %in% cntGeneSub$tx_id)
  }
  
}

names(InfRepsListCondB) <- sampsCondB






RATsObsRes <- RATsObsAnalysis(cntGene = cntGeneSub, cntDatasets = cntDatasetsToUse, key = key_to_use, infRepDataA = InfRepsListCondA, 
                              infRepDataB = InfRepsListCondB, fullReturn = TRUE, DRIMSeqFiltering = DRIMSeqFiltering)


#Get an idea of the memory usage
print(gc())

if(!dir.exists(save_dir)){dir.create(save_dir, recursive = TRUE)}

if(DRIMSeqFiltering==TRUE){
  save(RATsObsRes, file = paste0(save_dir, "RATsObsResDRIMSeqFiltering.RData"))
}else if(DRIMSeqFiltering==FALSE){
  save(RATsObsRes, file = paste0(save_dir, "RATsObsRes.RData"))
}

#rm(list = ls())
# gc()

#source('~/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/DRIMSeq.R')



# 
# #Load tx2gene object
# load(paste0(def_wd, "tx2gene.RData"))
# 
# setwd(paste0(def_wd, sub_dir1))
# 
# QuantFiles <- mixedsort(list.files(pattern = ".sf", recursive = TRUE, full.names = TRUE))
# QuantDirsTemp <- mixedsort(list.dirs(recursive = FALSE, full.names = TRUE))
# 
# QuantDirs <- QuantDirsTemp[!QuantDirsTemp %in% "./GibbsSamps"]
# 
# if(QuantDirs[1]=="./GibbsSamps"){
#   QuantDirs[1] <- NA
# }
# names(QuantFiles) <- paste0("Sample", 1:length(QuantFiles))
# names(QuantDirs) <- paste0("Sample", 1:length(QuantDirs))
# 
# 
# #Create key matrix that contains matches samples to conditions and identifiers
# key <- matrix(c("SRR950078", "A",
#                 "SRR950080", "A",
#                 "SRR950082", "A",
#                 "SRR950084", "A",
#                 "SRR950086", "A",
#                 "SRR950079", "B",
#                 "SRR950081", "B",
#                 "SRR950083", "B",
#                 "SRR950085", "B",
#                 "SRR950087", "B"), nrow = 10, ncol = 2, byrow = TRUE)
# key <- as.data.frame(key, stringsAsFactors = FALSE)
# key[3] <- paste0("Sample",1:nrow(key))
# colnames(key) <- c("Sample", "Condition", "Identifier")
# key$Condition <- relevel(as.factor(key$Condition), ref = 1)
# 
# samples_A <- QuantDirs[1:5]
# samples_B <- QuantDirs[6:10]
# 
# #Change tx2gene column names to match what RATs default values are, I couldn't make it work inputting
#   #custom ones
# colnames(tx2gene) <- c("parent_id", "target_id", "NTrans")
# 
# 
# mydata <- fish4rodents(A_paths = samples_A, B_paths = samples_B, annot = tx2gene, scaleto = 100000000)
# 
# mydtu <- call_DTU(annot= tx2gene, boot_data_A= mydata$boot_data_A, 
#                   boot_data_B= mydata$boot_data_B, verbose= FALSE, 
#                   scaling=c(25, 26, 23, 50, 45, 48, 52))
