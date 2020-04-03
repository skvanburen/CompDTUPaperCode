
library(DRIMSeq)
library(data.table)

func_loc <- "/Users/Scott/Documents/Dissertation/Code/CompFunctions.R"
source(func_loc)

ninfreps <- 100

#dataset <- "SQCC"
dataset <- "GEUV1"

TwentySamplesTotalAnalysis <- FALSE



useOtherGroups <- FALSE
DRIMSeqFiltering <- TRUE
if(DRIMSeqFiltering==TRUE){
  useOtherGroups <- FALSE
}

if(dataset=="GEUV1"){
  def_wd <- "/Users/Scott/Documents/Dissertation Data/GEUV1Data/"
  if(TwentySamplesTotalAnalysis==TRUE){
    save_dir <- "/Users/Scott/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/GEUV1TwentySamples/"
  }else{
    save_dir <- "/Users/Scott/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/GEUV1/"
  }

  npartstouse <- 100
}else if(dataset=="SQCC"){
  def_wd <- "/Users/Scott/Documents/Dissertation Data/SQCCDataReproduceOldResBeforeCommonCode/"
  save_dir <- "/Users/Scott/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/SQCC/"
  npartstouse <- 10
}
setwd(def_wd)

  #npartstouse <- 5
  # #Only load the ilrMeansCovs to get a list of genes to use that matches the Compositional Analysis
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
  
  if(DRIMSeqFiltering==TRUE){
    load("abDatasetsNoOtherGroupsFiltered.RData")
    abDatasets <- abDatasetsFiltered
    FilteringMethod <- "DRIMSeq"
  }else if(DRIMSeqFiltering==FALSE){
    load("abDatasets.RData")
    FilteringMethod <- "OtherGroups"
  }
  #genestouseT <- names(ilrMeansCovs)
  #genestouseT2 <- names(abDatasets)
  
  #genestouse <- intersect(genestouseT, genestouseT2)[1:500]

  genestouse <- names(abDatasets)


#Set countsFromAbundance parameter, which controls how the counts are estimated 
#(see tximport help documentation for more information on the differences this parameter makes)
  countsFromAbundance <- "scaledTPM"
  
  if(TwentySamplesTotalAnalysis==TRUE){
    load("/Users/Scott/Documents/Dissertation Data/GEUV1Data/FilesForPowerAnalysis1.RData")
    
    s_c1 <- as.numeric(rownames(sub_key)[sub_key$Condition=="CEU"][1:10])
    s_c2 <- as.numeric(rownames(sub_key)[sub_key$Condition=="GBR"][1:10])
    
    #Order the sample numbers from smallest to largest
    s_ordered <- sort(c(s_c1, s_c2))
    
    #sampstouse <- c(paste0("Sample", s_c1), paste0("Sample", s_c2))
    sampstouse <- paste0("Sample", s_ordered)
    
    
    
    #key_to_use <- subset(key, key$Identifier %in% sampstouse)
    
    #Remove the unneeded factor levels
    #key_to_use$Condition <- factor(key_to_use$Condition)
    #nsamp <- length(sampstouse)
    
    #Group <- key_to_use$Condition
  }else{
    sampstouse <- key$Identifier
  }


  #Only dmPrecision, dmFit, and dmTest statements are included in the computation time that is used to construct the table
    #Not any time to create/load in or reformat data, etc to keep this comparison fair between methods
  
  DRIMSeqObservedAnalysis(countsFromAbundance = countsFromAbundance, FilteringMethod = FilteringMethod, sampstouse = sampstouse,
                          save_dir = save_dir, genestouse = genestouse, fullReturn = TRUE)
  
  