#Code to run the compositional analysis one gene at a time, used to get the computation time tables
  #Saving the gene-specific files for the dataset will be necessary before running this code
library(compositions)
library(plyr)
library(tidyr)
library(data.table)
library(gtools)
library(Matrix)

GibbsSamps <- FALSE

useOtherGroups <- TRUE
DRIMSeqFiltering <- TRUE
if(DRIMSeqFiltering==TRUE){
  useOtherGroups <- FALSE
}

dataset <- "GEUV1"

#Set to TRUE to calculate results for the 20 sample only analysis, and FALSE for the 100 sample analysis
TwentySamplesTotalAnalysis <- FALSE

#dataset <- "SQCC"

#The updated power files in the UpdatedPowerFilesDirec for GibbsSamps=TRUE use Gibbs with Thin100, but needs
#to just be called "Gibbs" here for code to work as currently implemented
#But, all the files referenced here for "Gibbs" are using a thinning value of 100
if(GibbsSamps==TRUE){
  infReps <- "Gibbs"
}else{
  infReps <- "Boot"
}


ninfreps <- 100



  if(dataset=="GEUV1"){
    def_wd <- "/Users/Scott/Documents/Dissertation Data/GEUV1Data/"
    wd2 <- "/Users/Scott/Documents/Dissertation Data/GEUV1Data/"
    if(TwentySamplesTotalAnalysis==TRUE){
      save_dir <- "/Users/Scott/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/GEUV1TwentySamples/"
    }else{
      save_dir <- "/Users/Scott/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/GEUV1/"
    }
    npartstouse <- 100
  }else if(dataset=="SQCC"){
    def_wd <- "/Users/Scott/Documents/Dissertation Data/SQCCDataReproduceOldResBeforeCommonCode/"
    wd2 <- "/Users/Scott/Documents/Dissertation Data/SQCCDataReproduceOldResBeforeCommonCode/"
    save_dir <- "/Users/Scott/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/SQCC/"
    npartstouse <- 10
  }
  source("/Users/Scott/Documents/Dissertation/code/CompFunctions.R")

setwd(wd2)

if(!dir.exists(save_dir)){dir.create(save_dir)}


if(DRIMSeqFiltering==TRUE & TwentySamplesTotalAnalysis==TRUE){
  direc_to_save <- paste0(def_wd, "GeneLevelFilesDRIMSeqFilteringTwentySamples/")
}else if(DRIMSeqFiltering==TRUE & TwentySamplesTotalAnalysis==FALSE){
  direc_to_save <- paste0(def_wd, "GeneLevelFilesDRIMSeqFiltering/")
}else if(DRIMSeqFiltering==FALSE & TwentySamplesTotalAnalysis==TRUE){
  direc_to_save <- paste0(def_wd, "GeneLevelFilesTwentySamples/")
}else if(DRIMSeqFiltering==FALSE & TwentySamplesTotalAnalysis==FALSE){
  direc_to_save <- paste0(def_wd, "GeneLevelFiles/")
}


if(DRIMSeqFiltering==TRUE){
  #load("abGeneFiltered.RData")
  #load("cntGenecntsScaledTPMFiltered.RData")
  load("abDatasetsNoOtherGroupsFiltered.RData")
  #load("tx2gene.RData")
  
  
  #abGene <- abGeneFiltered
  #cntGene <- cntGeneFiltered
  abDatasets <- abDatasetsFiltered
}else if(DRIMSeqFiltering==FALSE){
  #load("abGene.RData")
  #load("cntGenecntsScaledTPM.RData")
  load("abDatasets.RData")
  #load("tx2gene.RData")
}



genestouse <- names(abDatasets)



func1 <- function(x, load_dir, gfiles, infReps){
  
  curr_gene <- x
  fil_to_load <- paste0(curr_gene, ".RData")
  if(!(fil_to_load %in% gfiles)){
    return(NULL)
  }
  load(paste0(load_dir, curr_gene, ".RData"))
  sr <- proc.time()
  if(infReps==TRUE){
    # res <- CompositionalObsAnalysis2(abDatasetsToUse = abDatasetsToUse, fullgenenames = curr_gene, Group = Group,
    #                                  samps = samps, nsamp = length(samps), NewModelingFeb2019 = TRUE,
    #                                  newAbDatasetsGibbsFinal = newAbDatasetsGibbsFinal, ilrMeansCovs = ilrMeansCovs,
    #                                  ninfreps = ninfreps, useOtherGroups = FALSE)
    
    res <- CompositionalObsAnalysis2(abDatasetsToUse = NULL, fullgenenames = curr_gene, Group = Group,
                                     samps = samps, nsamp = length(samps), NewModelingFeb2019 = TRUE,
                                     newAbDatasetsGibbsFinal = NULL, ilrMeansCovs = NULL,
                                     ninfreps = ninfreps, useOtherGroups = FALSE, Y = Y, YInfRep = YInfRep, mean.withinhat = mean.withinhat)
  }else{
    # res <- CompositionalObsAnalysis2(abDatasetsToUse = abDatasetsToUse, fullgenenames = curr_gene, Group = Group,
    #                                  samps = samps, nsamp = length(samps), NewModelingFeb2019 = FALSE,
    #                                  newAbDatasetsGibbsFinal = NULL, ilrMeansCovs = NULL,
    #                                  ninfreps = ninfreps, useOtherGroups = FALSE, PillaiOnly = T)
    
    res <- CompositionalObsAnalysis2(abDatasetsToUse = NULL, fullgenenames = curr_gene, Group = Group,
                                     samps = samps, nsamp = length(samps), NewModelingFeb2019 = FALSE,
                                     newAbDatasetsGibbsFinal = NULL, ilrMeansCovs = NULL,
                                     ninfreps = ninfreps, useOtherGroups = FALSE, Y = Y, YInfRep = YInfRep, mean.withinhat = mean.withinhat)
  }


  comptgene <- proc.time() - sr
  #print(paste0("current gene is ", curr_gene))
  #print(comptgene)
  if(is.null(res)){
    return(res)
  }
  #res$CompPvals$comptgene <- comptgene[3]
  #return(res$CompPvals)
  res$comptgene <- comptgene[3]
  return(res)
}

load_dir <- direc_to_save
gfiles <- list.files(load_dir)
gc()
stt <- proc.time()
#c(genestouse[-1], genestouse[1])
t1 <- lapply(genestouse, func1, load_dir =  load_dir, gfiles = gfiles, infReps = F)
TotalTimeCompSplit <- proc.time() - stt
gcCompositionalResSplit <- gc()
CompositionalResSplit <- rbindlist(t1, fill = T)

if(DRIMSeqFiltering==TRUE){
  save(CompositionalResSplit, TotalTimeCompSplit, gcCompositionalResSplit, file = paste0(save_dir, "CompositionalResSplitDRIMSeqFiltering.RData"))
}else if(DRIMSeqFiltering==FALSE){
  save(CompositionalResSplit, TotalTimeCompSplit, gcCompositionalResSplit, file = paste0(save_dir, "CompositionalResSplit.RData"))
}





gc()
stt2 <- proc.time()
t5 <- lapply(genestouse, func1, load_dir =  load_dir, gfiles = gfiles, infReps = T)
#t5 <- laply(genestouse, func1, load_dir =  load_dir, infReps = T, .inform = T)
TotalTimeModelingSplit <- proc.time() - stt2
gcCompositionalModelingResSplit <- gc()
CompositionalModelingResSplit <- rbindlist(t5, fill = T)
if(DRIMSeqFiltering==TRUE){
  save(CompositionalModelingResSplit, TotalTimeModelingSplit, gcCompositionalModelingResSplit, file = paste0(save_dir, "CompositionalModelingResSplitDRIMSeqFiltering.RData"))
}else if(DRIMSeqFiltering==FALSE){
  save(CompositionalModelingResSplit, TotalTimeModelingSplit, gcCompositionalModelingResSplit, file = paste0(save_dir, "CompositionalModelingResSplit.RData"))
}







