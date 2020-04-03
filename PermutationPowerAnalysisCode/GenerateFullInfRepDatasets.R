#Need to have run SaveGibbsAsRData to be able to load the necessary files, which have Gibbs sample information for each biological sample

Sys.sleep(sample(1:1000,1))

if(version$nickname=="Planting of a Tree"){
  .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
}else{
  if(onCluster==TRUE){source("~/.Rprofile")}
}

library(data.table)
library(plyr)
library(dplyr)
library(parallel)
library(gtools)

startTime <- proc.time()
#The number of array values needs to match the number of chunks the genestouse are split into (100 being the default)
  #ie the number of parts
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

onCluster <- TRUE #Change to True if running on the cluster, False if running on macbook

#Set to TRUE for Gibbs Samps, FALSE for Boot Samps
GibbsSamps <- FALSE

DRIMSeqFiltering <- TRUE


dir1 <- "~/res/GEUV1Data/"
func_loc <-  "~/code/CompFunctions.R"
load(paste0(dir1,"tx2gene.RData"))
save_dir <- "/pine/scr/s/k/skvanbur/GEUV1/"
if(GibbsSamps==TRUE){
    def_wd1 <- "/pine/scr/s/k/skvanbur/GEUV1/Salmon/"
    sd1 <- "GibbsabDatasets/"
    sd2 <- "GibbscntDatasets/"
    sd3 <- "GibbsFullinfRepDat/"
    type <- "GibbsReps"
    load(paste0(dir1,"Salmon/SalmonData.RData"))
}else{
    def_wd1 <- "/pine/scr/s/k/skvanbur/GEUV1/SalmonBootSamps/"
    sd1 <- "BootabDatasets/"
    sd2 <- "BootcntDatasets/"
    sd3 <- "BootFullinfRepDat/"
    type <- "BootSamps"
    load(paste0(dir1,"SalmonBootSamps/SalmonDataCountsScaledTPM.RData"))
}
  

#abDatasets are only needed to generate the (observed) bootstrap datasets, so can load filtered ones here  
if(DRIMSeqFiltering==TRUE){
  load(paste0(dir1,"tx2gene.RData"))
  if(countsFromAbundance=="scaledTPM"){
    #load(paste0(dir1, "cntGenecntsScaledTPMFiltered.RData"))
    load(paste0(dir1, "cntGenecntsScaledTPM.RData"))
  }else if(countsFromAbundance=="no"){
    #load(paste0(dir1, "cntGeneFiltered.RData"))
    load(paste0(dir1, "cntGene.RData"))
  }
  
  #load(paste0(dir1, "abDatasets.RData"))
  load(paste0(dir1, "abDatasetsNoOtherGroupsFiltered.RData"))
  
  filteredgenenames <- names(abDatasetsFiltered)
  #cntGene <- cntGeneFiltered
  #load(paste0(dir1, "abDatasetsNoOtherGroups.RData"))
  
}else if(DRIMSeqFiltering==FALSE){
  load(paste0(dir1,"tx2gene.RData"))
  load(paste0(dir1, "abDatasets.RData"))
  load(paste0(dir1, "abDatasetsNoOtherGroups.RData"))
  
  if(countsFromAbundance=="scaledTPM"){
    load(paste0(dir1, "cntGenecntsScaledTPM.RData"))
  }else if(countsFromAbundance=="no"){
    load(paste0(dir1, "cntGene.RData"))
  }
  filteredgenenames <- names(abDatasets)
}

#clust <- makeCluster(12) # Change this to match whatever number of cores requested on cluster or are on macbook

source(func_loc)  


setwd(def_wd1)

nparts <- 100

curr_part_num <- array_val %% nparts
if(curr_part_num==0){
  curr_part_num <- nparts
}


SaveFullinfRepDat(CreateabAndcntDatasets = FALSE)

proc.time() - startTime
print(gc())




