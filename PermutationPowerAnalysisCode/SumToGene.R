#Need to have run "SaveSalmonDataAsRData.R" for this to work
#Code to summarize the Salmon results to the gene level and create datasets and objects necessary 
  #for the Compositional analysis
  #Note that the regular quantification results don't differ depending on Gibbs/Bootsamps
  #But the results where the abundance values are computed using the mean/median of boot or gibbs samps do

array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
OnlyCalculateAbundanceFromInfReps <- TRUE
if(version$nickname=="Planting of a Tree"){
  .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
}else{
  if(onCluster==TRUE){source("~/.Rprofile")}
}

library(parallel)
library(plyr)
library(tidyr)
library(gtools)
library(Matrix)

onCluster <- TRUE #Change to True if running on the cluster, False if running on macbook

#Set countsFromAbundance parameter, which controls how the counts are estimated 
  #(see tximport help documentation for more information on the differences this parameter makes)

countsFromAbundance <- "no"
#countsFromAbundance <- "scaledTPM"

if(onCluster==TRUE){
  def_wd <- "~/res/GEUV1Data/"
  load_dir <- "/pine/scr/s/k/skvanbur/GEUV1/"
  setwd(load_dir)
  func_loc <-  "~/code/CompFunctions.R"
  clust <- makeCluster(1) # Change this to match whatever number of cores requested on cluster or macbook
  if(countsFromAbundance == "no"){
    load(paste0(load_dir,"Salmon/SalmonData.RData"))
  }else if(countsFromAbundance =="scaledTPM"){
    load(paste0(load_dir,"Salmon/SalmonDataCountsScaledTPM.RData"))
  }
}else{
  def_wd <- "/Users/Scott/Documents/Dissertation Data/GEUV1Data/"
  setwd(def_wd)
  func_loc <- "/Users/Scott/Documents/Dissertation/Code/CompFunctions.R"
  clust <- makeCluster(1) # Change this to match whatever number of cores requested on cluster or macbook
  if(countsFromAbundance =="no"){
    load("SalmonData.RData")
  }else if(countsFromAbundance =="scaledTPM"){
    load("SalmonDataCountsScaledTPM.RData")
  }

}

source(func_loc)
load(paste0(def_wd, "tx2gene.RData"))

#Export all loaded functions to the cluster object
clusterExport(clust, lsf.str())

if(OnlyCalculateAbundanceFromInfReps==FALSE){
  #Keep track of computing time
  startTime <- proc.time()
  
  sumToGene(QuantSalmon = QuantSalmon, tx2gene = tx2gene, clust = clust, key = key)
  
  proc.time() - startTime
}



rm(QuantSalmon)
gc()
print(gc())
#Now, summarize the results with abundances calculated also using the mean/medians of the bootstrap/gibbs samples also

#Results summarizing abundances using boot/gibbs samples has only been run for scaledTPM values

if(OnlyCalculateAbundanceFromInfReps==TRUE){
  #Summarize Gibbs Samples to gene
  if(array_val==1){
    if(countsFromAbundance=="scaledTPM"){
      QuantSalmon <- loadRData(paste0(load_dir,"Salmon/SalmonDataCountsScaledTPMAbRowMedianInfRep.RData"), objNameToGet = "QuantSalmonAbRowMedianInfRep")
    }else if(countsFromAbundance=="no"){
      QuantSalmon <- loadRData(paste0(load_dir,"Salmon/SalmonDataAbRowMedianInfRep.RData"), objNameToGet = "QuantSalmonAbRowMedianInfRep")
    }
    abFromInfRepFunc <- "median"
    GibbsSamps <- TRUE
    
    setwd(load_dir)
    sumToGene(QuantSalmon = QuantSalmon, tx2gene = tx2gene, clust = clust, key = key, abFromInfRepFunc = abFromInfRepFunc, GibbsSamps = GibbsSamps, countsFromAbundance = countsFromAbundance)
    
    rm(QuantSalmon)
    gc()
    print(gc())
  }

  #Summarize Gibbs Samples to gene
  if(array_val==2){
    if(countsFromAbundance=="scaledTPM"){
      QuantSalmon <- loadRData(paste0(load_dir,"Salmon/SalmonDataCountsScaledTPMAbRowMeanInfRep.RData"), objNameToGet = "QuantSalmonAbRowMeanInfRep")
    }else if(countsFromAbundance=="no"){
      QuantSalmon <- loadRData(paste0(load_dir,"Salmon/SalmonDataAbRowMeanInfRep.RData"), objNameToGet = "QuantSalmonAbRowMeanInfRep")
    }
    abFromInfRepFunc <- "mean"
    GibbsSamps <- TRUE
    
    setwd(load_dir)
    sumToGene(QuantSalmon = QuantSalmon, tx2gene = tx2gene, clust = clust, key = key, abFromInfRepFunc = abFromInfRepFunc, GibbsSamps = GibbsSamps, countsFromAbundance = countsFromAbundance)
    
    rm(QuantSalmon)
    gc()
    print(gc())
  }

  
  #Now, bootstrap samples
  if(array_val==3){
    if(countsFromAbundance=="scaledTPM"){
      QuantSalmon <- loadRData(paste0(load_dir,"SalmonBootSamps/SalmonDataCountsScaledTPMAbRowMedianInfRep.RData"), objNameToGet = "QuantSalmonAbRowMedianInfRep")
    }else if(countsFromAbundance=="no"){
      QuantSalmon <- loadRData(paste0(load_dir,"SalmonBootSamps/SalmonDataAbRowMedianInfRep.RData"), objNameToGet = "QuantSalmonAbRowMedianInfRep")
      
    }
    abFromInfRepFunc <- "median"
    GibbsSamps <- FALSE
    
    setwd(load_dir)
    sumToGene(QuantSalmon = QuantSalmon, tx2gene = tx2gene, clust = clust, key = key, abFromInfRepFunc = abFromInfRepFunc, GibbsSamps = GibbsSamps, countsFromAbundance = countsFromAbundance)
    
    rm(QuantSalmon)
    gc()
    print(gc())
  }

  #Now, bootstrap samples
  if(array_val==4){
    if(countsFromAbundance=="scaledTPM"){
      QuantSalmon <- loadRData(paste0(load_dir,"SalmonBootSamps/SalmonDataCountsScaledTPMAbRowMeanInfRep.RData"), objNameToGet = "QuantSalmonAbRowMeanInfRep")
    }else if(countsFromAbundance=="no"){
      QuantSalmon <- loadRData(paste0(load_dir,"SalmonBootSamps/SalmonDataAbRowMeanInfRep.RData"), objNameToGet = "QuantSalmonAbRowMeanInfRep")
    }
    abFromInfRepFunc <- "mean"
    GibbsSamps <- FALSE
    
    setwd(load_dir)
    sumToGene(QuantSalmon = QuantSalmon, tx2gene = tx2gene, clust = clust, key = key, abFromInfRepFunc = abFromInfRepFunc, GibbsSamps = GibbsSamps, countsFromAbundance = countsFromAbundance)
    
    rm(QuantSalmon)
    gc()
    print(gc())
  }

}







