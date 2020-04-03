#Save Gibbs or bootstrap data as an initial RData file for each sample that can later be loaded and manipulated to create abDatasets, etc.
  #Need to this file separately for each Biological sample differently because otherwise the file is too big
  #So, run as an array job indexing over all biological samples

Sys.sleep(sample(1:500,1))

if(version$nickname=="Planting of a Tree"){
  .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
}else{
  if(onCluster==TRUE){source("~/.Rprofile")}
}

library(data.table)
library(dplyr)
library(tximport)
library(rjson)
library(gtools)

onCluster <- TRUE #Change to True if running on the cluster, False if running on macbook

#Set to True if you've drawn Gibbs samples and FALSE is you've drawn bootstrap samples
GibbsSamps <- FALSE

#Note that countsFromAbundance is currently ignored because the Gibbs counts don't change using tximport, and this means
  #The resulting TPMs won't change either
#Set countsFromAbundance parameter, which controls how the counts are estimated 
  #(see tximport help documentation for more information on the differences this parameter makes)
  #countsFromAbundance <- "no"
  countsFromAbundance <- "scaledTPM"


if(onCluster==TRUE){
  if(GibbsSamps==TRUE){
    def_wd1 <- "/pine/scr/s/k/skvanbur/GEUV1/Salmon/"
  }else{
    def_wd1 <- "/pine/scr/s/k/skvanbur/GEUV1/SalmonBootSamps/"
  }
  
  dir1 <- "~/res/GEUV1Data/"
  func_loc <-  "~/code/CompFunctions.R"
  setwd(def_wd1)
  
  if(countsFromAbundance == "no" & GibbsSamps==TRUE){
    load(paste0(dir1,"Salmon/SalmonData.RData"))
  }else if(countsFromAbundance == "no" & GibbsSamps==FALSE){
    load(paste0(dir1,"SalmonBootSamps/SalmonData.RData"))
  }else if(countsFromAbundance =="scaledTPM" & GibbsSamps==TRUE){
    load(paste0(dir1,"Salmon/SalmonDataCountsScaledTPM.RData"))
  }else if(countsFromAbundance =="scaledTPM" & GibbsSamps==FALSE){
    load(paste0(dir1,"SalmonBootSamps/SalmonDataCountsScaledTPM.RData"))
  }
  
  
  load(paste0(dir1,"tx2gene.RData"))
  #clust <- makeCluster(12) # Change this to match whatever number of cores requested on cluster or are on macbook
}else{
#Probably won't ever run this on the macbook anyway
}
source(func_loc)

#val indexes the biological samples
val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

#Extract file location of all of the quantification files of interest, the .sf files
#Also need to carefully name each QuantFile corresponding to the identifiers created above
QuantFiles <- mixedsort(list.files(pattern = ".sf", recursive = TRUE, full.names = TRUE))
subnames <- dir(recursive=F, pattern = "ERR")
for(i in 1:length(subnames)){
  currname <- subnames[i]
  wherefound <- grep(currname, QuantFiles)
  names(QuantFiles)[wherefound] <- key$Identifier[key$ENARun==currname]
}

QuantFiles2 <- QuantFiles[mixedsort(names(QuantFiles))]

curr_samp <- names(QuantFiles2)[val]
curr_file_loc <- QuantFiles2[val]

startTime <- proc.time()
SaveGibbsDataAsRData(curr_samp = curr_samp, curr_file_loc = curr_file_loc, GibbsSamps = GibbsSamps)
proc.time() - startTime


  

