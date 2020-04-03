#Summarize Power Results

#array_val should index over the change values
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if(version$nickname=="Planting of a Tree"){
  .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
}else{
  if(onCluster==TRUE){source("~/.Rprofile")}
}

Sys.sleep(sample(1:100,1))
onCluster <- TRUE #Change to True if running on the cluster, False if running on macbook

print(paste0("Current Slurm Job ID is ", Sys.getenv("SLURM_ARRAY_JOB_ID")))
#Parameter that controls whether other groups are used in the analysis or not
useOtherGroups <- FALSE

if(array_val <= 1000){
  TwentySamplesOnlyCase <- TRUE
  val_to_use <- array_val
}else if(array_val > 1000 & array_val <= 2000){
  TwentySamplesOnlyCase <- FALSE
  val_to_use <- array_val - 1000
}else{
  stop("array_val is misspecified")
}


DRIMSeqFiltering <- TRUE
if(DRIMSeqFiltering==TRUE){
  useOtherGroups <- FALSE
}

#Parameter that controls if Genes with Low MeanTGE (Total Gene Expression) across all samples for TPMs should be filtered and specifies what the value should be
#Can change the FilterCutoff Value to be more/less strict
#This is set to FALSE because genes are prefiltered by the DRIMSeqFiltering and this option is thus not needed, as DRIMSeq's filters
  #are stricter and preferred
FilterGenes <- FALSE
if(FilterGenes==TRUE){
  FilterCutoff <- 0.50
}



#The updated power files in the UpdatedPowerFilesDirec for GibbsSamps=TRUE use Gibbs with Thin100, but needs
#to just be called "Gibbs" here for code to work as currently implemented
#But, all the files referenced here for "Gibbs" are using a thinning value of 100
#infReps <- "Gibbs"
infReps <- "Boot"

#Specify whether all methods pvalues need to me non missing to in order to use that gene/group combination grouping
#Note that only the usual pillai MANOVA pvalue is used to determine if the Comp Method has non missing pvalues
#(ie the likelihood and likelihood with gibbs replicates are not used)
  #The code will have to be verified that it is working correctly if this is set to TRUE
AllMethodsNonMissingpVals <- FALSE



#Specify the number of parts the genes are split into as well as the number of group combinations for the permutation style analysis
#Since the results are nt saved in different parts, set nparts to NA to get it all to work right
nparts <- NA
ncombos <- 100

if(version$nickname=="Planting of a Tree"){
  .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
}else{
  if(onCluster==TRUE){source("~/.Rprofile")}
}
library(compositions)
library(plyr)
library(tidyr)
library(data.table)

data_name <- "GEUV1"
if(onCluster==TRUE){
  def_wd <- "~/res/GEUV1Data/"
  setwd(def_wd)
  func_loc <-  "~/code/CompFunctions.R"
  save_dir <- "/pine/scr/s/k/skvanbur/GEUV1/"
  pvals_save_dir <- "/pine/scr/s/k/skvanbur/GEUV1/"
}else{
  
}
load("tx2gene.RData")

source(func_loc)
#load("changes.RData")
#nchange <- length(changes)


nchange <- 1
# changestouse <- c(1,4)
# curr_change <- changestouse[val_to_use]


methds <- c("CompDTU", "CompDTUme", "DRIMSeq", "DRIMSeqAdd1ToEveryCount", "CompMI", "CompMICombineCoefs", "RATs", "RATsInfReps", "CompDTUAbRowMeanInfRepBoot")
#Note, to run the last CompDTUAbRowMeanInfRepBoot method, set DRIMSeqFiltering to TRUE and infReps to Boot because this method doesn't depend on infReps
  #and this way you will always know where the results will be saved

#Current change value for this run
if(val_to_use < 100){
  curr_change <- 1
}else if(val_to_use >=101 & val_to_use < 200){
  curr_change <- 4
}else if(val_to_use >=201 & val_to_use < 300){
  curr_change <- 8
}else if(val_to_use >=301 & val_to_use < 400){
  curr_change <- 2
}else if(val_to_use >=401 & val_to_use < 500){
  curr_change <- 1.5
}else if(val_to_use >=501 & val_to_use < 600){
  curr_change <- 16
}

val <- val_to_use %% 100

curr_methd <- methds[val]

if(curr_methd=="DRIMSeq" | curr_methd== "DRIMSeqAdd1ToEveryCount"){
  save_dir <- "~/res/GEUV1Data/"
}

#For these methods, set DRIMSeqFiltering to TRUE and infReps to Boot because these methods don't depend on infReps (since the full set of infReps is not used)
  #and this way you will always know where the results will be saved
if(curr_methd %in% c("CompDTUAbRowMeanInfRepBoot")){
  infReps <- "Boot"
  DRIMSeqFiltering <- TRUE
  useOtherGroups <- FALSE
}

#Get list of the genes that were changed from one of the results
  #The list differs between DRIMSeqFiltering and Other Groups since the list of final genes used in the analysis differs and
  #This ensures that half are changed and half are not in both cases
load("FilesForPowerAnalysis1.RData")

#There is a list of genestochange within FilesForPowerAnalysis1 but this does not correspond to the DRIMSeqFiltering
#If DRIMSeqFiltering is true make sure this is removed and that the proper one is loaded
if(DRIMSeqFiltering==TRUE){
  rm(genestochange)
  geneschanged <- loadRData(paste0(def_wd, "geneschangedDRIMSeqFilteringPowerAnalysis1.RData"))
}

print(paste0("Length of geneschanged is ", length(geneschanged)))



SummarizePowerRes(def_wd = def_wd, save_dir = save_dir, curr_change = curr_change, curr_methd = curr_methd,
                  geneschanged = geneschanged, data_name = data_name, nparts = nparts, 
                  ncombos = ncombos, FilterGenes = FilterGenes, FilterCutoff = FilterCutoff, 
                  useOtherGroups = useOtherGroups, AllMethodsNonMissingpVals = AllMethodsNonMissingpVals, 
                  tx2gene = tx2gene, pvals_save_dir = pvals_save_dir, cntsScaledTPM = TRUE,
                  infReps = infReps, DRIMSeqFiltering = DRIMSeqFiltering,
                  GEUV1Data = TRUE, useInfRV = TRUE, TwentySamplesOnlyCase = TwentySamplesOnlyCase)

print(gc())



