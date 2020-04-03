#Filter the genes for the analysis and save a list of all genes to use
  #Filters here based on the filters that DRIMSeq has (which are based on counts, not TPM)
  #See Mike's paper Swimming downstream: statistical analysis of differential transcript usage following Salmon quantification
  #Note that the DRIMSeq vignette only implements 2 of the 3 at once, while Mike's paper uses all 3 so we will use all 3 here
StartTime <- proc.time()
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if(version$nickname=="Planting of a Tree"){
  .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
}else{
  if(onCluster==TRUE){source("~/.Rprofile")}
}

library(data.table)
library(plyr)
library(parallel)
library(gtools)
library(compositions)
library(DRIMSeq)

onCluster <- TRUE #Change to True if running on the cluster, False if running on macbook



if(onCluster==TRUE){
  def_wd1 <- "/pine/scr/s/k/skvanbur/GEUV1/Salmon/"
  save_dir <- "/pine/scr/s/k/skvanbur/GEUV1/"
  dir1 <- "~/res/GEUV1Data/"
  func_loc <-  "~/code/CompFunctions.R"
  setwd(dir1)
}else{
  def_wd <- "/Users/Scott/Documents/Dissertation Data/GEUV1Data/"
  def_wd2 <- "/Users/Scott/Documents/Dissertation/res/GEUV1Data"
  setwd(def_wd)
  func_loc <- "/Users/Scott/Documents/Dissertation/code/CompFunctions.R"
}

load("abGene.RData")
load("abDatasets.RData")
load("tx2gene.RData")

countsFromAbundance <- "scaledTPM"
#countsFromAbundance <- "no"
source(func_loc)

if(countsFromAbundance=="scaledTPM"){
  load("cntGenecntsScaledTPM.RData")
  load("cntDatasetsNoOtherGroupscntsScaledTPM.RData")
  load("failedgibbssampsCountsScaledTPM.RData")
}else if(countsFromAbundance=="no"){
  load("cntGene.RData")
  load("cntDatasetsNoOtherGroups.RData")
  load("failedgibbssamps.RData")
}





#cntGeneToUse <- cntGene

#Code to use DRIMSeq to get the list of transcripts/genes that pass the filtering criterion
#Set parameters used by DRIMSeq for its filtering
n <- nrow(key)
n.small <- min(table(key$Condition))
min_samps_feature_expr=n.small
min_feature_expr=10
min_samps_feature_prop=n.small
min_feature_prop=0.1
min_samps_gene_expr=n
min_gene_expr=10
sampstouse <- key$Identifier

if(!(array_val %in% c(1,2,3,4))){
  DRIMSeqFilter(abGene = abGene, cntGene = cntGene)
  
  total_t <- proc.time() - StartTime
  print(paste0("Total time needed to complete DRIMSeq Filtering is ", total_t))
}
#DRIMFiltering()


rm(cntGene)
setwd(save_dir)


if(array_val==1){
  runDRIMSeqFilterOnAbFromInfReps(abFromInfRepFuncT = "median", GibbsSampsT = TRUE, abGeneRegularAbValues = abGene)
}else if(array_val==2){
  runDRIMSeqFilterOnAbFromInfReps(abFromInfRepFuncT = "mean", GibbsSampsT = TRUE, abGeneRegularAbValues = abGene)
}else if(array_val==3){
  runDRIMSeqFilterOnAbFromInfReps(abFromInfRepFuncT = "median", GibbsSampsT = FALSE, abGeneRegularAbValues = abGene)
}else if(array_val==4){
  runDRIMSeqFilterOnAbFromInfReps(abFromInfRepFuncT = "mean", GibbsSampsT = FALSE, abGeneRegularAbValues = abGene)
}









