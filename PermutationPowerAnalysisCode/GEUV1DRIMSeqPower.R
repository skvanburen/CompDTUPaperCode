#DRIMSeq Code to run power analysis changing the factor that is multiplied by the counts
Sys.sleep(sample(1:1000,1))


val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if(version$nickname=="Planting of a Tree"){
  .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
}else{
  #source("~/.RProfile")
  .libPaths("/nas/longleaf/home/skvanbur/bin/R")
}
#library(parallel)
library(DRIMSeq)
library(data.table)

#useOtherGroups <- TRUE
DRIMSeqFiltering <- TRUE
if(DRIMSeqFiltering==TRUE){
  useOtherGroups <- FALSE
}else{
  useOtherGroups <- TRUE
}

#Add 1 to every count in DRIMSeq to greatly stabilize results
Add1ToEveryCount <- TRUE

TwentySamplesTotalAnalysis <- FALSE

onCluster <- TRUE #Change to True if running on the cluster, False if running on macbook

if(onCluster==TRUE){
  def_wd <- "~/res/GEUV1Data/"
  setwd(def_wd)
  ncores <- 1
  #clust <- makeCluster(ncores) # Change this to match whatever number of cores requested on cluster or are on macbook
  func_loc <-  "~/code/CompFunctions.R"
  source(func_loc)
}else{
  def_wd <- "/Users/Scott/Documents/Dissertation Data/GEUV1Data"
  setwd(def_wd)
  
  func_loc <-  "/Users/Scott/Documents/Dissertation/code/CompFunctions.R"
  source(func_loc)
}


#load("AllGroupCombinations.RData")

if(DRIMSeqFiltering==TRUE){
  #cntGenecntsScaledTPMFiltering is only used to get the transcripts/genes that pass filtering
  load("cntGenecntsScaledTPMFiltered.RData")
  load("cntGene.RData")
  load("abGene.RData")
  load("abDatasetsNoOtherGroupsFiltered.RData")
  load("cntDatasetscntsScaledTPMNoOtherGroupsFiltered.RData")
  #load("changes.RData")
  load("tx2gene.RData")
  
  abDatasets <- abDatasetsFiltered
  cntDatasets <- cntDatasetsFiltered
}else if(DRIMSeqFiltering==FALSE){
  load("cntGene.RData")
  load("abGene.RData")
  load("abDatasets.RData")
  load("cntDatasets.RData")
  #load("changes.RData")
  load("tx2gene.RData")
}


#Export all loaded functions to the cluster object
#clusterExport(clust, lsf.str())

#Load DRIMSeq package to the cluster since it is not in the default location
# if(onCluster==TRUE){
#   clusterEvalQ(clust, source("~/.Rprofile"))
# }
# clusterEvalQ(clust, library(DRIMSeq))

#ChangeNum <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
ncombos <- 100
GroupNum <- val %% ncombos
if(GroupNum==0){
  GroupNum <- ncombos
}

if(val <= ncombos){
  change <- 4
}else if(val >=ncombos + 1 & val <=2 * ncombos){
  change <- 1 
}else if(val >=(2 * ncombos) + 1 & val <= (3 * ncombos)){
  change <- 2
}else if(val >=(3 * ncombos) + 1 & val <= (4 * ncombos)){
  change <- 8
}else if(val >=(4 * ncombos) + 1 & val <= (5 * ncombos)){
  change <- 1.5
}else if(val >=(5 * ncombos) + 1 & val <=(6 * ncombos)){
  change <- 16
}

if(TwentySamplesTotalAnalysis==TRUE){
  rm(Group)
  load("grpcombosAndsub_keyDRIMSeqFilteringPowerAnalysis1TwentySamples.RData")
  grpcombos <- grpcombos_TwentySamples
}else{
  rm(Group)
  load("DRIMSeqFilesForPowerAnalysis1.RData")
}








Group <- as.factor(as.character(grpcombos[GroupNum,]))

curr_key <- data.frame(sub_key$Sample, Group, sub_key$Identifier, stringsAsFactors = F)
colnames(curr_key) <- c("Sample", "Condition", "Identifier")

curr_key$Condition <- relevel(as.factor(curr_key$Condition), ref = 1)


curr_key_to_use <- curr_key

genestochange <- loadRData("geneschangedDRIMSeqFilteringPowerAnalysis1.RData")

# if(TwentySamplesTotalAnalysis==TRUE){
#   #Take the first 10 samples in each condition from the 100 samples used from power analysis 1 based on the true condition assignment and only use data from those 
#   #s_c1 and s_c2 will be the sample numbers for the first and second conditions respectively
#   rm(Group)
#   rm(curr_key)
#   
#   
#   
#   s_c1 <- as.numeric(rownames(sub_key)[sub_key$Condition=="CEU"][1:10])
#   s_c2 <- as.numeric(rownames(sub_key)[sub_key$Condition=="GBR"][1:10])
#   
#   #Order the sample numbers from smallest to largest
#   s_ordered <- sort(c(s_c1, s_c2))
#   
#   #samps <- mixedsort(c(paste0("Sample", s_c1), paste0("Sample", s_c2)))
#   samps <- paste0("Sample", s_ordered)
#   
#   #sub_key2 <- subset(sub_key, sub_key$Identifier %in% samps)
#   #curr_cond_temp <- sub_key2$Condition
#   
#   curr_cond <- factor(rep(NA, 20), levels = levels(sub_key$Condition))
#   
#   #Set this seed for the generation of the group combinations for the 20 sample case to ensure the group combinations are the same as RATs and the CompDTU methods
#     #Need to generate a new group combination because simply subsetting from the existing one would not ensure each condition had 10 samples exactly
#   set.seed(53 * GroupNum)
#   
#   pos <- sample(1:20, size = 10, replace = FALSE)
#   curr_cond[pos] <- levels(curr_cond)[1]
#   curr_cond[-pos] <- levels(curr_cond)[2]
#   
#   curr_key_to_use <- subset(curr_key, curr_key$Identifier %in% samps)
#   
#   #Update the group variables used in the analysis to correspond to the new randomly generated condition variable
#   curr_key_to_use$Condition <- curr_cond
#   Group <- curr_key_to_use$Condition
#   print(Group)
# }else{
#   
#   curr_key_to_use <- curr_key
# }


startTime <- proc.time() 
#Must set a seed for a power analysis to ensure each Gibbs replicate has the same list of genes changed to be able to compute ROC curve properly
assign(paste0("DRIMSeqROCResGroupCombo", GroupNum), DRIMSeqPower(y = Group, change = change, cntGene = cntGene, abGene = abGene, cntGeneFiltered = cntGeneFiltered, 
                                                                 cntDatasets = cntDatasets, key = curr_key_to_use, seed = NULL, tx2gene = tx2gene,
                                                                 genestochange = genestochange, useOtherGroups = useOtherGroups,
                                                                 Add1ToEveryCount = Add1ToEveryCount))
assign(paste0("DRIMSeqROCResGroupCombo", GroupNum, "Time"), proc.time() - startTime)
name <- paste0("DRIMSeqROCResGroupCombo", GroupNum)
name2 <- paste0("DRIMSeqROCResGroupCombo", GroupNum, "Time")

if(TwentySamplesTotalAnalysis==TRUE){
  direc_modifier <- "TwentySamples"
}else{
  direc_modifier <- ""
}

if(DRIMSeqFiltering==TRUE & Add1ToEveryCount==TRUE){
  direc <- paste0(def_wd, "DRIMSeqROCResDRIMSeqFilteringAdd1ToEveryCount", direc_modifier, "/", "Change", change, "/")
}else if(DRIMSeqFiltering==TRUE & Add1ToEveryCount==FALSE){
  direc <- paste0(def_wd, "DRIMSeqROCResDRIMSeqFiltering", direc_modifier, "/", "Change", change, "/")
}else if(DRIMSeqFiltering==FALSE & Add1ToEveryCount==TRUE){
  direc <- paste0(def_wd, "DRIMSeqROCResAdd1ToEveryCount", direc_modifier, "/", "Change", change, "/")
}else if(DRIMSeqFiltering==FALSE & Add1ToEveryCount==FALSE){
  direc <- paste0(def_wd, "DRIMSeqROCRes", direc_modifier, "/", "Change", change, "/")
  
}

if(!dir.exists(direc)) {dir.create(direc, recursive = TRUE)}
save(list = c(name, name2), file = paste0(direc, "DRIMSeqROCResGroupCombo", GroupNum, ".RData"))

print(gc())

#So, given above comment, just run it on all changes instead
#load("changes.RData")
#changes <- as.numeric(c(0.1, 0.25))

#This will run analysis for observed group arrangement for all changes
#parApply(clust, as.matrix(changes), 1, DRIMSeqPower2, key = key)
