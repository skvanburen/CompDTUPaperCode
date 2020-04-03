
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if(array_val <= 1000){
  TwentySamplesTotalAnalysis <- FALSE
  val <- array_val
}else if(array_val > 1000 & array_val <= 2000){
  TwentySamplesTotalAnalysis <- TRUE
  val <- array_val - 1000
}

ModifyAbundancesFromInfRepsForPowerAnalysis <- TRUE

onCluster <- TRUE #Change to True if running on the cluster, False if running on macbook
#Need to source .Rprofile on cluster before running to code so it loads other install directories for compositions packages
#Add the proper libpaths to R so it knows where to load MVN package from
if(version$nickname=="Planting of a Tree"){
  .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
}else{
  if(onCluster==TRUE){source("~/.Rprofile")}
}

library(compositions)
library(parallel)
library(plyr)
library(tidyr)
library(data.table)
library(gtools)

#Parameter that controls whether other groups are used in the analysis or not
DRIMSeqFiltering <- TRUE
useOtherGroups <- FALSE
#useilrMeansCovs <- TRUE

infReps <- "GibbsThin100"
#infReps <- "Boot"

if(DRIMSeqFiltering==TRUE){
  useOtherGroups <- FALSE
}

if(infReps=="GibbsThin100"){
  load_dir <- "/pine/scr/s/k/skvanbur/GEUV1/GibbsFullinfRepDat/"
}else if(infReps=="Boot"){
  load_dir <- "/pine/scr/s/k/skvanbur/GEUV1/BootFullinfRepDat/"
}

ninfreps <- 100
nparts <- 100
if(onCluster==TRUE){
  def_wd <- "~/res/GEUV1Data/"
  setwd(def_wd)
  ncores <- 1
  #clust <- makeCluster(ncores) # Change this to match whatever number of cores requested on cluster or are on macbook
  func_loc <-  "~/code/CompFunctions.R"
  wd2 <- "/pine/scr/s/k/skvanbur/GEUV1/"
  source(func_loc)
}else{
}

#Load the group combinations for this power analysis
#There is a list of genestochange within FilesForPowerAnalysis1 but this does not correspond to the DRIMSeqFiltering
  #IF DRIMSeqFiltering is true make sure this is removed and that will be recalculated below
load("FilesForPowerAnalysis1.RData")

#Load the changes (quantities the major transcript counts are multiplied by)
#load("changes.RData") 
#These counts to be modified for the power analysis should not use the cntsScaledTPM values because the regular, counts unscaled counts have to be updated
  #and converted to TPMs and use of the scaledTPM counts to start with would result in an incorrect total library size, giving incorrect TPM values
  #The object cntGenecntsScaledTPMFiltered is only loaded to get the list of genes/transcripts that pass filtering, and
  #The counts that will be updated come from the usual, unscaled counts located in cntGene.RData
if(DRIMSeqFiltering==TRUE){
  load("abGene.RData")
  load("cntGenecntsScaledTPMFiltered.RData")
  load("cntGene.RData")
  load("abDatasetsNoOtherGroupsFiltered.RData")
  load("tx2gene.RData")
}else{
  load("abGene.RData")
  load("cntGenecntsScaledTPM.RData")
  load("abDatasets.RData")
  load("tx2gene.RData")
}

GroupNum <- val %% 100
if(GroupNum==0){
  GroupNum <- 100
}
#ChangeNum <- floor((val-1)/252) + 1

if(TwentySamplesTotalAnalysis==TRUE){
  rm(grpcombos)
  rm(sub_key)
  load("grpcombosAndsub_keyDRIMSeqFilteringPowerAnalysis1TwentySamples.RData")
  grpcombos <- grpcombos_TwentySamples
  
}

Group <- as.factor(as.character(grpcombos[GroupNum,]))
y <- Group




#This line is very important.  sub_key has the real observed group combination but needs to be changed for each 
  #separate groupcombo that the updated datasets are to be generated for
sub_key$Condition <- y


if(DRIMSeqFiltering==TRUE){
  rm(genestochange)
}

if(TwentySamplesTotalAnalysis==TRUE & length(y)!=20){
  stop("Something is wrong with the group specification")
}

#Use this to have change values from 0.1 to 100, a total of 10 change values
# So, want want array matrix to be indexed from 1 to 9(need array job to be 1:2520)
#change <- changes[ChangeNum + 2]



if(val <=100){
  change <- 4
}else if(val >100 & val <=200){
  change <- 1 
}else if(val >200 & val <=300){
  change <- 2 
}else if(val >300 & val <=400){
  change <- 8 
}



seedtoset <- 3776
half <- TRUE
#Set genestochange to null for now then generate this based on the seed
#Set to null first just to match the previous code that worked as closely as possible
genestochange <- NULL

if(!is.null(seedtoset)){
  set.seed(seedtoset)
}

#Choose genes that are changed in the power analysis based only on the genes that pass filtering, ie the ones that are in cntGeneFiltered
CompCntGene <- subset(cntGeneFiltered, (cntGeneFiltered$NTrans!=1 & cntGeneFiltered$SumTGE!=0))
fullgenenames <- sort(unique(CompCntGene$gene_id))


#Only change for half of the genes to be able to generate a ROC curve
if((half==TRUE & is.null(genestochange))){
  genestochange <- sort(sample(fullgenenames, ceiling(length(fullgenenames)/2)))
}

if(val==1 & DRIMSeqFiltering==TRUE & !file.exists(paste0(def_wd, "geneschangedDRIMSeqFilteringPowerAnalysis1.RData"))){
  save(genestochange, file = paste0(def_wd, "geneschangedDRIMSeqFilteringPowerAnalysis1.RData"), version = 2)
}

#Set generateNewIlrMeansCovs to false because the newer structure has these being calculated as the CompDTUme results are being run, not before
updateilrMeansCovsAndabDatasets(subset_data = TRUE, generateNewIlrMeansCovs = FALSE, InfRepSpecificPartFullinfRepDatFilesExist = TRUE,
                                TwentySamplesTotalAnalysis = TwentySamplesTotalAnalysis, ModifyAbundancesFromInfRepsForPowerAnalysis = ModifyAbundancesFromInfRepsForPowerAnalysis,
                                AbundancesFromInfRepsDir = wd2)

#This code can help you check after the jobs are all done if the the files all saved successfully
wd2 <- "/pine/scr/s/k/skvanbur/GEUV1/"
ModifyAbundancesFromInfRepsForPowerAnalysis <- TRUE
for(i in 1:300){
  #print(paste0("i is ", i))
  GroupNum <- i %% 100
  if(GroupNum==0){
    GroupNum <- 100
  }

  if(i <=100){
    change <- 4
  }else if(i >100 & i <=200){
    change <- 1
  }else if(i >200 & i <=300){
    change <- 2
  }else if(i >300 & i <=400){
    change <- 8
  }
  if(ModifyAbundancesFromInfRepsForPowerAnalysis==TRUE){
    curr_direc <- paste0(wd2, "UpdatedAbDatasetsAbFromInfRepsDRIMSeqFiltering/Change", change, "/", "GroupCombo", GroupNum, "/")
    #curr_direc <- paste0(wd2, "UpdatedAbDatasetsAbFromInfRepsDRIMSeqFilteringTwentySamples/Change", change, "/", "GroupCombo", GroupNum, "/")
    fil1 <- paste0(curr_direc, "UpdatedabDatasetsChange", change, "GroupCombo", GroupNum, "AbRowMedianInfRepBoot",  ".RData")
    fil2 <- paste0(curr_direc, "UpdatedabDatasetsChange", change, "GroupCombo", GroupNum, "AbRowMeanInfRepBoot",  ".RData")
    fil3 <- paste0(curr_direc, "UpdatedabDatasetsChange", change, "GroupCombo", GroupNum, "AbRowMedianInfRepGibbs",  ".RData")
    fil4 <- paste0(curr_direc, "UpdatedabDatasetsChange", change, "GroupCombo", GroupNum, "AbRowMeanInfRepGibbs",  ".RData")

    if(!file.exists(fil1) | !file.exists(fil2) | !file.exists(fil3) | !file.exists(fil4)){
      print(paste0("Job did not work properly for GroupCombo Num ", GroupNum, " and Change ", change))
    }

  }else{
    #curr_direc <- paste0(wd2, "UpdatedPowerData/Change", change, "/", "GroupCombo", GroupNum, "/")
    #curr_direc <- paste0(wd2, "UpdatedPowerDataBootDRIMSeqFiltering/Change", change, "/", "GroupCombo", GroupNum, "/")
    #curr_direc <- paste0(wd2, "UpdatedPowerDataBootDRIMSeqFilteringTwentySamples/Change", change, "/", "GroupCombo", GroupNum, "/")
    curr_direc <- paste0(wd2, "UpdatedAbDatasetsAbFromInfRepsDRIMSeqFiltering/Change", change, "/", "GroupCombo", GroupNum, "/")
    fil1 <- paste0(curr_direc, "UpdatedabDatasetsChange", change, "GroupCombo", GroupNum, ".RData")
    if(!file.exists(fil1)){
      print(paste0("Job did not work properly for GroupCombo Num ", GroupNum, " and Change ", change))
    }
    for(part in 1:nparts){
      fil2 <- paste0(curr_direc, "NewabDatasetsGibbsChange", change, "GroupCombo", GroupNum, "Part", part, ".RData")
      ###fil3 <- paste0(curr_direc, "NewilrMeansCovsChange", change, "GroupCombo", GroupNum, "Part", part, ".RData")
      if(!file.exists(fil2)){
        print(paste0("Job did not work properly for GroupCombo Num ", GroupNum, " and Change ", change))
      }
      ## if((!file.exists(fil2) | !file.exists(fil3))){
      ##   print(paste0("Job did not work properly for GroupCombo Num ", GroupNum))
      ## }
    }
  }
  rm(change)
}
  
  





