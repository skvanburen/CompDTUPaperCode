Sys.sleep(sample(1:1000,1))

if(version$nickname=="Planting of a Tree"){
  .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
}else{
  #source("~/.RProfile")
  .libPaths("/nas/longleaf/home/skvanbur/bin/R")
}
library(gtools)
library(rats)
library(data.table)
library(tidyr)

#infReps <- "GibbsThin100"
infReps <- "Boot"

#Filter transcripts like is done for the DRIMSeq and CompDTU analyses to make comparisons easier and more fair
  #This is despite the fact that the RATs draft does not directly do this
DRIMSeqFiltering <- TRUE

TwentySamplesTotalAnalysis <- FALSE

#val is the current group combination number
val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

ncombos <- 100
GroupNum <- val %% ncombos
if(GroupNum==0){
  GroupNum <- ncombos
}

if(val <= ncombos){
  change <- 4
}else if(val >=ncombos + 1 & val <=2 * ncombos){
  change <- 1 
}else if(val >=(2 * ncombos) + 1 & val <= 3 * ncombos){
  change <- 2
}else if(val >=(3 * ncombos) + 1 & val <= 4 * ncombos){
  change <- 8
}else if(val >=(4 * ncombos) + 1 & val <= 5 * ncombos){
  change <- 1.5
}else if(val >=(5 * ncombos) + 1 & val <=(6 * ncombos)){
  change <- 16
}

onCluster <- TRUE #Change to True if running on the cluster, False if running on macbook

if(onCluster==TRUE){
  def_wd <- "~/res/GEUV1Data/"
  def_wd2 <- "/pine/scr/s/k/skvanbur/GEUV1/"
  setwd(def_wd)
  func_loc <-  "~/code/CompFunctions.R"
  source(func_loc)
  #RATs_InfRep_save_dir <- paste0(def_wd, "RATsInfReps/")
}

ninfreps <- 100
#ninfreps <- 500

if(infReps=="Boot" & DRIMSeqFiltering==TRUE & TwentySamplesTotalAnalysis==FALSE){
  sub_dir1 <- "SalmonBootSamps/"
  RATs_InfRep_save_dir <- paste0(def_wd2, "RATsInfRepsBootDRIMSeqFiltering/")
  save_direc <- paste0(def_wd2, "RATsPowerResBootDRIMSeqFiltering/", "Change", change, "/")
}else if(infReps=="Boot" & DRIMSeqFiltering==TRUE & TwentySamplesTotalAnalysis==TRUE){
  sub_dir1 <- "SalmonBootSamps/"
  RATs_InfRep_save_dir <- paste0(def_wd2, "RATsInfRepsBootDRIMSeqFiltering/")
  save_direc <- paste0(def_wd2, "RATsPowerResBootDRIMSeqFilteringTwentySamples/", "Change", change, "/")
}

if(DRIMSeqFiltering==TRUE){
  #cntGenecntsScaledTPMFiltering is only used to get the transcripts/genes that pass filtering
  load("cntGenecntsScaledTPMFiltered.RData")
  load("cntGene.RData")
  load("cntDatasetscntsScaledTPMNoOtherGroupsFiltered.RData")
  load("tx2gene.RData")
  
  cntDatasets <- cntDatasetsFiltered
}else if(DRIMSeqFiltering==FALSE){
  load("cntGene.RData")
  load("cntDatasets.RData")
  load("tx2gene.RData")
}


if(TwentySamplesTotalAnalysis==TRUE){
  rm(Group)
  load("grpcombosAndsub_keyDRIMSeqFilteringPowerAnalysis1TwentySamples.RData")
  grpcombos <- grpcombos_TwentySamples
}else{
  rm(Group)
  load("FilesForPowerAnalysis1.RData")
}





Group <- as.factor(as.character(grpcombos[GroupNum,]))

curr_group <- Group

curr_key <- data.frame(sub_key$Sample, curr_group, sub_key$Identifier, stringsAsFactors = F)
colnames(curr_key) <- c("Sample", "Condition", "Identifier")

curr_key$Condition <- relevel(as.factor(curr_key$Condition), ref = 1)

curr_key_to_use <- curr_key

# if(TwentySamplesTotalAnalysis==TRUE){
#   #Take the first 10 samples in each condition from the 100 samples used from power analysis 1 based on the true condition assignment and only use data from those 
#   #s_c1 and s_c2 will be the sample numbers for the first and second conditions respectively
#   rm(Group)
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
#   #Need to generate a new group combination because simply subsetting from the existing one would not ensure each condition had 10 samples exactly
#   set.seed(53 * GroupNum)
#   
#   pos <- sample(1:20, size = 10, replace = FALSE)
#   curr_cond[pos] <- levels(curr_cond)[1]
#   curr_cond[-pos] <- levels(curr_cond)[2]
#   
#   curr_key_to_use <- subset(curr_key, curr_key$Identifier %in% samps)
#   curr_key_to_use$Condition <- curr_cond
#   Group <- curr_key_to_use$Condition
#   print(Group)
# }else{
#   
#   curr_key_to_use <- curr_key
# }


genestochange <- loadRData("geneschangedDRIMSeqFilteringPowerAnalysis1.RData")



cntGeneTemp <- subset(cntGene, (cntGene$MajorTrans==1 & cntGene$gene_id %in% genestochange))
transtochange <- cntGeneTemp$tx_id


sampsCondA <- subset(curr_key_to_use$Identifier, curr_key_to_use$Condition==levels(curr_key_to_use$Condition)[1])
sampsCondB <- subset(curr_key_to_use$Identifier, curr_key_to_use$Condition==levels(curr_key_to_use$Condition)[2])

InfRepsListCondA <- vector(mode = "list", length = length(sampsCondA))
InfRepsListCondB <- vector(mode = "list", length = length(sampsCondB))
print(paste0("Number of samples in Condition A is ", length(sampsCondA)))
print(paste0("Number of samples in Condition B is ", length(sampsCondB)))

#Load the infReps for the samples for condition "A" and update them by the change value for the MajorTrans of those
#genes that are changed
stt <- proc.time()
for(k in 1:length(sampsCondA)){
  curr_samp <- sampsCondA[k]
  print(paste0("Currently loading infRep data for Sample ", curr_samp))
  curr_dat <- loadRData(paste0(RATs_InfRep_save_dir, "RATsInfRepData", curr_samp, ".RData"))
  
  print(paste0("Nrows of the current infRep file is ", nrow(curr_dat)))
  colstochange <- colnames(curr_dat)[colnames(curr_dat)!="tx_id" & colnames(curr_dat)!="feature_id"]
  if("tx_id" %in% colnames(curr_dat)){
    rowstochange <- rownames(curr_dat)[curr_dat$tx_id %in% transtochange]
  }else if("feature_id" %in% colnames(curr_dat)){
    rowstochange <- rownames(curr_dat)[curr_dat$feature_id %in% transtochange]
  }else if("target_id" %in% colnames(curr_dat)){
    rowstochange <- rownames(curr_dat)[curr_dat$target_id %in% transtochange]
  }else{
    stop("colnames of curr_dat of the transcript id needs to be tx_id, feature_id, or target_id")
  }
  
  if(length(colstochange)!=ninfreps){stop("Something goes wrong when trying to update infRep data for power analysis")}
  curr_dat_changed <- data.frame(curr_dat)
  colnames(curr_dat_changed) <- colnames(curr_dat)
  curr_dat_changed[rowstochange,colstochange] <- curr_dat_changed[rowstochange,colstochange] * change
  
  curr_len_col <- paste0(curr_samp, "Len")
  curr_len1 <- cntGene[,c("tx_id", curr_len_col), drop = FALSE]
  curr_len2 <- curr_len1[order(curr_len1$tx_id),]
  
  if(sum(curr_dat_changed$feature_id!=curr_len2$tx_id) !=0){stop()}
  
  curr_dat_changed2 <- curr_dat_changed
  for(o in 2:ncol(curr_dat_changed2)){
    if(is.character(curr_dat_changed2[,1])==FALSE){stop("First column needs to have the transcript names")}
    curr_dat_changed2[,o] <- calcScaledCountsRATsPower(curr_col_curr_dat_changed = curr_dat_changed[,o], curr_len2 = curr_len2, curr_samp = curr_samp)
  }
  
  
  colnames(curr_dat_changed2)[colnames(curr_dat_changed2)=="tx_id"] <- "target_id"
  colnames(curr_dat_changed2)[colnames(curr_dat_changed2)=="feature_id"] <- "target_id"
  
  #Restrict list of transcripts to those that pass filtering
  curr_dat_changed3 <- subset(curr_dat_changed2, curr_dat_changed2$target_id %in% cntGeneFiltered$tx_id)
  
  InfRepsListCondA[[k]] <- data.table(curr_dat_changed3)
  
  rm(curr_dat)
  rm(curr_dat_changed)
  rm(curr_dat_changed2)
  rm(curr_dat_changed3)
  #print(gc())
  gc()
}

names(InfRepsListCondA) <- sampsCondA
ct <- proc.time() - stt
print(ct)

#Load the infReps for the samples for condition "B" but these are not changed since only vals for condition "A" are to be changed
#for the power analysis
stt2 <- proc.time()

for(k in 1:length(sampsCondB)){
  curr_samp <- sampsCondB[k]
  curr_dat <- data.table(loadRData(paste0(RATs_InfRep_save_dir, "RATsInfRepData", curr_samp, ".RData")))
  colnames(curr_dat)[colnames(curr_dat)=="tx_id"] <- "target_id"
  colnames(curr_dat)[colnames(curr_dat)=="feature_id"] <- "target_id"
  
  curr_dat_2 <- subset(curr_dat, curr_dat$target_id %in% cntGeneFiltered$tx_id)
  InfRepsListCondB[[k]] <- curr_dat_2
  rm(curr_dat)
  gc()
}

names(InfRepsListCondB) <- sampsCondB
ct2 <- proc.time() - stt2
print(ct2)


#cntGeneSub <- subset(cntGene, cntGene$gene_id %in% names(cntDatasets))
#cntGeneSubUpdated <- cntGeneSub

cntGeneUpdated <- cntGene
cntGeneUpdatedColsToChange <- paste0(sampsCondA, "Cnt")
cntGeneUpdatedRowsToChange <- rownames(cntGeneUpdated)[cntGeneUpdated$tx_id %in% transtochange]

cntGeneUpdated[cntGeneUpdatedRowsToChange, cntGeneUpdatedColsToChange] <- cntGeneUpdated[cntGeneUpdatedRowsToChange,cntGeneUpdatedColsToChange] * change


#Nsamp is loaded in and will always be equal to 462 for the GEUV1 data even though the power analysis is run on fewer samples
  #This is fine because the code will subset to only take columns that are needed for the power analysis
  #And any updating of counts and calculation os ScaledTPM counts is done sample by sample such that doing extra for now will not be a problem
cntGeneUpdated2 <- cntGeneUpdated[,c("gene_id", "tx_id", paste0("Sample", 1:nsamp, "Cnt"), paste0("Sample", 1:nsamp, "Len"))]

cntGeneUpdatedcntsScaledTPM <- cntsToScaledTPMCounts(cnts = cntGeneUpdated2, nsamp = nsamp)
cntGeneUpdatedcntsScaledTPM2 <- merge(cntGeneUpdatedcntsScaledTPM, tx2gene, by = "tx_id")
cntGeneUpdatedcntsScaledTPM3 <- subset(cntGeneUpdatedcntsScaledTPM2, cntGeneUpdatedcntsScaledTPM2$tx_id %in% cntGeneFiltered$tx_id)
rm(cntGeneUpdated)
rm(cntGeneUpdated2)
rm(cntGeneUpdatedcntsScaledTPM)
rm(cntGeneUpdatedcntsScaledTPM2)
gc()
St2 <- proc.time()
res <- RATsObsAnalysis(cntGene = cntGeneUpdatedcntsScaledTPM3, cntDatasets = cntDatasets, key = curr_key_to_use, infRepDataA = InfRepsListCondA, 
                       infRepDataB = InfRepsListCondB, DRIMSeqFiltering = DRIMSeqFiltering)

CompTime <- proc.time() - St2

assign(paste0("RATsPowerResGroupCombo", GroupNum), res)
assign(paste0("RATsPowerResGroupCombo", GroupNum, "CompTime"), CompTime)


if(!dir.exists(save_direc)){dir.create(save_direc, recursive = TRUE)}
#geneschanged <- genestochange
save(list = c(paste0("RATsPowerResGroupCombo", GroupNum), paste0("RATsPowerResGroupCombo", GroupNum, "CompTime")),
     file = paste0(save_direc, "RATsPowerResGroupCombo", GroupNum, ".RData"))

print(gc())