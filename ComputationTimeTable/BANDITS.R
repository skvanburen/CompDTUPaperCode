#Code to calcuate computation time results for BANDITS


#Set onCluster to TRUE to run for all 462 samples (just to get an idea of memory requirements/time required)
  #Results from the paper were run on macbook to ensure more consistency between computation time comparison for different methods
onCluster <- FALSE
#Only consider the Twenty Samples analysis for BANDITS since the computation time is so much higher than DRIMSeq even
TwentySamplesTotalAnalysis <- FALSE

BANDITSDataAlreadyCreated <- TRUE

#To get the Salmon equivalence class files for the necessary samples, sun the following code on the cluster
  #to copy the files to another directory for easier download and to have them all in one place
  #For now only consider the 20 sample analysis
  #This list comes from the sampstouse list below
  # samps_to_keep <- c(4,7,8,12,19,20,21,37,39,41,45,47,50,51,53,56,64,66,79,80)
  # base_dir <- "/pine/scr/s/k/skvanbur/GEUV1/SalmonBootSamps/"
  # #Load these counts to get the key
  # load("~/res/GEUV1Data/cntGenecntsScaledTPMFiltered.RData")
  # key_to_use <- subset(key, key$Identifier %in% paste0("Sample", samps_to_keep))
  # dir_modifiers <- key_to_use$ENARun
  # 
  # for(di in 1:length(samps_to_keep)){
  #   print(paste0("current number is ", di))
  #   curr_eq_class_fil <- paste0(base_dir, dir_modifiers[di], "/", "aux_info/", "eq_classes.txt")
  #   save_dir <- paste0("/pine/scr/s/k/skvanbur/GEUV1/EqClassFilesForBANDITS/", "Sample", samps_to_keep[di], "/")
  #   if(!dir.exists(save_dir)){dir.create(save_dir, recursive = TRUE)}
  #   system(paste0("cp ", curr_eq_class_fil, " ", paste0(save_dir, "EqClassFile.txt")))
  # }

  #Now, create a directory that has all samples
  # base_dir <- "/pine/scr/s/k/skvanbur/GEUV1/SalmonBootSamps/"
  # #Load these counts to get the key
  # load("~/res/GEUV1Data/cntGenecntsScaledTPMFiltered.RData")
  # samps_to_keep <- 1:462
  # key_to_use <- key
  # dir_modifiers <- key_to_use$ENARun
  # 
  # for(di in 1:length(samps_to_keep)){
  #   print(paste0("current number is ", di))
  #   curr_eq_class_fil <- paste0(base_dir, dir_modifiers[di], "/", "aux_info/", "eq_classes.txt")
  #   save_dir <- paste0("/pine/scr/s/k/skvanbur/GEUV1/EqClassFilesForBANDITSAllSamps/", "Sample", samps_to_keep[di], "/")
  #   if(!dir.exists(save_dir)){dir.create(save_dir, recursive = TRUE)}
  #   system(paste0("cp ", curr_eq_class_fil, " ", paste0(save_dir, "EqClassFile.txt")))
  # }

if(onCluster==TRUE){
  .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
  func_loc <- "~/code/CompFunctions.R"
  
  #Only consider the Twenty Samples analysis for BANDITS since the computation time is so much higher than DRIMSeq even
  TwentySamplesTotalAnalysis <- FALSE
  
  load("~/res/GEUV1Data/tx2gene.RData")
  
  def_wd <- "~/res/GEUV1Data/"
  
  EQFilesDir <- "/pine/scr/s/k/skvanbur/GEUV1/EqClassFilesForBANDITSAllSamps/"
  
  BANDITS_Data_File <- paste0("/pine/scr/s/k/skvanbur/GEUV1/", "BANDITSInputData462Samples.RData")
  
  if(!dir.exists("/pine/scr/s/k/skvanbur/GEUV1/")){dir.create("/pine/scr/s/k/skvanbur/GEUV1/", recursive = TRUE)}
}else{
  func_loc <- "/Users/Scott/Documents/Dissertation/Code/CompFunctions.R"

  
  def_wd <- "/Users/Scott/Documents/Dissertation Data/GEUV1Data/"
  
  load("/Users/Scott/Documents/Dissertation Data/GEUV1Data/tx2gene.RData")
  
  EQFilesDir <- "/Users/Scott/Documents/Dissertation Data/GEUV1Data/EqClassFilesForBANDITS/"
  
  if(TwentySamplesTotalAnalysis==TRUE){
    BANDITS_Data_File <- paste0("/Users/Scott/Documents/Dissertation Data/GEUV1Data/", "BANDITSInputData20Samples.RData")
  }else{
    BANDITS_Data_File <- paste0("/Users/Scott/Documents/Dissertation Data/GEUV1Data/", "BANDITSInputData462Samples.RData")
  }

  if(!dir.exists("/Users/Scott/Documents/Dissertation Data/GEUV1Data/")){dir.create("/Users/Scott/Documents/Dissertation Data/GEUV1Data/", recursive = TRUE)}
}

library(BANDITS)

setwd(def_wd)

source(func_loc)

ninfreps <- 100

dataset <- "GEUV1"

if(onCluster==TRUE){
  save_dir <- "~/res/GEUV1Data/BANDITSFullDataComputationTime/"
}else{
  if(TwentySamplesTotalAnalysis==TRUE){
    save_dir <- "/Users/Scott/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/GEUV1TwentySamples/"
  }else{
    save_dir <- "/Users/Scott/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/GEUV1/"
  }
}

  
if(!dir.exists(save_dir)){dir.create(save_dir, recursive = TRUE)}
  #These counts (which have been scaled to the TPM values) are only loaded to get the list of transcripts
    #that pass filtering (based on DRIMSeq's filtering) and the effective lengths of the transcripts for each sample
    #This will ensure the transcripts used matche the other methods
    #These counts are not actually used by BANDITS, but the effective lengths from this file are needed
  load("cntGenecntsScaledTPMFiltered.RData")
  tx_to_use <- cntGeneFiltered$tx_id


genestouse <- unique(cntGeneFiltered$gene_id)


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

transcript_counts <- cntGeneFiltered[tx_to_use, paste0(sampstouse, "Cnt")]


tx2gene2 <- tx2gene[,-3]
gene_tr_id <- subset(tx2gene2, tx2gene2$tx_id %in% tx_to_use)

key_to_use <- subset(key, key$Identifier %in% sampstouse)
samples_design <- key_to_use[,c("Identifier", "Condition")]
colnames(samples_design)[colnames(samples_design)=="Identifier"] <- "sample_id"
colnames(samples_design)[colnames(samples_design)=="Condition"] <- "group"
levels(samples_design$group) <- relevel(samples_design$group, ref = "CEU")

#Ensure the order of files in the EQClassFiles list is the same as the samples in samples_design
EQClassFilesT <- list.files(EQFilesDir, recursive = TRUE, full.names = TRUE)
EQClassFiles <- gtools::mixedsort(EQClassFilesT)

print(paste0("The number of equivalence class files found is ", length(EQClassFiles)))

#Compute the median effective lengths of the transcripts
length_vals <- as.matrix(cntGeneFiltered[,paste0(sampstouse, "Len")])
medians_length_vals <- matrixStats::rowMedians(length_vals)
names(medians_length_vals) <- rownames(length_vals)
medians_length_vals_to_use <- medians_length_vals[tx_to_use]

#Create BANDITS data object
if(BANDITSDataAlreadyCreated==TRUE){
  load(BANDITS_Data_File)
}else{
  CDT1 <- proc.time()
  input_data <- create_data(salmon_or_kallisto = "salmon",
                            gene_to_transcript = gene_tr_id,
                            salmon_path_to_eq_classes = EQClassFiles,
                            eff_len = medians_length_vals_to_use, 
                            n_cores = 1,
                            transcripts_to_keep = tx_to_use)
  TimeToCreateBANDITSData <- proc.time() - CDT1
  
  ObjSizeInput_Data <- object.size(input_data)
  
  save(input_data, TimeToCreateBANDITSData,  ObjSizeInput_Data, file = BANDITS_Data_File)
}




#Run the (highly recommended) step to calculate an informative prior for the precision parameter
CDT2 <- proc.time()
BANDITSprecision <-  prior_precision(gene_to_transcript = gene_tr_id,
                            transcript_counts = transcript_counts, n_cores = 1,
                            transcripts_to_keep = tx_to_use)
TimeToCreateBANDITSPrecision <- proc.time() - CDT2
print(paste0("Time to Create BANDITS Precision is ", TimeToCreateBANDITSPrecision[3]))
save(BANDITSprecision, TimeToCreateBANDITSPrecision, file = "/Users/Scott/BANDITSPrecision.RData")


#Run the BANDITS DTU step
CDT3 <- proc.time()
set.seed(87203)
results <-  test_DTU(BANDITS_data = input_data,
                   precision = BANDITSprecision$prior,
                   samples_design = samples_design,
                   group_col_name = "group", n_cores = 1,
                   gene_to_transcript = gene_tr_id)
TimeToRunBANDITS <- proc.time() - CDT3
print(paste0("Time to run BANDITS is ", TimeToRunBANDITS[3]))


GeneLevelSignRes <- top_genes(results)

if(onCluster==TRUE){
  save_fil <- paste0(save_dir, "BANDITSComputationTimeResFull.RData")
  print(gc())
  save(TimeToCreateBANDITSData, TimeToCreateBANDITSPrecision, TimeToRunBANDITS, ObjSizeInput_Data, GeneLevelSignRes, file = save_fil)
}else{
  if(TwentySamplesTotalAnalysis==TRUE){
    save_fil <- "/Users/Scott/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/GEUV1TwentySamples/BANDITSResDRIMSeqFiltering.RData"
  }else{
    save_fil <- "/Users/Scott/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/GEUV1/BANDITSResDRIMSeqFiltering.RData"
  }
  
  if(BANDITSDataAlreadyCreated==TRUE){
    save(TimeToCreateBANDITSPrecision, TimeToRunBANDITS, GeneLevelSignRes, results, file = save_fil)
  }else{
    save(TimeToCreateBANDITSData, TimeToCreateBANDITSPrecision, TimeToRunBANDITS, GeneLevelSignRes, results, file = save_fil)
  }

}





