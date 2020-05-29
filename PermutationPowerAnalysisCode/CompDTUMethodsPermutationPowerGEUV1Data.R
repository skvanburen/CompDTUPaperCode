#This file calculates the permutation power analysis results for the CompDTU methods

onCluster <- TRUE



#Calculate results for data where abundance values come from the mean/median of bootstrap samples
  #If this option is true only CompDTU calculated on mean/median vals of boot/gibbs samples will be run 
  #(Along with CompDTU on regular point estimates just for comparison)
CalcCompDTUResForAbFromInfReps <- FALSE


#Set to TRUE to run on the subset of 20 samples, FALSE to run on 100 sample analysis
TwentySamplesTotalAnalysis <- TRUE

BootSampsDRIMSeqFiltering <- TRUE
DRIMSeqFiltering <- TRUE
BootSamps <- TRUE

#Array value controls the group combination, part number and change value
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  
  if(array_val < 1000){
    #Sys.sleep(array_val)
    Sys.sleep(sample(1:100,1))
  }else{
    Sys.sleep(sample(1:100,1))
  }

  
  #Add the proper libpaths to R so it knows where to load packages from
  if(version$nickname=="Planting of a Tree"){
    .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
  }


library(compositions)
library(gtools)
library(data.table)



#Setting curr_change to 1 gives results corresponding to the real observed data, while change 4 represents modified data
if(array_val <=10000){
  curr_change <- 4
  
  groupcombo <- ceiling(array_val/100)
  
  
  PartNum <- array_val %% 100
  if(PartNum==0){
    PartNum <- 100
  }
  
}else if(array_val >10000 & array_val <=20000){
  curr_change <- 2
  
  t_val <- array_val - 10000
  
  groupcombo <- ceiling(t_val/100)
  
  
  PartNum <- t_val %% 100
  if(PartNum==0){
    PartNum <- 100
  }
  
}else if(array_val >20000 & array_val <=30000){
  curr_change <- 1
  
  t_val <- array_val - 20000
  
  groupcombo <- ceiling(t_val/100)
  
  
  PartNum <- t_val %% 100
  if(PartNum==0){
    PartNum <- 100
  }
  
}else if(array_val >30000 & array_val <=40000){
  curr_change <- 8
  
  t_val <- array_val - 30000
  
  groupcombo <- ceiling(t_val/100)
  
  
  PartNum <- t_val %% 100
  if(PartNum==0){
    PartNum <- 100
  }
  
}


CalculateCompDTUmeRes <- TRUE
  
  if(onCluster==TRUE){
    source("/nas/longleaf/home/skvanbur/code/CompFunctions.R")
    source("/nas/longleaf/home/skvanbur/code/CompDTUMethodsPermutationPowerFunctions.R")
    
    #For now load the Updated power results that exist to compare to
      #These are other groups and not DRIMSeqFiltered
    
    if(TwentySamplesTotalAnalysis==TRUE){
      dir_modifier <- "TwentySamples"
    }else{
      dir_modifier <- ""
    }
    
    if(BootSampsDRIMSeqFiltering==TRUE){
      BaseFilesDir <- paste0("/pine/scr/s/k/skvanbur/GEUV1/UpdatedPowerDataBootDRIMSeqFiltering", dir_modifier, "/Change", curr_change, "/GroupCombo", groupcombo, "/")
      BootDataT <- loadRData(paste0("/pine/scr/s/k/skvanbur/GEUV1/UpdatedPowerDataBootDRIMSeqFiltering", dir_modifier, "/Change", curr_change, "/GroupCombo", groupcombo, "/NewabDatasetsGibbsChange", curr_change, "GroupCombo", groupcombo, "Part", PartNum, ".RData"))
      NonBootDataT <- loadRData(paste0("/pine/scr/s/k/skvanbur/GEUV1/UpdatedPowerDataBootDRIMSeqFiltering", dir_modifier, "/Change", curr_change, "/GroupCombo", groupcombo, "/UpdatedabDatasetsChange", curr_change, "GroupCombo", groupcombo, ".RData"))
    }else if(DRIMSeqFiltering==TRUE & BootSamps==FALSE){
      BootDataT <- loadRData(paste0("/pine/scr/s/k/skvanbur/GEUV1/UpdatedPowerDataDRIMSeqFiltering", dir_modifier, "/Change", curr_change, "/GroupCombo", groupcombo, "/NewabDatasetsGibbsChange", curr_change, "GroupCombo", groupcombo, "Part", PartNum, ".RData"))
      NonBootDataT <- loadRData(paste0("/pine/scr/s/k/skvanbur/GEUV1/UpdatedPowerDataDRIMSeqFiltering", dir_modifier, "/Change", curr_change, "/GroupCombo", groupcombo, "/UpdatedabDatasetsChange", curr_change, "GroupCombo", groupcombo, ".RData"))
    }else{
      BootDataT <- loadRData(paste0("/pine/scr/s/k/skvanbur/GEUV1/UpdatedPowerData", dir_modifier, "/Change", curr_change, "/GroupCombo", groupcombo, "/NewabDatasetsGibbsChange", curr_change, "GroupCombo", groupcombo, "Part", PartNum, ".RData"))
      NonBootDataT <- loadRData(paste0("/pine/scr/s/k/skvanbur/GEUV1/UpdatedPowerData", dir_modifier, "/Change", curr_change, "/GroupCombo", groupcombo, "/UpdatedabDatasetsChange", curr_change, "GroupCombo", groupcombo, ".RData"))
    }
    
    if(CalcCompDTUResForAbFromInfReps==TRUE){
      wd2 <- "/pine/scr/s/k/skvanbur/GEUV1/"
      UpdatedAbDatasetsAbFromInfRepsDir <- paste0(wd2, "UpdatedAbDatasetsAbFromInfRepsDRIMSeqFiltering", dir_modifier, "/Change", curr_change, "/", "GroupCombo", groupcombo, "/")
      
      fil_mod <- returnFilMod(abFromInfRepFunc = "mean", GibbsSamps = FALSE)
      
      DatMeanBoot <- loadRData(paste0(UpdatedAbDatasetsAbFromInfRepsDir, "UpdatedabDatasetsChange", curr_change, "GroupCombo", groupcombo, fil_mod, ".RData"))
      
    }else{
      UpdatedAbDatasetsAbFromInfRepsDir <- NULL
      DatMeanBoot <- NULL
      
    }
    
    load("~/res/GEUV1Data/abGeneFiltered.RData")
    
    #Restrict to only genes that pass DRIMSeq filtering (For now the transcripts come from other groups, but at least don't include
    #tons of genes that don't pass these filtering criteria)
    genes_to_use <- unique(abGeneFiltered$gene_id)
    txsfiltered <- unique(abGeneFiltered$tx_id)
    BootData <- subset(BootDataT, names(BootDataT) %in% genes_to_use)
    NonBootData <- subset(NonBootDataT, names(NonBootDataT) %in% genes_to_use)
    rm(BootDataT)
    rm(NonBootDataT)
    gc()

    rm(nsamp)

    if(TwentySamplesTotalAnalysis==TRUE){
      load("/nas/longleaf/home/skvanbur/res/GEUV1Data/grpcombosAndsub_keyDRIMSeqFilteringPowerAnalysis1TwentySamples.RData")
      AllGroupCombinations <- grpcombos_TwentySamples
        
      samps <- sub_key$Identifier
      
      nsamp <- 20
      nboot <- 100
      ngrpcombos <- 100
    }else{
      load("/nas/longleaf/home/skvanbur/res/GEUV1Data/FilesForPowerAnalysis1.RData")
      AllGroupCombinations <- grpcombos
      
      samps <- sub_key$Identifier
      
      nsamp <- 100
      nboot <- 100
      ngrpcombos <- 100
    }
    
  }


curr_cond <- as.factor(AllGroupCombinations[groupcombo,])

if(length(curr_cond)!=nsamp | length(curr_cond[levels(curr_cond)[1]])!=length(curr_cond[levels(curr_cond)[2]])){
  stop("The condition variable is misspecified")
}
print(curr_cond)



st1 <- proc.time()
#[1:5] restricts the run to the first 5 genes to enable testing
#Takes about 200 seconds per gene if running with CompDTUme, more like 30-45 if not
  #And, for this dataset with this structure there are no more than 300 genes per part
genes <- names(BootData)

if(is.null(genes)){
  curr_res <- NULL
}else if(length(genes)==0){
  curr_res <- NULL
}else{
  curr_res <- CompDTUMethodsPowerAnalysis(genes = genes, curr_cond = curr_cond, onCluster = onCluster, 
                                        AllGroupCombinations = AllGroupCombinations, curr_change = curr_change, CalculateCompDTUmeRes = CalculateCompDTUmeRes, nsamp = nsamp,
                                        nboot = nboot, ngrpcombos = ngrpcombos, samps = samps, BootSamps = BootSamps, GEUV1Data = TRUE,
                                        txsfiltered = txsfiltered, CalcCompDTUResForAbFromInfReps = CalcCompDTUResForAbFromInfReps, 
                                        UpdatedAbDatasetsAbFromInfRepsDir = UpdatedAbDatasetsAbFromInfRepsDir,DatMeanBoot = DatMeanBoot)
}

comptime <- proc.time() - st1
print(comptime)

if(CalcCompDTUResForAbFromInfReps==TRUE){
  direc_modifier <- "CompDTUResForAbFromInfReps"
}else{
  direc_modifier <- ""
}

if(onCluster==TRUE){
  assign(paste0("PermuteImputeRes", PartNum), curr_res)
    if(CalcCompDTUResForAbFromInfReps==TRUE){
      if(TwentySamplesTotalAnalysis==TRUE){
        direc <- paste0("/pine/scr/s/k/skvanbur/GEUV1/CompDTUMethodsPowerAnalysisResCompDTUResForAbFromInfRepsTwentySamples/","ActualData", direc_modifier, "/", "Change", curr_change, "/", "GroupCombo", groupcombo, "/")
      }else{
        direc <- paste0("/pine/scr/s/k/skvanbur/GEUV1/CompDTUMethodsPowerAnalysisResCompDTUResForAbFromInfReps/","ActualData", direc_modifier, "/", "Change", curr_change, "/", "GroupCombo", groupcombo, "/")
      }
    }else if(BootSampsDRIMSeqFiltering==TRUE){
      if(TwentySamplesTotalAnalysis==TRUE){
        direc <- paste0("/pine/scr/s/k/skvanbur/GEUV1/CompDTUMethodsPowerAnalysisResBootDRIMSeqFilteringTwentySamples/","ActualData", direc_modifier, "/", "Change", curr_change, "/", "GroupCombo", groupcombo, "/")
      }else{
        direc <- paste0("/pine/scr/s/k/skvanbur/GEUV1/CompDTUMethodsPowerAnalysisResBootDRIMSeqFiltering/","ActualData", direc_modifier, "/", "Change", curr_change, "/", "GroupCombo", groupcombo, "/")
      }

    }else if(DRIMSeqFiltering==TRUE & BootSamps==FALSE){
      if(TwentySamplesTotalAnalysis==TRUE){
        direc <- paste0("/pine/scr/s/k/skvanbur/GEUV1/CompDTUMethodsPowerAnalysisResGibbsDRIMSeqFilteringTwentySamples/","ActualData", direc_modifier, "/", "Change", curr_change, "/", "GroupCombo", groupcombo, "/")
      }else{
        direc <- paste0("/pine/scr/s/k/skvanbur/GEUV1/CompDTUMethodsPowerAnalysisResGibbsDRIMSeqFiltering/","ActualData", direc_modifier, "/", "Change", curr_change, "/", "GroupCombo", groupcombo, "/")
      }

    }else{
      if(TwentySamplesTotalAnalysis==TRUE){
        direc <- paste0("/pine/scr/s/k/skvanbur/GEUV1/CompDTUMethodsPowerAnalysisResTwentySamples/","ActualData", direc_modifier, "/", "Change", curr_change, "/", "GroupCombo", groupcombo, "/")
      }else{
        direc <- paste0("/pine/scr/s/k/skvanbur/GEUV1/CompDTUMethodsPowerAnalysisRes/","ActualData", direc_modifier, "/", "Change", curr_change, "/", "GroupCombo", groupcombo, "/")
      }
    }
  
  
  if(!dir.exists(direc)) {dir.create(direc, recursive = TRUE)}
  print(paste0("Results saved to ", direc))
  save(list = paste0("PermuteImputeRes", PartNum), file = paste0(direc, "GroupCombo", groupcombo, "PermuteImputeRes", PartNum, ".RData"), version = 2)
  print(gc())
}

#This code can help you confirm every job ran properly (and find those that did not run properly- usually due to random slurm error)
#Then, just resubmit the ones that did not work
#In general comment this code out, but put this is.numeric check to make sure it doesn't get executed accidently if you forget to comment it out
if(!is.numeric(array_val)){

  # direc <- "/pine/scr/s/k/skvanbur/GEUV1/CompDTUMethodsPowerAnalysisRes/ActualData/"
  #direc <- "/pine/scr/s/k/skvanbur/GEUV1/CompDTUMethodsPowerAnalysisResBootDRIMSeqFiltering/ActualData/"
  #direc <- "/pine/scr/s/k/skvanbur/GEUV1/CompDTUMethodsPowerAnalysisResBootDRIMSeqFilteringTwentySamples/ActualData/"
  # direc <- "/pine/scr/s/k/skvanbur/GEUV1/CompDTUMethodsPowerAnalysisResGibbsDRIMSeqFiltering/ActualData/"
  # #direc <- "/pine/scr/s/k/skvanbur/GEUV1/CompDTUMethodsPowerAnalysisResGibbsDRIMSeqFilteringTwentySamples/ActualData/"
  #direc <- paste0("/pine/scr/s/k/skvanbur/GEUV1/CompDTUMethodsPowerAnalysisResCompDTUResForAbFromInfReps/","ActualDataCompDTUResForAbFromInfReps/")
  #direc <- paste0("/pine/scr/s/k/skvanbur/GEUV1/CompDTUMethodsPowerAnalysisResCompDTUResForAbFromInfRepsTwentySamples/","ActualDataCompDTUResForAbFromInfReps/")
  # #All

  #changes_to_use <- c(1,2,4) #Use this one for the 100 sample case
#changes_to_use <- c(1,2,4,8) #Use this one for the twenty sample case
changes_to_use <- c(1,2,4)
  vals_needed_to_resubmit <- c()
  for(chv in changes_to_use){
    change <- chv
    print(paste0("Current change is ", change))
    if(change==4){
      amt_to_add <- 0
    }else if(change==2){
      amt_to_add <- 10000
    }else if(change==1){
      amt_to_add <- 20000
    }else if(change==8){
      amt_to_add <- 30000
    }
    ncombos <- 100
    nparts <- 100
    #
    for(j in 1:ncombos){
      print(paste0("j is ", j))
      for(k in 1:nparts){


        if(!file.exists(paste0(direc, "Change", change, "/GroupCombo", j, "/", "GroupCombo", j, "PermuteImputeRes", k, ".RData"))){
          print(paste0("File does not exist for GroupCombo ", j, " and Part ", k))

          vals_needed_to_resubmit <- c(vals_needed_to_resubmit, (j-1)*100 + k + amt_to_add)
        }else if(file.size(paste0(direc, "Change", change, "/GroupCombo", j, "/", "GroupCombo", j, "PermuteImputeRes", k, ".RData"))==0){
          print(paste0("Job did not work properly for for GroupCombo ", j, " and Part ", k))

          vals_needed_to_resubmit <- c(vals_needed_to_resubmit, (j-1)*100 + k + amt_to_add)
        }
      }
    }

    # print(vals_needed_to_resubmit)
    # for(g in 1:(length(vals_needed_to_resubmit))){
    #   curr <- vals_needed_to_resubmit[g]
    #   system(paste0("module load r/3.6.0"))
    #   system(paste0("sbatch --array=", curr, " ~/res/GEUV1Data/ExamineImputeandPermuteApproachGEUV1Data.sh"))
    # }

  }

}