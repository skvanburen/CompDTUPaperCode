
.libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
library(DRIMSeq)
library(data.table)
library(pryr)

func_loc <- "~/code/CompFunctions.R"
source(func_loc)

ninfreps <- 100

DRIMSeqFiltering <- TRUE

def_wd <- "/pine/scr/s/k/skvanbur/TCGABRCAAnalysis/"

npartstouse <- 100

setwd(def_wd)
  
load("abDatasetsNoOtherGroupsFiltered.RData")
abDatasets <- abDatasetsFiltered
FilteringMethod <- "DRIMSeq"

genestouse <- names(abDatasets)


#Set countsFromAbundance parameter, which controls how the counts are estimated 
#(see tximport help documentation for more information on the differences this parameter makes)
countsFromAbundance <- "scaledTPM"

array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if(array_val==1){
  BasalVsNonBasalComp <- TRUE
}else if(array_val==2){
  BasalVsNonBasalComp <- FALSE
}


if(BasalVsNonBasalComp==TRUE){
  custom_key <- key
  #"A" is for Basal and "B" is for non-Basal
  custom_key$Condition <- relevel(as.factor(ifelse(custom_key$Condition=="Basal", "A", "B")), ref = "A")
  save_dir <- "~/res/TCGABRCAAnalysis/DRIMSeqResultsBasalVsNonBasal/"
  if(!dir.exists(save_dir)){dir.create(save_dir)}
}else{
  custom_key <- key
  save_dir <- "~/res/TCGABRCAAnalysis/DRIMSeqResults/"
  if(!dir.exists(save_dir)){dir.create(save_dir)}
}

print("Table of condition values is given below")
print(table(custom_key$Condition))

sampstouse <- custom_key$Identifier


  #Only dmPrecision, dmFit, and dmTest statements are included in the computation time that is used to construct the table
    #Not any time to create/load in or reformat data, etc to keep this comparison fair between methods
  
  DRIMSeqObservedAnalysis(countsFromAbundance = countsFromAbundance, FilteringMethod = FilteringMethod, sampstouse = sampstouse,
                          save_dir = save_dir, genestouse = genestouse, fullReturn = TRUE, custom_key = custom_key)
  
  