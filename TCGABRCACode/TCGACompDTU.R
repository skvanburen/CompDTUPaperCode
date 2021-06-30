#Sample code to run the compositional regression models for DTU analysis

library(CompDTUReg)
library(pryr)

#Specify the same directory where the gene level files to be used for the analysis were saved in file (3)
def_wd <- "/pine/scr/s/k/skvanbur/TCGABRCAAnalysis/"

#Directory where the gene level files were saved in file (3)SaveNecessaryDatasetsForCompDTUReg.R  (or whatever directory the GeneLevelFiles to be used in the analysis are saved in)
GeneLevelFilesSaveDir <- paste0(def_wd, "GeneLevelFiles/")

#Generate list of all gene level files
GeneFiles <- list.files(GeneLevelFilesSaveDir, full.names = TRUE)

save_dir <- "~/res/TCGABRCAAnalysis/CompDTURes/"
if(!dir.exists(save_dir)){dir.create(save_dir)}

setwd(def_wd)
#Now, examples providing custom specified design matrices under the null and alternative hypotheses

array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if(array_val==1){
  BasalVsNonBasalComp <- TRUE
}else if(array_val==2){
  BasalVsNonBasalComp <- FALSE
}

#Load the first gene-level file to get the group (cond) information
load("key.RData")
if(BasalVsNonBasalComp==TRUE){
  file_to_save <- paste0(save_dir, "CompDTUObsResBasalVsNonBasal.RData")
  custom_key <- key
  #"A" is for Basal and "B" is for non-Basal
  custom_key$Condition <- relevel(as.factor(ifelse(custom_key$Condition=="Basal", "A", "B")), ref = "A")
  cond <- custom_key$Condition
  
  key_temp <- CompDTUReg::loadRData(GeneFiles[1], objNameToGet = "key")
  if(sum(custom_key$Identifier!=key_temp$Identifier)!=0){
    stop("There is something wrong with the ordering of the condition variable")
  }
}else{
  file_to_save <- paste0(save_dir, "CompDTUObsRes.RData")
  cond <- key$Condition
  key_temp <- CompDTUReg::loadRData(GeneFiles[1], objNameToGet = "key")
  if(sum(key$Identifier!=key_temp$Identifier)!=0){
    stop("There is something wrong with the ordering of the condition variable")
  }
}

print("The condition variable used is printed below")
print(cond)


#Specify null and alternative design matrices.  The rows must be in the same order as key$Identifier, where key is extracted above
  #These results are for comparison of all five types
NullDesign <- model.matrix(~1, data = cond)
AltDesign <- model.matrix(~cond)

#These results test for DTU (ie the significance of the Group (cond) variable) controlling for pred1 and pred2
st1 <- proc.time()
mem_change_CompDTU <- mem_change(CompDTUResults <- rbindlist(lapply(GeneFiles, startCompDTUReg, runWithME = FALSE, customHypTest = TRUE, NullDesign = NullDesign, AltDesign = AltDesign)))
CompDTUTime <- proc.time() - st1

st2 <- proc.time()
mem_change_CompDTUme <- mem_change(CompDTUmeResults <- rbindlist(lapply(GeneFiles, startCompDTUReg, runWithME = TRUE, customHypTest = TRUE, NullDesign = NullDesign, AltDesign = AltDesign)))
CompDTUmeTime <- proc.time() - st2

save(mem_change_CompDTU, mem_change_CompDTUme, CompDTUResults, CompDTUmeResults, CompDTUTime, CompDTUmeTime, file = file_to_save)
print(gc())