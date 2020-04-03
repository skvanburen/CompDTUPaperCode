#Calculate InfRV value for each gene as described in Nonparametric expression analysis using inferential replicate counts (Anqi Zhu et al 2019)

#Directory where the file "SaveGibbsDataAsRData" saves files of interest
#For now use the bootstrap samples

array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

dir1 <- "/pine/scr/s/k/skvanbur/GEUV1/"

if(array_val==1){
  GibbsSamps <- TRUE
  sub_dir1 <- "Salmon/GibbsSamps/"
}else if(array_val==2){
  GibbsSamps <- FALSE
  sub_dir1 <- "SalmonBootSamps/BootSamps/"
}

if(version$nickname=="Planting of a Tree"){
  .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
}else{
  source("~/.RProfile")
}

library(compositions)
library(reshape2)
source("~/code/CompFunctions.R")

dir2 <- "~/res/GEUV1Data/"

func_loc <-  "~/code/CompFunctions.R"
def_wd1 <- paste0(dir1, sub_dir1)
setwd(def_wd1)

load(paste0(dir2,"tx2gene.RData"))
load(paste0(dir2,"abGeneFiltered.RData"))

load("/nas/longleaf/home/skvanbur/res/GEUV1Data/FilesForPowerAnalysis1.RData")
samps <- sub_key$Identifier
nsamp <- length(samps)

source(func_loc)


#Set type = to either "Cnt" or "TPM"

InfRVRes <- calcInfRVRes(type = "Cnt", nsamp = nsamp, samps = samps, abGeneFiltered = abGeneFiltered, tx2gene = tx2gene, GibbsSamps = GibbsSamps, def_wd1 = def_wd1)
InfRVResTPM <- calcInfRVRes(type = "TPM", nsamp = nsamp, samps = samps, abGeneFiltered = abGeneFiltered, tx2gene = tx2gene, GibbsSamps = GibbsSamps, def_wd1 = def_wd1)

st1 <- proc.time()
InfRepVarilr <- calcInfRVRes(type = "ilr", nsamp = nsamp, samps = samps, abGeneFiltered = abGeneFiltered, tx2gene = tx2gene, GibbsSamps = GibbsSamps, def_wd1 = def_wd1)
ct <- proc.time() - st1
print(ct)

colnames(InfRVRes) <- c("gene_id", "MeanGeneMaxInfRVCnt", "MeanGeneMedianInfRVCnt", "MeanGeneMeanInfRVCnt")
colnames(InfRVResTPM) <- c("gene_id", "MeanGeneMaxInfRVTPM", "MeanGeneMedianInfRVTPM", "MeanGeneMeanInfRVTPM")

InfRVResAllT <- merge(InfRVRes, InfRVResTPM, by = "gene_id")
InfRVResAll <- merge(InfRVResAllT, InfRepVarilr, by = "gene_id")


if(GibbsSamps==TRUE){
  InfRVResAllGibbsSamps <- InfRVResAll
  save(InfRVResAllGibbsSamps, file = paste0(dir2, "InfRVResAllForGibbsSamps.RData"), version = 2)
}else{
  InfRVResAllBootSamps <- InfRVResAll
  save(InfRVResAllBootSamps, file = paste0(dir2, "InfRVResAllForBootSamps.RData"), version = 2)
}


#save(InfRVRes, file = paste0(dir2, "InfRVRes.RData"), version = 2)
#save(InfRVResTPM, file = paste0(dir2, "InfRVResTPM.RData"), version = 2)