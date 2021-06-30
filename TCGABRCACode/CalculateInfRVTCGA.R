#Calculate InfRV value for each gene as described in Nonparametric expression analysis using inferential replicate counts (Anqi Zhu et al 2019)

#Directory where the file "SaveGibbsDataAsRData" saves files of interest
#For now use the bootstrap samples

dir1 <- "/pine/scr/s/k/skvanbur/TCGABRCAAnalysis/"

GibbsSamps <- FALSE
sub_dir1 <- "SalmonBootSamps/BootSamps/"

.libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")

library(compositions)
library(reshape2)
source("~/code/CompFunctions.R")

dir2 <- "~/res/TCGABRCAAnalysis/"

func_loc <-  "~/code/CompFunctions.R"
def_wd1 <- paste0(dir1, sub_dir1)
setwd(def_wd1)

load(paste0(dir1,"tx2gene.RData"))
load(paste0(dir1,"abGeneFiltered.RData"))

load(paste0(dir1, "key.RData"))
samps <- key$Identifier
nsamp <- length(samps)

source(func_loc)


#Set type = to either "Cnt" or "TPM"

InfRVRes <- calcInfRVRes(type = "Cnt", nsamp = nsamp, samps = samps, abGeneFiltered = abGeneFiltered, tx2gene = tx2gene, GibbsSamps = GibbsSamps, def_wd1 = def_wd1)
InfRVResTPM <- calcInfRVRes(type = "TPM", nsamp = nsamp, samps = samps, abGeneFiltered = abGeneFiltered, tx2gene = tx2gene, GibbsSamps = GibbsSamps, def_wd1 = def_wd1)


colnames(InfRVRes) <- c("gene_id", "MeanGeneMaxInfRVCnt", "MeanGeneMedianInfRVCnt", "MeanGeneMeanInfRVCnt")
colnames(InfRVResTPM) <- c("gene_id", "MeanGeneMaxInfRVTPM", "MeanGeneMedianInfRVTPM", "MeanGeneMeanInfRVTPM")

InfRVResAll <- merge(InfRVRes, InfRVResTPM, by = "gene_id")


if(GibbsSamps==TRUE){
  InfRVResAllGibbsSamps <- InfRVResAll
  save(InfRVResAllGibbsSamps, file = paste0(dir2, "InfRVResAllForGibbsSamps.RData"))
}else{
  InfRVResAllBootSamps <- InfRVResAll
  save(InfRVResAllBootSamps, file = paste0(dir2, "InfRVResAllForBootSamps.RData"))
}

print(gc())
#save(InfRVRes, file = paste0(dir2, "InfRVRes.RData"), version = 2)
#save(InfRVResTPM, file = paste0(dir2, "InfRVResTPM.RData"), version = 2)