#Generate ROC Curves Using results generated in SummarizePowerResults
library(data.table)
library(tikzDevice)

onCluster <- FALSE #Change to True if running on the cluster, False if running on macbook

if(onCluster==TRUE){
  def_wd <- "~/res/GEUV1Data/"
  setwd(def_wd)
  func_loc <-  "~/code/CompFunctions.R"
  #clust <- makeCluster(12) # Change this to match whatever number of cores requested on cluster or are on macbook
}else{
  def_wd <- "/Users/Scott/Documents/Dissertation Data/GEUV1Data/"
  def_wd2 <- "/Users/Scott/Documents/Dissertation/res/GEUV1Data"
  save_dir <- "/Users/Scott/Documents/Dissertation Data/GEUV1Data/"
  setwd(def_wd)
  func_loc <- "/Users/Scott/Documents/Dissertation/Code/CompFunctions.R"
}

source(func_loc)


geneGroupings <- c("Full", "genesMeanGeneMaxVarilrScaleLow","genesMeanGeneMaxVarilrScaleMedium","genesMeanGeneMaxVarilrScaleHigh","genesMeanGeneMaxVarilrScaleGr80","genesMeanGeneMaxVarilrScaleGr90",
                   "LowerThirdOverlap", "HighThirdOverlap", "genesMeanGeneMaxInfRVTPMGrEighty", "genesMeanGeneMaxInfRVTPMGrNinety")


#Main Paper Plots
#ROC Type Curve 1
GeneratePlots(infReps = "Boot", DRIMSeqFiltering = TRUE, SuppPlots = FALSE, UseRealSensRes = FALSE, 
              AllMethodsNonMissingpVals = FALSE, FilterGenes = FALSE,
              changestouse = 1, PlotTypes = 2, stdal = F, TwentySamplesOnly = FALSE, def_wd = def_wd, 
              def_wd2 = def_wd2, save_dir = save_dir, abDatasets = abDatasets, geneGroupings = geneGroupings,
              saveCompiledPDFs = TRUE)
rm(list = ls(pattern = "SensFpr"))
gc()

#ROC Type Curve 2,3,4
GeneratePlots(infReps = "Boot", DRIMSeqFiltering = TRUE, SuppPlots = FALSE, UseRealSensRes = FALSE, 
              AllMethodsNonMissingpVals = FALSE, FilterGenes = FALSE,
              changestouse = 2, PlotTypes = c(10,11), stdal = F, TwentySamplesOnly = FALSE, 
              def_wd = def_wd, def_wd2 = def_wd2, save_dir = save_dir, abDatasets = abDatasets, geneGroupings = geneGroupings,
              saveCompiledPDFs = TRUE)
rm(list = ls(pattern = "SensFpr"))
gc()
#End main paper plots


#Supplementary figures

#ROC Type Curve 5
GeneratePlots(infReps = "Boot", DRIMSeqFiltering = TRUE, SuppPlots = TRUE, UseRealSensRes = FALSE, 
              AllMethodsNonMissingpVals = FALSE, FilterGenes = FALSE,
              changestouse = 1, PlotTypes = 2, stdal = F, TwentySamplesOnly = TRUE, 
              def_wd = def_wd, def_wd2 = def_wd2, save_dir = save_dir, abDatasets = abDatasets, geneGroupings = geneGroupings,
              saveCompiledPDFs = TRUE)
rm(list = ls(pattern = "SensFpr"))
gc()

#ROC Type Curve 6
GeneratePlots(infReps = "Boot", DRIMSeqFiltering = TRUE, SuppPlots = TRUE, UseRealSensRes = FALSE, 
              AllMethodsNonMissingpVals = FALSE, FilterGenes = FALSE,
              changestouse = 2, PlotTypes = 10, stdal = F, TwentySamplesOnly = TRUE, 
              def_wd = def_wd, def_wd2 = def_wd2, save_dir = save_dir, abDatasets = abDatasets, geneGroupings = geneGroupings,
              saveCompiledPDFs = TRUE)
rm(list = ls(pattern = "SensFpr"))
gc()

#ROC Type Curve 7
GeneratePlots(infReps = "Boot", DRIMSeqFiltering = TRUE, SuppPlots = TRUE, UseRealSensRes = FALSE, 
              AllMethodsNonMissingpVals = FALSE, FilterGenes = FALSE,
              changestouse = 1, PlotTypes = 2, stdal = F, TwentySamplesOnly = FALSE, 
              def_wd = def_wd, def_wd2 = def_wd2, save_dir = save_dir, abDatasets = abDatasets, geneGroupings = geneGroupings,
              saveCompiledPDFs = TRUE)
rm(list = ls(pattern = "SensFpr"))
gc()

#ROC Type Curve 10, 11
GeneratePlots(infReps = "Boot", DRIMSeqFiltering = TRUE, SuppPlots = FALSE, UseRealSensRes = FALSE, 
              AllMethodsNonMissingpVals = FALSE, FilterGenes = FALSE,
              changestouse = 2, PlotTypes = c(9,10,11), stdal = F, TwentySamplesOnly = FALSE, 
              def_wd = def_wd, def_wd2 = def_wd2, save_dir = save_dir, abDatasets = abDatasets, geneGroupings = geneGroupings,
              saveCompiledPDFs = TRUE)
rm(list = ls(pattern = "SensFpr"))
gc()


#ROC Type Curve 12
GeneratePlots(infReps = "Boot", DRIMSeqFiltering = TRUE, SuppPlots = FALSE, UseRealSensRes = FALSE, 
              AllMethodsNonMissingpVals = FALSE, FilterGenes = FALSE,
              changestouse = 1, PlotTypes = 2, stdal = F, TwentySamplesOnly = TRUE, 
              def_wd = def_wd, def_wd2 = def_wd2, save_dir = save_dir, abDatasets = abDatasets, geneGroupings = geneGroupings,
              saveCompiledPDFs = TRUE)
rm(list = ls(pattern = "SensFpr"))
gc()


#ROC Type Curve 13,14,15
GeneratePlots(infReps = "Boot", DRIMSeqFiltering = TRUE, SuppPlots = FALSE, UseRealSensRes = FALSE, 
              AllMethodsNonMissingpVals = FALSE, FilterGenes = FALSE,
              changestouse = 2, PlotTypes = c(10,11), stdal = F, TwentySamplesOnly = TRUE, 
              def_wd = def_wd, def_wd2 = def_wd2, save_dir = save_dir, abDatasets = abDatasets, geneGroupings = geneGroupings,
              saveCompiledPDFs = TRUE)
rm(list = ls(pattern = "SensFpr"))
gc()

geneGroupings <- c("Full", "genesMeanGeneMaxVarilrScaleLow","genesMeanGeneMaxVarilrScaleMedium","genesMeanGeneMaxVarilrScaleHigh","genesMeanGeneMaxVarilrScaleGr80","genesMeanGeneMaxVarilrScaleGr90")
GeneratePlots(infReps = "Boot", DRIMSeqFiltering = TRUE, SuppPlots = FALSE, UseRealSensRes = FALSE, 
              AllMethodsNonMissingpVals = FALSE, FilterGenes = FALSE,
              changestouse = 2, PlotTypes = c(9,10,11), stdal = T, TwentySamplesOnly = FALSE, 
              def_wd = def_wd, def_wd2 = def_wd2, save_dir = save_dir, abDatasets = abDatasets, AbFromInfReps = TRUE, geneGroupings = geneGroupings,
              saveCompiledPDFs = TRUE)
rm(list = ls(pattern = "SensFpr"))

GeneratePlots(infReps = "Boot", DRIMSeqFiltering = TRUE, SuppPlots = FALSE, UseRealSensRes = FALSE, 
              AllMethodsNonMissingpVals = FALSE, FilterGenes = FALSE,
              changestouse = 2, PlotTypes = c(9,10,11), stdal = T, TwentySamplesOnly = TRUE, 
              def_wd = def_wd, def_wd2 = def_wd2, save_dir = save_dir, abDatasets = abDatasets, AbFromInfReps = TRUE, geneGroupings = geneGroupings,
              saveCompiledPDFs = TRUE)
rm(list = ls(pattern = "SensFpr"))


