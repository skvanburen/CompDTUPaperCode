#Code to Create Observed/Real Data Analysis Table for Paper/Project 1


def_wd <- "/Users/Scott/Documents/Dissertation/"
def_wd2 <- "/Users/Scott/Documents/Dissertation Data/"

save_dir <- "/Users/Scott/Documents/Dissertation/Paper1/Tables"


setwd(def_wd)


tab_func <- function(fil1, fil2, fil3 = NULL, fil4 = NULL, fil5 = NULL, IncludeHardcodedBANDITSRes462Samples = FALSE){
  if(is.null(fil3) & is.null(fil4)){
    da <- data.frame(Method = c("CompDTU", "DRIMSeq"), stringsAsFactors = FALSE)
  }else if(is.null(fil4)){
    if(IncludeHardcodedBANDITSRes462Samples==TRUE){
      da <- data.frame(Method = c("CompDTU", "CompDTUme", "DRIMSeq", "BANDITS"), stringsAsFactors = FALSE)
    }else{
      da <- data.frame(Method = c("CompDTU", "CompDTUme", "DRIMSeq"), stringsAsFactors = FALSE)
    }

  }else{
    da <- data.frame(Method = c("CompDTU", "CompDTUme", "DRIMSeq", "RATsNoBoot", "RATsBoot", "BANDITS"), stringsAsFactors = FALSE)
  }

  rownames(da) <- da$Method

  #da$dataset <- rep(dataset, nrow(da))
  #da$filtering <- rep(filtering, nrow(da))

  #da$NSamp<- 462
  #da$NCond<- 5

  #Total number of genes that pass filtering can be found by looking at the length of abDatasetsNoOtherGroupsFiltered for each dataset
    
    #Computation time only includes time to fit the model and compute the pvalues, not the time it takes to load in datasets or return results, etc
    #This is accomplished for the Compositional methods by use of sum(CompObsPvals$comptgene) as the total computation time
    #This is true for all methods to keep comparisons as fair as possible
  load(fil1)
  #CompObsPvals <- CompositionalResSplit$CompPvals
  CompObsPvals <- CompositionalResSplit
  numg <- length(CompObsPvals$gene_id) - sum(is.na(CompObsPvals$pval_pillai))
  da["CompDTU", "$G$"] <- numg
  da["CompDTU", "$SG$"] <- sum(p.adjust(CompObsPvals$pval_pillai, method = "fdr") <=0.05, na.rm=T)
  #da["CompDTU", "NumNAPvals"] <-
  da["CompDTU", "$\\%SG$"] <- 100*sum(p.adjust(CompObsPvals$pval_pillai, method = "fdr") <=0.05, na.rm=T)/numg
  da["CompDTU", "Run Time (s)"] <- sum(CompObsPvals$comptgene)
  da["CompDTU", "Run Time Per Gene (s)"] <- mean(CompObsPvals$comptgene)

  #Load Necessary Files

  #The use of "DRIMSeqRunTime" to determine the total amount of computation time means and preprocessing time to load data/get the data into
    #the correct format or time for R functions written by SVB to return the results, etc are not included
    #For DRIMSeq, this time only includes time to calculate precisions with dmPrecision, fit the models with dmFit, and calculate tests with
    #dmTest
    #This keeps comparisons more fair between all methods
  load(fil2)
  numg <- length(DRIMPvals$pvalue) - sum(is.na(DRIMPvals$pvalue))
  da["DRIMSeq", "$G$"] <- numg
  da["DRIMSeq", "$SG$"] <- sum(p.adjust(DRIMPvals$pvalue, method = "fdr") <=0.05, na.rm=T)
  #da["DRIMSeq", "NumNAPvals"] <- sum(is.na(DRIMPvals$pvalue))
  da["DRIMSeq", "$\\%SG$"] <- 100*sum(p.adjust(DRIMPvals$pvalue, method = "fdr") <=0.05, na.rm=T)/numg
  da["DRIMSeq", "Run Time (s)"] <- DRIMSeqRunTime[3]
  da["DRIMSeq", "Run Time Per Gene (s)"] <- DRIMSeqRunTime[3]/numg

  if(is.null(fil3)){
    return(da)
  }
  load(fil3)
  
  CompModelingObsPvals <- CompositionalModelingResSplit
  #CompModelingObsPvals <- CompositionalModelingResSplit$CompPvals
  numg <- length(CompModelingObsPvals$gene_id) - sum(is.na(CompModelingObsPvals$pval_pillai))
  da["CompDTUme", "$G$"] <- numg
  da["CompDTUme", "$SG$"] <- sum(p.adjust(CompModelingObsPvals$pval_pillai, method = "fdr") <=0.05, na.rm=T)
  da["CompDTUme", "$\\%SG$"] <- 100*sum(p.adjust(CompModelingObsPvals$pval_pillai, method = "fdr") <=0.05, na.rm=T)/numg
  #da["CompDTUme", "NumNAPvals"] <- sum(is.na(CompModelingObsPvals$pval_pillai))
  da["CompDTUme", "Run Time (s)"] <- sum(CompModelingObsPvals$comptgene)
  da["CompDTUme", "Run Time Per Gene (s)"] <- mean(CompModelingObsPvals$comptgene)


  if(IncludeHardcodedBANDITSRes462Samples==TRUE){
    numg <- 7522
    da["BANDITS", "$G$"] <- NA
    
    #Include the calculation of the precision factor in the compuation time as well as the time to run the test_DTU call because the precision 
    #step comes "highly recommended" by the vignette
    BANDITSTotalTime <- 302400
    da["BANDITS", "Run Time (s)"] <- BANDITSTotalTime
    da["BANDITS", "Run Time Per Gene (s)"] <- BANDITSTotalTime/numg
  }

  #fil4 will be used for RATs if it is added
  if(is.null(fil4)){
    return(da)
  }
  load(fil4)
  RATsObsPvalsNoInfReps <- RATsObsRes$rats_res
  RATsObsPvalsInfReps <- RATsObsRes$rats_res_infReps

  ValsTemp2 <- RATsObsRes$rats_res$Genes
  
  numg <- length(ValsTemp2$pval) -  sum(is.na(ValsTemp2$pval))
  da["RATsNoBoot", "$G$"] <- numg
  
  #Modify RATs Pvalues to incorporate the "elig_fx" and "rep_reprod" post-hoc filters the paper recommends to make comparison as fair as possible to the full
    #recommendations from RATs 
    #This code/procedure was taken from Love (2018) "Swimdown" paper
    #In particular see the file "dtu_analysis" from swimdown code available at https://github.com/mikelove/swimdown/blob/master/dtu/dtu_analysis.R
  ValsTemp2$p <- ValsTemp2$pval
  ValsTemp2$p[!(ValsTemp2$elig_fx & ValsTemp2$rep_reprod)] <- 1
  ValsTemp2$padjfdr <- p.adjust(ValsTemp2$pval, method = "fdr")
  ValsTemp2$padjfdr[!(ValsTemp2$elig_fx & ValsTemp2$rep_reprod)] <- 1

  
  #Computation time only includes time to fit the model and compute the pvalues, not the time it takes to load in datasets or return results, etc
    #This is accomplished for RATs by only measuring time it takes to run call_DTU statement
    #This is true for all methods to keep comparisons as fair as possible
  da["RATsNoBoot", "$SG$"] <- sum(ValsTemp2$padjfdr <=0.05, na.rm=T)
  #da["RATsNoBoot", "NumNAPvals"] <- length(subset(RATsObsRes$rats_res$Genes$pval, is.na(RATsObsRes$rats_res$Genes$pval)))
  da["RATsNoBoot", "$\\%SG$"] <- 100*sum(ValsTemp2$padjfdr <=0.05, na.rm=T)/numg
  da["RATsNoBoot", "Run Time (s)"] <- RATsObsRes$RatsCompTime[3]
  da["RATsNoBoot", "Run Time Per Gene (s)"] <- RATsObsRes$RatsCompTime[3]/numg

  
  
  ValsTemp3 <- RATsObsRes$rats_res_infReps$Genes

  numg <- length(ValsTemp3$pval) -  sum(is.na(ValsTemp3$pval))
  da["RATsBoot", "$G$"] <- numg
  
  #Modify RATs Pvalues to incorporate the "elig_fx" and "rep_reprod" post-hoc filters the paper recommends to make comparison as fair as possible to the full
  #recommendations from RATs 
  #This code/procedure was taken from Love (2018) "Swimdown" paper
  #In particular see the file "dtu_analysis" from swimdown code available at https://github.com/mikelove/swimdown/blob/master/dtu/dtu_analysis.R
  ValsTemp3$p <- ValsTemp3$pval
  ValsTemp3$p[!(ValsTemp3$elig_fx & ValsTemp3$rep_reprod)] <- 1
  ValsTemp3$padjfdr <- p.adjust(ValsTemp3$pval, method = "fdr")
  ValsTemp3$padjfdr[!(ValsTemp3$elig_fx & ValsTemp3$rep_reprod)] <- 1
  #da["RATsBoot", "NumNAPvals"] <- length(subset(RATsObsRes$rats_res_infReps$Genes$pval, is.na(RATsObsRes$rats_res_infReps$Genes$pval)))
  
  #Computation time only includes time to fit the model and compute the pvalues, not the time it takes to load in datasets or return results, etc
    #This is accomplished for RATs by only measuring time it takes to run call_DTU statement
    #This is true for all methods to keep comparisons as fair as possible
  da["RATsBoot", "$SG$"] <- sum(ValsTemp3$padjfdr <=0.05, na.rm=T)
  da["RATsBoot", "$\\%SG$"] <- 100*sum(ValsTemp3$padjfdr <=0.05, na.rm=T)/numg
  da["RATsBoot", "Run Time (s)"] <- RATsObsRes$RatsCompTimeInfReps[3]
  da["RATsBoot", "Run Time Per Gene (s)"] <- RATsObsRes$RatsCompTimeInfReps[3]/numg
  
  if(is.null(fil5)){
    return(da)
  }
  
  load(fil5)
  numg <- sum(!is.na(GeneLevelSignRes$adj.p.values))
  da["BANDITS", "$G$"] <- numg
  da["BANDITS", "$SG$"] <- sum(GeneLevelSignRes$adj.p.values <=0.05, na.rm=T)
  da["BANDITS", "$\\%SG$"] <- 100*sum(GeneLevelSignRes$adj.p.values <=0.05, na.rm=T)/numg
  
  #Include the calculation of the precision factor in the compuation time as well as the time to run the test_DTU call because the precision 
    #step comes "highly recommended" by the vignette
  BANDITSTotalTime <- TimeToCreateBANDITSPrecision[3] + TimeToRunBANDITS[3]
  da["BANDITS", "Run Time (s)"] <- BANDITSTotalTime
  da["BANDITS", "Run Time Per Gene (s)"] <- BANDITSTotalTime/numg
  
  return(da)
}

fil_dir <- "/Users/Scott/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/SQCC/"
fil1 <- paste0(fil_dir, "CompositionalResSplit.RData")
fil2 <- paste0(fil_dir, "DRIMSeqRes.RData")
fil3 <- paste0(fil_dir, "CompositionalModelingResSplit.RData")
fil4 <- paste0(fil_dir, "RATsObsRes.RData")
ObsAnalysisTableSQCC <- tab_func(fil1 = fil1, fil2 = fil2, fil3 = fil3, fil4 = fil4)


save(ObsAnalysisTableSQCC, file = paste0(save_dir, "/", "Paper1CompTimeTableSQCC.RData"))

fil_dir <- "/Users/Scott/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/GEUV1/"
fil1 <- paste0(fil_dir, "CompositionalResSplit.RData")
fil2 <- paste0(fil_dir, "DRIMSeqRes.RData")
fil3 <- paste0(fil_dir, "CompositionalModelingResSplit.RData")
ObsAnalysisTableGEUV1 <- tab_func(fil1 = fil1, fil2 = fil2, fil3 = fil3)

save(ObsAnalysisTableGEUV1, file = paste0(save_dir, "/", "Paper1CompTimeTableGEUV1.RData"))


fil_dir <- "/Users/Scott/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/SQCC/"
fil1 <- paste0(fil_dir, "CompositionalResSplitDRIMSeqFiltering.RData")
fil2 <- paste0(fil_dir, "DRIMSeqResDRIMSeqFiltering.RData")
fil3 <- paste0(fil_dir, "CompositionalModelingResSplitDRIMSeqFiltering.RData")
fil4 <- paste0(fil_dir, "RATsObsResDRIMSeqFiltering.RData")
ObsAnalysisTableSQCCDRIMSeqFiltering <- tab_func(fil1 = fil1, fil2 = fil2, fil3 = fil3, fil4 = fil4)


save(ObsAnalysisTableSQCCDRIMSeqFiltering, file = paste0(save_dir, "/", "Paper1CompTimeTableSQCCDRIMSeqFiltering.RData"))

fil_dir <- "/Users/Scott/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/GEUV1/"
fil1 <- paste0(fil_dir, "CompositionalResSplitDRIMSeqFiltering.RData")
fil2 <- paste0(fil_dir, "DRIMSeqResDRIMSeqFiltering.RData")
fil3 <- paste0(fil_dir, "CompositionalModelingResSplitDRIMSeqFiltering.RData")
ObsAnalysisTableGEUV1DRIMSeqFiltering <- tab_func(fil1 = fil1, fil2 = fil2, fil3 = fil3, IncludeHardcodedBANDITSRes462Samples = TRUE)

save(ObsAnalysisTableGEUV1DRIMSeqFiltering, file = paste0(save_dir, "/", "Paper1CompTimeTableGEUV1DRIMSeqFiltering.RData"))



fil_dir <- "/Users/Scott/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/GEUV1TwentySamples/"
fil1 <- paste0(fil_dir, "CompositionalResSplitDRIMSeqFiltering.RData")
fil2 <- paste0(fil_dir, "DRIMSeqResDRIMSeqFiltering.RData")
fil3 <- paste0(fil_dir, "CompositionalModelingResSplitDRIMSeqFiltering.RData")
fil4 <- paste0(fil_dir, "RATsObsResDRIMSeqFiltering.RData")
fil5 <- paste0(fil_dir, "BANDITSResDRIMSeqFiltering.RData")
ObsAnalysisTableGEUV1DRIMSeqFilteringTwentySamples <- tab_func(fil1 = fil1, fil2 = fil2, fil3 = fil3, fil4 = fil4, fil5 = fil5)

save(ObsAnalysisTableGEUV1DRIMSeqFilteringTwentySamples, file = paste0(save_dir, "/", "Paper1CompTimeTableGEUV1DRIMSeqFilteringTwentySamples.RData"))


ObsAnalysisTableGEUV1DRIMSeqFiltering$SampleSize <- 462
ObsAnalysisTableGEUV1DRIMSeqFilteringTwentySamples$SampleSize <- 20

#, c("", NA, NA, NA, NA, NA)
ObsAnalysisTableGEUV1DRIMSeqFilteringTotalT <- rbind(ObsAnalysisTableGEUV1DRIMSeqFilteringTwentySamples, ObsAnalysisTableGEUV1DRIMSeqFiltering)

#Results including G and SG columns
#ObsAnalysisTableGEUV1DRIMSeqFilteringTotal <- ObsAnalysisTableGEUV1DRIMSeqFilteringTotalT[,c(1,7,2,3,4,5,6)]
ObsAnalysisTableGEUV1DRIMSeqFilteringTotal <- ObsAnalysisTableGEUV1DRIMSeqFilteringTotalT[,c(1,7,2,5,6)]

save(ObsAnalysisTableGEUV1DRIMSeqFilteringTotal, file = paste0(save_dir, "/", "Paper1CompTimeTableGEUV1DRIMSeqFilteringTotal.RData"))




# tab_func2 <- function(fil1, fil2){
#   da <- data.frame(method = c("Comp", "CompME"))
#
#   rownames(da) <- da$method
#
#   #da$dataset <- rep(dataset, nrow(da))
#   #da$filtering <- rep(filtering, nrow(da))
#
#   #da$NSamp<- 462
#   #da$NCond<- 5
#
#
#   load(fil1)
#   CompObsPvals <- CompositionalResSplit
#   numg <- length(CompObsPvals$gene_id) - sum(is.na(CompObsPvals$pval_pillai))
#   da["Comp", "$G$"] <- numg
#   #da["Comp", "NumNAPvals"] <-
#   da["Comp", "RunTime (s)"] <- TotalTimeCompSplit[3]
#   da["Comp", "RunTime Per Gene (s)"] <- TotalTimeCompSplit[3]/numg
#
#
#   load(fil2)
#   CompModelingObsPvals <- CompositionalModelingResSplit
#   numg <- length(CompModelingObsPvals$gene_id) - sum(is.na(CompModelingObsPvals$pval_pillai))
#   da["CompME", "$G$"] <- numg
#   #da["CompME", "NumNAPvals"] <- sum(is.na(CompModelingObsPvals$pval_pillai))
#   da["CompME", "RunTime (s)"] <- TotalTimeModelingSplit[3]
#   da["CompME", "RunTime Per Gene (s)"] <- TotalTimeModelingSplit[3]/numg
#
#
#   return(da)
# }
#
#
# fil_dir <- "/Users/Scott/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/SQCC/"
# fil1 <- paste0(fil_dir, "CompositionalResSplit.RData")
# fil2 <- paste0(fil_dir, "CompositionalModelingResSplit.RData")
# ObsAnalysisTableSQCCSplit <- tab_func2(fil1 = fil1, fil2 = fil2)
#
#
# fil_dir <- "/Users/Scott/Documents/Dissertation/Paper1/Tables/ComputationTimeTable/GEUV1/"
# fil1 <- paste0(fil_dir, "CompositionalResSplit.RData")
# fil2 <- paste0(fil_dir, "CompositionalModelingResSplit.RData")
# ObsAnalysisTableGEUV1Split <- tab_func2(fil1 = fil1, fil2 = fil2)
#
#
# ComparisonTableSQCC <- rbind(ObsAnalysisTableSQCC[c("Comp", "CompME"),],ObsAnalysisTableSQCCSplit[c("Comp", "CompME"),])
# ComparisonTableSQCC[,"Run1GeneAtATime"] <- c("No", "No", "Yes", "Yes")
# ComparisonTableSQCC$Dataset <- "SQCC"
#
# ComparisonTableGEUV1 <- rbind(ObsAnalysisTableGEUV1[c("Comp"),],ObsAnalysisTableGEUV1Split[c("Comp", "CompME"),])
# ComparisonTableGEUV1[,"Run1GeneAtATime"] <- c("No", "Yes", "Yes")
# ComparisonTableGEUV1$Dataset <- "GEUV1"
#
# ComparisonTableT <- rbind(ComparisonTableSQCC, ComparisonTableGEUV1)
# rownames(ComparisonTableT) <- NULL
#
# ComparisonTable <- ComparisonTableT[c(1,3,2,4,5,6,7),]
# rownames(ComparisonTable) <- NULL


