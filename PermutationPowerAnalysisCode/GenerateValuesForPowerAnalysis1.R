#Generate the subset of samples used for Power Analysis 1 for the GEUV1Data and group combinations

if(version$nickname=="Planting of a Tree"){
  .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
}else{
  if(onCluster==TRUE){source("~/.Rprofile")}
}

library(data.table)
library(dplyr)
onCluster <- FALSE #Change to True if running on the cluster, False if running on macbook

if(onCluster==TRUE){
  def_wd <- "~/res/GEUV1Data/"
  setwd(def_wd)
  func_loc <-  "~/code/CompFunctions.R"
  #clust <- makeCluster(1) # Change this to match whatever number of cores requested on cluster or macbook
}else{
  def_wd <- "/Users/Scott/Documents/Dissertation Data/GEUV1Data/"
  setwd(def_wd)
  func_loc <- "/Users/Scott/Documents/Dissertation/Code/CompFunctions.R"
  #func_loc <- "/Users/Scott/Documents/Dissertation/res/SQCCData/SQCCFunctions.R"
  #clust <- makeCluster(1) # Change this to match whatever number of cores requested on cluster or macbook
  
}


load("abDatasets.RData")

source(func_loc)
load("tx2gene.RData")

#Now, generate the data frame that will contain the Part and Group Combination Values needed for the Compositional Power Analysis
total_samp <- 100
num_parts <- 100
num_combos <- 100

#Set the seed for the file
set.seed(3776)


#Now, generate the genes that will be changed
#Randomly choose half of all genes in the annotation that have >=2 Trans to change


TransNameCol <- which(colnames(tx2gene)=="tx_id")
tx2geneUniqueTemp <- tx2gene[,-TransNameCol]
tx2geneUniqueTemp2 <- distinct(tx2geneUniqueTemp)
tx2geneUniqueToUse <- subset(tx2geneUniqueTemp2, tx2geneUniqueTemp2$NTrans > 1)
gnames <- tx2geneUniqueToUse$gene_id
ngenes <- length(gnames)


sizgtochange <- ceiling(ngenes/2)
probsgtochange <- rep(sizgtochange/ngenes, ngenes)

genestochangeInd <- sample(ngenes, size = sizgtochange, replace = FALSE, prob = probsgtochange)

genestochange <- sort(gnames[genestochangeInd])




#For this test, I will choose 50 samples from each of two groups (CEU and GBR)
#levels(key$Condition)
sub_keyCEU <- subset(key, key$Condition=="CEU")

randnum1 <- sample(1:nrow(sub_keyCEU), 50, replace = FALSE)
sampsCEU <- sub_keyCEU$Identifier[randnum1]

sub_keyGBR <- subset(key, key$Condition=="GBR")

randnum2 <- sample(1:nrow(sub_keyGBR), 50, replace = FALSE)
sampsGBR <- sub_keyGBR$Identifier[randnum2]

sub_key <- subset(key, key$Identifier %in% sampsCEU | key$Identifier %in% sampsGBR)

sub_key$Condition <- relevel(factor(sub_key$Condition), ref=1)


save(sub_key, file = "KeyForPowerAnalysis1.RData")

#Specify change values to be considered use in the power analysis
  #This structure is depreciated such that these values do not actually control the change values but this is left in to ensure 
  #all code can still work
changes <- c(0.01, 0.02, 0.1, 0.25, 0.50, 1, 2, 4, 10, 50, 100)

temp1 <- rbindlist(lapply(1:num_combos, function(x, num_parts){return(data.frame(as.matrix(rep(x, num_parts), ncol=1)))}, num_parts = num_parts))

#Part, Group Combo
PowerAnalysisValuesTemp <- data.frame(rep(1:num_parts, num_combos), temp1)
colnames(PowerAnalysisValuesTemp) <- c("Part", "GroupNum")

poweranalyvaluesfunc <- function(x, PowerAnalysisValuesTemp){
  PowerAnalysisValuesTemp$change <- x
  return(PowerAnalysisValuesTemp)
}

PowerAnalysisValues <- rbindlist(lapply(changes, poweranalyvaluesfunc, PowerAnalysisValuesTemp = PowerAnalysisValuesTemp))



grpcombost <- matrix(0, nrow = num_combos, ncol = total_samp)
grpcombos <- data.frame(grpcombost)

for(i in 1:num_combos){
  grpcombos[i,] <- sample(sub_key$Condition, total_samp, prob = rep(1/total_samp, total_samp))
}

colnames(grpcombos) <- NULL

save(sub_key, PowerAnalysisValues, grpcombos, genestochange, file = "FilesForPowerAnalysis1.RData")

save(changes, file = "changes.RData")
#Reload the values from above to then save the DRIMSeq ones to make sure group combinations are the same

#Now, generate the necessary DRIMSeq values that don't need to vary based on Part
PowerAnalysisValues2 <- PowerAnalysisValues[,c("GroupNum", "change")]
DRIMSeqPowerAnalysisValues <- unique(PowerAnalysisValues2)
save(sub_key, DRIMSeqPowerAnalysisValues, grpcombos, file = "DRIMSeqFilesForPowerAnalysis1.RData")






