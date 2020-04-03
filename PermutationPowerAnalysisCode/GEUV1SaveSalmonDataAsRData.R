#Read in all gibbs replicates and quant data and save as R object with Gibbs reps and tximport object
#It seems to be faster than reading in using read_tsv, and is much faster reading this in than reading
# the raw data each time
  #I have just been running this code in interactive mode for now, no need to submit it as a job for now

if(version$nickname=="Planting of a Tree"){
  .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
}else{
  if(onCluster==TRUE){source("~/.Rprofile")}
}

library(gtools)
library(readr)
library(tximport)
library(rjson)

onCluster <- TRUE
if(onCluster==TRUE){
  dir1 <- "~/res/GEUV1Data/"
  dir2 <- "/pine/scr/s/k/skvanbur/GEUV1/"
}else{
  dir1 <- "/Users/Scott/Documents/Dissertation/res/GEUV1Data/"
}
  
#Set countsFromAbundance parameter, which controls how the counts are estimated 
  #(see tximport help documentation for more information on the differences this parameter makes)
countsFromAbundance <- "no"
#countsFromAbundance <- "scaledTPM"

#Set to True if you've drawn Gibbs samples and FALSE if you've drawn bootstrap samples
GibbsSamps <- TRUE

#Set this to TRUE if all samples should have been included in the quantification
  #Set to FALSE if you only ran some of the samples intentionally (as an example)ie if you only ran first 10 samples intentionally, etc
AllSampsShouldHaveBeRun <- TRUE

if(GibbsSamps==TRUE){
  def_wd <- "/pine/scr/s/k/skvanbur/GEUV1/Salmon"
}else{
  def_wd <- "/pine/scr/s/k/skvanbur/GEUV1/SalmonBootSamps"
}
setwd(def_wd)

#base_dir <- "/pine/scr/s/k/skvanbur/GEUV1/Data/"

#Create key matrix that contains matches samples to conditions and identifiers
#Load file that has the identifying information
E_GEUV_1_sdrf <- read_delim(paste0(dir1, "E-GEUV-1.sdrf.txt"), 
                            "\t", escape_double = FALSE, trim_ws = TRUE)

#Take true sample name, condition, as well as the ENA-RUN., which is the unique idenfitier used to download data and do quantification
cols1 <- which(colnames(E_GEUV_1_sdrf) %in% c("Comment[ENA_SAMPLE]", "Comment[ENA_RUN]", "Characteristics[population]"))
temp1 <- E_GEUV_1_sdrf[,cols1]

#Are 2 identical rows for each sample (1 for each of the 2 paired end fastq files), so only take one of each
#Should give a total of 462 rows for this dataset
temp2 <- unique(temp1)
temp3 <- temp2[mixedorder(temp2$`Comment[ENA_RUN]`),]
#Create simpler identifier that is just Sample1, ..., Sample462
#Samples are ordered by ERA run number, Sample Number because this is how files are downloaded
temp3$Identifier <- paste0("Sample", 1:nrow(temp3))

colnames(temp3) <- c("Sample", "Condition", "ENARun", "Identifier")
key <- as.data.frame(temp3, stringsAsFactors = FALSE)

#Set Condition to be a factor to help code downstream
key$Condition <- relevel(as.factor(key$Condition), ref = 1)

#Extract file location of all of the quantification files of interest, the .sf files
#Also need to carefully name each QuantFile corresponding to the identifiers created above
QuantFiles <- mixedsort(list.files(pattern = ".sf", recursive = TRUE, full.names = TRUE))
subnames <- dir(recursive=F, pattern = "ERR")
for(i in 1:length(subnames)){
  currname <- subnames[i]
  wherefound <- grep(currname, QuantFiles)
  names(QuantFiles)[wherefound] <- key$Identifier[key$ENARun==currname]
}

QuantFiles2 <- QuantFiles[mixedsort(names(QuantFiles))]

#Check and make sure there are the right number of .sf files, if not something went wrong with the 
  #download/quantification, etc and you need to check on that before continuing
if(AllSampsShouldHaveBeRun==TRUE){
  stopifnot(length(QuantFiles)==length(key$Sample))
}

#For this data, don't save the Gibbs samples to QuantSalmon because it would take too much space and some of them have issues
  #Use countsFromAbundance = "scaledTPM" as these are the counts recommended for use with DRIMSeq
  #See Love, Soneson, Patro Sept 2018 (doi: 10.12688/f1000research.15398.2)


QuantSalmon <- tximport(QuantFiles2, type = "salmon", txOut = TRUE, ignoreTxVersion = FALSE, 
                        countsFromAbundance = countsFromAbundance, dropInfReps = T)
fulltransnames <- rownames(QuantSalmon$abundance) #transcript names

QuantSalmonAbRowMedianInfRep <- tximport(QuantFiles2, type = "salmon", txOut = TRUE, ignoreTxVersion = FALSE, 
                        countsFromAbundance = countsFromAbundance, dropInfReps = F, infRepStat = matrixStats::rowMedians)
QuantSalmonAbRowMedianInfRep$infReps <- NULL
gc()
fulltransnames <- rownames(QuantSalmonAbRowMedianInfRep$abundance) #transcript names




QuantSalmonAbRowMeanInfRep <- tximport(QuantFiles2, type = "salmon", txOut = TRUE, ignoreTxVersion = FALSE, 
                                       countsFromAbundance = countsFromAbundance, dropInfReps = F, infRepStat = rowMeans)
QuantSalmonAbRowMeanInfRep$infReps <- NULL
gc()
fulltransnames2 <- rownames(QuantSalmonAbRowMeanInfRep$abundance) #transcript names
sum(fulltransnames2==fulltransnames)

#Now, also check which gibbs samples failed and save an object containing that information for use downstream

failedgibbs <- c()
for(i in 1:length(QuantFiles2)){
  QuantSalmon2 <- tryCatch(tximport(QuantFiles2[i], type = "salmon", txOut = TRUE, ignoreTxVersion = FALSE, 
                                    countsFromAbundance = countsFromAbundance, dropInfReps =F),error=function(e){})
  if(is.null(QuantSalmon2)){
    failedgibbs <- c(failedgibbs, i)
  }
}

if(GibbsSamps==FALSE){
  failedboots <- failedgibbs
}
if(countsFromAbundance =="scaledTPM"){

  if(GibbsSamps==TRUE){
    save(QuantSalmon, QuantFiles, key, fulltransnames, failedgibbs, countsFromAbundance,
         file = paste0(dir2,"Salmon/SalmonDataCountsScaledTPM.RData"))
    
    save(QuantSalmonAbRowMedianInfRep, QuantFiles, key, fulltransnames, failedgibbs, countsFromAbundance,
         file = paste0(dir2,"Salmon/SalmonDataCountsScaledTPMAbRowMedianInfRep.RData"))
    
    save(QuantSalmonAbRowMeanInfRep, QuantFiles, key, fulltransnames, failedgibbs, countsFromAbundance,
         file = paste0(dir2,"Salmon/SalmonDataCountsScaledTPMAbRowMeanInfRep.RData"))
    
    save(failedgibbs, countsFromAbundance, file = paste0(dir2, "failedgibbssampsCountsScaledTPM.RData"))
  }else{
    save(QuantSalmon, QuantFiles, key, fulltransnames, failedboots, countsFromAbundance,
         file = paste0(dir2,"SalmonBootSamps/SalmonDataCountsScaledTPM.RData"))
    
    save(QuantSalmonAbRowMedianInfRep, QuantFiles, key, fulltransnames, failedgibbs, countsFromAbundance,
         file = paste0(dir2,"SalmonBootSamps/SalmonDataCountsScaledTPMAbRowMedianInfRep.RData"))
    
    save(QuantSalmonAbRowMeanInfRep, QuantFiles, key, fulltransnames, failedgibbs, countsFromAbundance,
         file = paste0(dir2,"SalmonBootSamps/SalmonDataCountsScaledTPMAbRowMeanInfRep.RData"))
    
    save(failedboots, countsFromAbundance, file = paste0(dir2, "failedbootsampsCountsScaledTPM.RData"))
  }
  
}else if(countsFromAbundance =="no"){
  
  if(GibbsSamps==TRUE){
    save(QuantSalmon, QuantFiles, key, fulltransnames, failedgibbs, countsFromAbundance, file = paste0(dir2,"Salmon/SalmonData.RData"))
    
    save(QuantSalmonAbRowMedianInfRep, QuantFiles, key, fulltransnames, failedgibbs, countsFromAbundance,
         file = paste0(dir2,"Salmon/SalmonDataAbRowMedianInfRep.RData"))
    
    save(QuantSalmonAbRowMeanInfRep, QuantFiles, key, fulltransnames, failedgibbs, countsFromAbundance,
         file = paste0(dir2,"Salmon/SalmonDataAbRowMeanInfRep.RData"))
    
    save(failedgibbs, countsFromAbundance, file = paste0(dir2, "failedgibbssamps.RData"))
  }else{
    save(QuantSalmon, QuantFiles, key, fulltransnames, failedboots, countsFromAbundance, file = paste0(dir2,"SalmonBootSamps/SalmonData.RData"))
    
    save(QuantSalmonAbRowMedianInfRep, QuantFiles, key, fulltransnames, failedgibbs, countsFromAbundance,
         file = paste0(dir2,"SalmonBootSamps/SalmonDataAbRowMedianInfRep.RData"))
    
    save(QuantSalmonAbRowMeanInfRep, QuantFiles, key, fulltransnames, failedgibbs, countsFromAbundance,
         file = paste0(dir2,"SalmonBootSamps/SalmonDataAbRowMeanInfRep.RData"))
    
    
    save(failedboots, countsFromAbundance, file = paste0(dir2, "failedbootsamps.RData"))
  }
  
}



