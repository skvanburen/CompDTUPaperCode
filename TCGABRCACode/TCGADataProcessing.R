#Make sure to use the module version of R
#Package conflicts with the module on the cluster prevent it from working properly
.libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")

#Code to process data quantified from Salmon
#Make sure to run this file with the custom R installation since the custom packages wouldn't install properly otherwise
library(CompDTUReg)
#library(TCGAutils)
library(readr)


#def_wd is the top level directory where files will be saved
  #Modify it to wherever you would like the files to be saved
#def_wd <- "~/res/TCGABRCAAnalysis/"
def_wd <- "/pine/scr/s/k/skvanbur/TCGABRCAAnalysis/"
if(!dir.exists(def_wd)){
  dir.create(def_wd)
}
setwd(def_wd)

#SalmonFilesDir is the directory where the Salmon quantification results have already been saved
#See the Github repo CompDTURegSampleData located at https://github.com/skvanburen/CompDTURegSampleData to download sample data for
  #10 replicates with 100 bootstrap samples each
SalmonFilesDir <- paste0(def_wd, "SalmonBootSamps/")

#Specify location of annotation to use in maketx2gene
  #GENCODE annotations can be downloaded from https://www.gencodegenes.org/human/
  #We used the annotations for the reference chromosomes only, and used version 27 for the results in the paper
txdb_loc <- "~/gencode.v27.annotation.gtf.gz"



#Construct a cluster from parallel package for possible use later.  For no parallelization, use makeCluster(1)
clust <- parallel::makeCluster(1)

#Build tx2gene data.frame matching transcripts to genes using GENCODE if it doesn't exist
#Will take no more than a few minutes to run
if(!file.exists("tx2gene.RData")){
  maketx2gene(txdb_loc = txdb_loc)
}

#Load tx2gene object that was just created
load("tx2gene.RData")


#Read in Salmon Files and save results (without inferential replicates for now, as trying to save all of them in one file
  #can get prohibitively large)
setwd(SalmonFilesDir)

#Set value for countsFromAbundnace parameter for use with tximport
  #Love (2018) (Swimming downstream: statistical analysis of differential transcript usage following Salmon quantification [version 3])
  #recommends "scaledTPM" for DTU analysis
  #See tximport for further options
countsFromAbundance <- "scaledTPM"

file_manifest <- data.frame(read_delim("~/res/TCGABRCAAnalysis/FilesForGDCDownloader/final_manifest_all.txt", delim = "\t"))
file_manifest$case_submitter_id <- as.character(lapply(file_manifest$barcode, function(x){substr(x, 1, 12)}))

file_manifest$NumCase <- as.numeric(lapply(file_manifest$case_id, function(x){sum(file_manifest$case_id==x)}))

bar_to_remove <- c()
uniq_cases <- unique(file_manifest$case_id)
for(i in 1:length(uniq_cases)){
  curr_case <- uniq_cases[i]
  curr_r <- file_manifest[file_manifest$case_id==curr_case, , drop = FALSE]
  
  if(nrow(curr_r)==1){
    next
  }else{
    print(paste0("Duplicate files found for case ", curr_case))
    map_rates <- c()
    for(j in 1:nrow(curr_r)){
      curr_s_r <- curr_r[j, , drop = FALSE]
      log_fil_loc <- paste0("/pine/scr/s/k/skvanbur/TCGABRCAAnalysis/SalmonBootSamps/", curr_s_r$barcode, "/logs/salmon_quant.log")
      map_r_t <- grep("Mapping rate", readLines(log_fil_loc), value = TRUE)
      map_r_t2 <- substr(map_r_t, nchar(map_r_t) - 7, nchar(map_r_t) - 1)
      map_rates <- c(map_rates, as.numeric(map_r_t2))
    }
    
    bar_to_keep <- curr_r$barcode[which(map_rates==max(map_rates))]
    print(paste0("Mapping rates found are: "))
    print(map_rates)
    print(paste0("The highest mapping rate is ", max(map_rates), " found for barcode ", bar_to_keep))
    
    bar_to_remove <- c(bar_to_remove, subset(curr_r$barcode, curr_r$barcode!=bar_to_keep))
  }
}

file_manifest2 <- subset(file_manifest, !(file_manifest$barcode%in%bar_to_remove))

#Recalculate the NumCase and make sure there are still no duplicates
file_manifest2$NumCaseAfterRemoval <- as.numeric(lapply(file_manifest2$case_id, function(x){sum(file_manifest2$case_id==x)}))

clinical_info <- data.frame(read_delim("~/res/TCGABRCAAnalysis/clinical.tsv", delim = "\t"))
clinical_info2 <- clinical_info[,c("case_id", "primary_diagnosis")]


clinical_info3 <- unique(clinical_info2)
#Merge in the case_submitter_ids for the samples used
#Case id 57a1604c-60b7-4b30-a75e-f70939532c5c has no clinical information and there are 7 participants with clinical into 
  #And no RNA Seq available
#Final result then is 1090 total results
file_manifest3 <- merge(file_manifest2, clinical_info3, by.x = "case_id", by.y = "case_id")
#file_manifest$CaseID <- UUIDtoUUID(file_manifest$id, to_type = "case_id")
BarcodesUsedT <- gtools::mixedsort(file_manifest3$barcode)

file_manifest4 <- subset(file_manifest3, file_manifest3$barcode!="TCGA-A2-A0EX-01A-21R-A034-07")
#The BAM file for barcode was corrupted and couldn't be used when the data was initially created
#Make sure to remove this sample from the full annotation and add it to the end of the Sample list if you ever are able to get results working for it
#BarcodesUsed <- subset(BarcodesUsedT, BarcodesUsedT!="TCGA-A2-A0EX-01A-21R-A034-07")


#Now, restrict to only the most common primary diagnosis type "Infiltrating duct carcinoma, NOS"
file_manifest5 <- subset(file_manifest4, file_manifest4$primary_diagnosis=="Infiltrating duct carcinoma, NOS")

#Now, remove and patients with the normal subtype (only 22 out of the 762 have a normal subtype and it often isn't really used)
file_manifest6 <- subset(file_manifest5, file_manifest5$BRCA_Subtype_PAM50!="Normal")

BarcodesUsed <- file_manifest6$barcode
Cond_to_use <- file_manifest6$BRCA_Subtype_PAM50

#Create key matrix that contains matches samples to conditions and identifiers
#Code later will be expecting key to have columns "Sample", "Condition", and "Identifier"
#"Sample" has a unique identifier for the particular sample/replicate
#"Condition" is a factor variable that corresponds to the different groups/conditions to conduct DTU analysis on
#With "Identifier" containing names as "Sample1", "Sample2", etc even if data isn't corresponding to unique biological samples
#because this is how future code will expect key to be constructed
key <- data.frame(BarcodesUsed, Cond_to_use)
key <- as.data.frame(key, stringsAsFactors = FALSE)
key[,3] <- paste0("Sample", 1:nrow(key))
colnames(key) <- c("Sample", "Condition", "Identifier")
key$Condition <- relevel(as.factor(key$Condition), ref = 1)

#List of Salmon quantification files for data to be used in the analysis, (files end in .sf)
  #Results may have been quantified for other samples but restricting to the primary diagnosis being 
  #Infiltrating duct carcinoma, NOS and using the subtype as the condition

QuantFiles <- paste0(SalmonFilesDir, key$Sample, "/", "quant.sf")
#QuantFiles <- gtools::mixedsort(list.files(pattern = ".sf", recursive = TRUE, full.names = TRUE))


#Names of each element in QuantFiles must be set to "Sample1", "Sample2", etc even if they are not unique biological samples
  #This is because this is how the code will expect the columns to be named
  #tximport will name the columns in its created Salmon output object with these names and the code will expect them to be there
  #in the format "Sample1", "Sample2", etc
names(QuantFiles) <- key$Identifier





#Use tximport package to load in the results that have been pre-quantified by salmon
  #Drop inferential replicates for now, as they will be read in in file (2) if used

if(!file.exists("SalmonData.RData")){
  QuantSalmon <- tximport::tximport(QuantFiles, type = "salmon", txOut = TRUE, ignoreTxVersion = FALSE,
                                    countsFromAbundance = countsFromAbundance, dropInfReps = TRUE)
  fulltransnames <- rownames(QuantSalmon$abundance) #transcript names
  
  
  #Save the tximport object that contains the results of the salmon quantification as an r quantification file
  save(QuantSalmon, QuantFiles, key, fulltransnames, countsFromAbundance,
       file = paste0(SalmonFilesDir, "SalmonData.RData"))
}else{
  load("SalmonData.RData")
}

print(gc())

#Reset working directory to value specified by def_wd
setwd(def_wd)

#Save file containing the key
save(key, file = "key.RData")


  
  #This code will save data into two formats that will be needed later
  #See the help file for sumToGene for more information
  #Files will save to the current working directory
sumToGene(QuantSalmon = QuantSalmon, tx2gene = tx2gene, clust = clust, key = key,
          countsFromAbundance = countsFromAbundance)

rm(QuantSalmon)
print(gc())

##################################################################
#Filter Genes using DRIMSeq's Filtering Approach
#See DRIMSeq Paper for more information
##################################################################
#Load dataframe that contains quantification results for each transcript(rows) and sample (column) along with other information
  #Output by sumToGene above
load("abGene.RData")

#Load list of dataframes that separates data by gene, with element in the list being a data frame of expression for that gene
  #Output by sumToGene above
  #Within this data frame, rows are samples and columns are transcripts
load("abDatasets.RData")



if(countsFromAbundance=="scaledTPM" | countsFromAbundance=="lengthScaledTPM"){
  load("cntGenecntsScaledTPM.RData")
  load("cntDatasetsNoOtherGroupscntsScaledTPM.RData")
}else if(countsFromAbundance=="no"){
  load("cntGene.RData")
  load("cntDatasetsNoOtherGroups.RData")
}



#The filtering values below can be modified to make the filtering more or less strict
  #These are the default values used in
  #Love et al (2018) (Swimming downstream: statistical analysis of differential transcript usage following Salmon quantification [version 3])
  #See the help file for the dmFilter function for more information about each filtering parameter
  #Even if no filtering is desired, make sure to run the DRIMSeqFilter function to create all expected datasets in the proper places
  #In particular for minimal filtering consider setting all values to 0, corresponding to only removing features (transcripts) with zero expression across all samples 
  #and genes with only one non-zero feature (since DTU analysis cannot be performed if a gene has only one transcript).


#Sample size of smallest condition
n.small <- min(table(key$Condition))
n <- nrow(key)

min_samps_feature_expr <- n.small
min_feature_expr <- 10

min_samps_feature_prop <- n.small
min_feature_prop <- 0.10

min_samps_gene_expr <- n
min_gene_expr <- 10

#Set the list of samples to use to determine which genes/transcripts pass filtering based on
  #Set to a character vector in the form of key$Identifier, for example ("Sample1", "Sample5", etc)
sampstouse <- key$Identifier

#Should take somewhere around 5 minutes to run with 10 samples
DRIMSeqFilter(abGene = abGene, cntGene = cntGene, key = key, min_samps_feature_expr = min_samps_feature_expr, min_feature_expr = min_feature_expr,
              min_samps_feature_prop = min_samps_feature_prop, min_feature_prop = min_feature_prop, min_samps_gene_expr = min_samps_gene_expr,
              min_gene_expr = min_gene_expr, sampstouse = sampstouse, tx2gene = tx2gene, countsFromAbundance = countsFromAbundance)

print(gc())


