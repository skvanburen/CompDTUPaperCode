Sys.sleep(sample(1:100,1))

#Make sure to use the module version of R
#Package conflicts with the module on the cluster prevent it from working properly
.libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")

#Save Inferential Replicates to R with one file per biological sample
  #These files will be saved within the SalmonFilesDir directory
  #If inferential replicates are not used you may skip the file and go straight to file (3)
  #Especially if the number of samples is large it will be very difficult to generate datasets for all inferential replicates without some form of
  #computing cluster available since the total amount of data gets quite large
  #In this case you may be forced to skip use of the inferential replicates and use CompDTU instead of CompDTUme
library(CompDTUReg)

#Set def_wd location to be the same from file (1)DataProcessing.R
def_wd <- "/pine/scr/s/k/skvanbur/TCGABRCAAnalysis/"

#Load tx2gene object that will be needed by future code
load(paste0(def_wd,"tx2gene.RData"))

#Load abGene file to get the key, which is used to construct list of files
load(paste0(def_wd, "key.RData"))

#SalmonFilesDir is the directory where the Salmon quantification results are saved
  #and where the sample specific inferential replicates will be saved

#SalmonFilesDir is the directory where the Salmon quantification results have already been saved
SalmonFilesDir <- paste0(def_wd, "SalmonBootSamps/")
setwd(SalmonFilesDir)

#Code needs to loop over the biological samples in some form, we provide sample
  #code for a slurm array but a loop could be used as well as long as curr_samp and curr_file_loc get assigned properly
  #Array val here needs to match the number of biological samples/replicates in the analysis
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))


#Set the number of observations used in the analysis (called nsamp even if they don't correspond to unique samples)
#Should match the number of rows in key from (1)
#(These results don't include normal subtypes so have 740 samples instead of 762)
nsamp <- 740

QuantFiles <- paste0(SalmonFilesDir, key$Sample, "/", "quant.sf")

#Names of each element in QuantFiles must be set to "Sample1", "Sample2", etc even if they are not unique biological samples
#This is because this is how the code will expect the columns to be named, just like for the key object
#tximport will name the columns in its created Salmon output object with these names and the code will expect them to be there
#in the format "Sample1", "Sample2", etc
names(QuantFiles) <- key$Identifier

#QuantFiles2 <- QuantFiles[gtools::mixedsort(names(QuantFiles))]
QuantFiles2 <- QuantFiles

array_val2 <- array_val

curr_samp <- names(QuantFiles2)[array_val2]
curr_file_loc <- QuantFiles2[array_val2]



#Set to TRUE if using Gibbs samples as inferential replicates, FALSE if using bootstrap samples
GibbsSamps <- FALSE

#Set value for countsFromAbundnace parameter for use with txImport
#Love (2018) (Swimming downstream: statistical analysis of differential transcript usage following Salmon quantification [version 3])
  #recommends "scaledTPM" for DTU analysis
  #See tximport for further options
countsFromAbundance <- "scaledTPM"




#Files will be saved to SalmonFilesDir/BootSamps for bootstrap samples or /GibbsSamps for Gibbs samples
#Should take something like 15 minutes to run per sample, leading to the recommendation that each sample
  #be run separately (such as part of an array job) is possible
SaveInfRepDataAsRData(curr_samp = curr_samp, tx2gene = tx2gene, curr_file_loc = curr_file_loc, GibbsSamps = GibbsSamps,
                     countsFromAbundance = countsFromAbundance, direc_to_save_res = SalmonFilesDir)
