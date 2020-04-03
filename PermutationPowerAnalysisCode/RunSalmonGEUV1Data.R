#Put in a sleep statement just to prevent too many jobs from trying to access the index file at the same time as this can cause the jobs to fail
Sys.sleep(sample(1:100, 1))

#Val should run from 1 to 462 in this file, since there are paired end reads for each  of 462 samples that will be run together in Salmon, resulting in a total of 462*2=924 rows
val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

reRunToGenerateEqClassFiles <- TRUE

if(version$nickname=="Planting of a Tree"){
  .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
}else{
  if(onCluster==TRUE){source("~/.Rprofile")}
}

library(readr)
onCluster <- TRUE

#Set to True to draw Gibbs samples and FALSE to draw Bootstrap samples
#Don't forget to change the log file location in RunSalmonGEUV1Data.sh file also
GibbsSamps <- FALSE

#Set the number of inf reps (either Gibbs Samples or Bootstrap samples, but not both)
ninf <- 100

#Set thinning factor (only used if drawing Gibbs samples)
thinf <- 100

if(onCluster==TRUE){
  dir1 <- "~/res/GEUV1Data/"
}else{
  dir1 <- "/Users/Scott/Documents/Dissertation/res/GEUV1Data/"
}

base_dir <- "/pine/scr/s/k/skvanbur/GEUV1/Data/"


if(GibbsSamps==TRUE){
  output_dir <- "/pine/scr/s/k/skvanbur/GEUV1/Salmon/"
}else{
  output_dir <- "/pine/scr/s/k/skvanbur/GEUV1/SalmonBootSamps/"
}

if(!dir.exists(output_dir)){
  dir.create(output_dir)
}
setwd(output_dir)
E_GEUV_1_sdrf <- read_delim(paste0(dir1, "E-GEUV-1.sdrf.txt"), 
                            "\t", escape_double = FALSE, trim_ws = TRUE)



#This is not the column ENA_SAMPLE, which are the true sample names, but this is the unique identifier used to download the data
#Sort it such that the samples are run in order based on this id to make any troubleshooting from problems with Salmon easier
samps <- sort(unique(E_GEUV_1_sdrf$`Comment[ENA_RUN]`))

curr_samp <- samps[val]


#logf <- paste0(dir1, "Salmon/quantlog", curr_samp, ".txt")
logf <- paste0(output_dir, curr_samp, "/", "quantlog", curr_samp, ".txt")
file1 <- paste0(base_dir, curr_samp, "/", curr_samp, "_1.fastq.gz")
file2 <- paste0(base_dir, curr_samp, "/", curr_samp, "_2.fastq.gz")

outfiles <- paste0(output_dir, curr_samp)



if(!dir.exists(outfiles)){
  dir.create(outfiles) 
}



index <- paste0(output_dir, "transcripts_index")

if(!dir.exists(index) & GibbsSamps==FALSE){
  stop("Ensure the Index exists in the proper location before continuing")
}

# if(!file.exists("/pine/scr/s/k/skvanbur/GEUV1/SalmonBootSamps/transcripts_index/hash.bin")){
#   cpcmd <- "cp -a /pine/scr/s/k/skvanbur/GEUV1/Salmon/transcripts_index/. /pine/scr/s/k/skvanbur/GEUV1/SalmonBootSamps/transcripts_index/"
#   system(cpcmd)
# }

# p1 <- "sbatch -N 1 -n 4 --mem=32g --time=24:00:00 -o "
# p2 <- " --wrap ' ~/bin/Salmon0.8.2/bin/salmon quant -p 4 -i "
# p3 <-  " -l A -1 "
# p4 <- " -2 " 
# p5 <- " -o "
# p6 <- "  --seqBias --gcBias --numGibbsSamples 100'"
# 
# cmd <- paste0(p1, logf, p2, index, p3, file1, p4, file2, p5, outfiles, p6)




p2 <- "~/bin/Salmon0.11.3/bin/salmon quant -p 1 -i "
#Add --dumpEq option to output equivalence class information
#p3 <-  "-l A -1 "
p3 <-  " --dumpEq -l A -1 "
p4 <- " -2 " 
p5 <- " -o "

if(GibbsSamps==TRUE){
  p6 <- paste0(" --seqBias --gcBias --numGibbsSamples ", ninf, " --thinningFactor ", thinf)
}else{
  if(reRunToGenerateEqClassFiles==TRUE){
    #Don't regenerate the bootstrap samples to save time
      #Existing bootstrap samples from a previous run will still remain
    p6 <- paste0(" --seqBias --gcBias")
  }else{
    p6 <- paste0(" --seqBias --gcBias --numBootstraps ", ninf)
  }
  #Bootstrap samples have been generated
}


if(reRunToGenerateEqClassFiles==TRUE){
  if(!file.exists(paste0(outfiles, "/aux_info/eq_classes.txt"))){
    cmd <- paste0(p2, index, p3, file1, p4, file2, p5, outfiles, p6)
    system(cmd)
  }else{
    #Don't need to rerun if the file already exists
  }
}else{
  cmd <- paste0(p2, index, p3, file1, p4, file2, p5, outfiles, p6)
  system(cmd)
}



print(gc())
# sbatch -N 1 -n 4 --mem=32g --time=24:00:00 -o ~/res/GEUV1/Salmon/quant1log.txt --wrap "~/bin/Salmon0.8.2/bin/salmon quant -p 4 -i /pine/scr/s/k/skvanbur/GEUV1/Salmon/transcripts_index -l A -1 /pine/scr/s/k/skvanbur/GEUV1/Data/ERR188465/ERR188465_1.fastq.gz -2 /pine/scr/s/k/skvanbur/GEUV1/Data/ERR188465/ERR188465_2.fastq.gz  -o ~/res/GEUV1/Salmon/transcripts_quantRep1  --seqBias --gcBias --numGibbsSamples 100"
