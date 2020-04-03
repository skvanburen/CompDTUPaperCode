
#Once this code is done for a specific file, need to manually check the DownloadGEUV1Dataval.log file for each x and
  #make sure the download actually worked- sometimes the server or cluster have weird errors/failures that are
  #fixed by just rerunning the file for that value

#When trying to generate the eq class files to use with BANDITS, it became apparent that some of the FASTQ files were broken, which resulted in the eq class files
  #Not generating.  The data had been previously redownloaded after the full set of results were run, so those results were unaffected
  #This code will identify the broken FASTQ files (based on samples where the Salmon run failed to generate an eq class) and delete them
  #Then, following this, simply submit the array job for all values with OverrideFiles set to FALSE such that it only downloads
  #files thata were broken (and have now since been deleted)
# load("~/res/GEUV1Data/cntGenecntsScaledTPMFiltered.RData")
# base_dir <- "/pine/scr/s/k/skvanbur/GEUV1/SalmonBootSamps/"
# dir_modifiers <- key$ENARun
# samps_to_keep <- 1:462
# 
# broken_files <- c()
# broken_dir_modifiers <- c()
# 
# for(di in 1:462){
#   print(paste0("current number is ", di))
#   curr_eq_class_fil <- paste0(base_dir, dir_modifiers[di], "/", "aux_info/", "eq_classes.txt")
#   if(!file.exists(curr_eq_class_fil)){
#     broken_files <- c(broken_files, di)
#     broken_dir_modifiers <- c(broken_dir_modifiers, dir_modifiers[di])
#   }
# }

# for(jj in 1:length(broken_dir_modifiers)){
#   curr_direc_to_remove <- paste0("/pine/scr/s/k/skvanbur/GEUV1/Data/", broken_dir_modifiers[jj], "/")
#   system(paste0("rm -r ", curr_direc_to_remove))
# }

#Val should run from 1 to 924 in this file, since there 924 rows in the .txt file, which each row containing a file that needs to be downloaded
val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

#Put in a sleep statement just to prevent too many jobs from trying to access files at the same time
#Sys.sleep(sample(1:100, 1))
onCluster <- TRUE
if(version$nickname=="Planting of a Tree"){
  .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
}else{
  if(onCluster==TRUE){source("~/.Rprofile")}
}

library(readr)

#Set OverrideFiles to be True if you want to override files that are already there
  #ie if you want to replace corrupted files, etc
OverrideFiles <- FALSE
dir1 <- "~/res/GEUV1Data/"
base_dir <- "/pine/scr/s/k/skvanbur/GEUV1/Data/"

E_GEUV_1_sdrf <- read_delim(paste0(dir1, "E-GEUV-1.sdrf.txt"), 
                            "\t", escape_double = FALSE, trim_ws = TRUE)




cols <- colnames(E_GEUV_1_sdrf) %in% c("Comment[ENA_SAMPLE]", "Characteristics[population]", "Comment[ENA_RUN]")
sub1 <- E_GEUV_1_sdrf[,cols]

#This is not the column ENA_SAMPLE, which are the true sample names, but this is the unique identifier used to download the data so I use this as the
  #"sample name"
  names <- E_GEUV_1_sdrf$`Comment[ENA_RUN]`


  curr <- names[val]
  if(!dir.exists(paste0(base_dir, curr))){
    dir.create(paste0(base_dir, curr)) 
  }

  direc <- paste0(base_dir, as.character(E_GEUV_1_sdrf[val, "Comment[ENA_RUN]"]))
  setwd(direc)
  
  if(val%%2==1){
    curr_file <- paste0(as.character(E_GEUV_1_sdrf[val, "Comment[ENA_RUN]"]), "_1.fastq.gz")
  }else{
    curr_file <- paste0(as.character(E_GEUV_1_sdrf[val, "Comment[ENA_RUN]"]), "_2.fastq.gz")
  }
  
  if(!file.exists(curr_file) | OverrideFiles==TRUE){
    curr <- as.character(E_GEUV_1_sdrf[val,"Comment[FASTQ_URI]"])
    print(paste0("Current download link is ", curr))
    p1 <- paste0("wget ", curr, "")
    cmd <- paste("cd", direc, ";", p1)
    
    #Put in a sleep statement just to prevent too many jobs from trying to access the download at the same time
    #Sys.sleep(sample(1:10000, 1))
    
    system(cmd)
  }
  
  
  #Code to check which files exist
  # files_missing <- c()
  # for(val in 1:924){
  #   curr <- names[val]
  #   direc <- paste0(base_dir, as.character(E_GEUV_1_sdrf[val, "Comment[ENA_RUN]"]))
  #   if(!dir.exists(direc)){
  #     files_missing <- c(files_missing, val)
  #     next
  #   }
  #   setwd(direc)
  #   if(val%%2==1){
  #     curr_file <- paste0(as.character(E_GEUV_1_sdrf[val, "Comment[ENA_RUN]"]), "_1.fastq.gz")
  #   }else{
  #     curr_file <- paste0(as.character(E_GEUV_1_sdrf[val, "Comment[ENA_RUN]"]), "_2.fastq.gz")
  #   }
  # 
  # 
  # 
  #   if(!file.exists(curr_file)){
  #     files_missing <- c(files_missing, val)
  #   }
  # }
  
  