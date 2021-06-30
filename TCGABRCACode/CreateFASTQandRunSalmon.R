.libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
#The file for array value 359 (/pine/scr/s/k/skvanbur/TCGABRCAAnalysis/Data/6d135cc4-4eca-48be-8f87-bf94daf83e7a/108dc6d1-5612-4276-a682-0ba4d324cb00_gdc_realn_rehead.bam)
  #is clearly a corrupt/bad BAM file
  #I tried redownloading BAM file several times and the conversion to fastq always failed with the error below
  #All other samples worked without any problems
#So, we will have to skip this sample
# [E::bgzf_uncompress] Inflate operation failed: progress temporarily not possible, or in() / out() returned an error
# [E::bgzf_read] Read block operation failed with error 1 after 96 of 142 bytes
# [bam2fq_mainloop] Failed to read bam record.
# [bam2fq_mainloop] Error writing to FASTx files.: No such file or directory
# [M::bam2fq_mainloop] discarded 0 singletons
# [M::bam2fq_mainloop] processed 102641709 reads
# samtools bam2fq: error closing "/pine/scr/s/k/skvanbur/TCGABRCAAnalysis/Data/6d135cc4-4eca-48be-8f87-bf94daf83e7a/108dc6d1-5612-4276-a682-0ba4d324cb00_gdc_realn_rehead.bam": -1

#6cb6f179-defd-4661-af0a-c353b74c0c49
#8785012f-f73e-4d68-87cf-1d804af32782
#f130f376-5801-40f9-975d-a7e2f7b5670d
library(readr)
base_dir <- "/pine/scr/s/k/skvanbur/TCGABRCAAnalysis/"

array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

file_manifest <- data.frame(read_delim("~/res/TCGABRCAAnalysis/FilesForGDCDownloader/final_manifest_all.txt", delim = "\t"))

VerifyCommandCompleted <- function(cmd){
  #Old code that has command specified as a string
  # command_val <- get(cmd)
  # if(command_val!=0){
  #   stop(paste0("The command ", cmd, " failed to complete successfully, giving exit code ", command_val))
  # }else{
  #   print(paste0("The command ", cmd, " completed successfully, giving exit code ", command_val))
  # }
  cmd_name <- deparse(substitute(cmd))
  if(cmd==0){
    print(paste0("The command ", cmd_name, " completed successfully, giving exit code ", "0"))
  }else{
    stop(paste0("The command ", cmd_name, " failed to complete successfully, giving exit code ", cmd))
  }
}

curr_r <- file_manifest[array_val, , drop = F]
curr_barcode <- curr_r[1, "barcode"]
print(paste0("Current barcode is ", curr_barcode))

#This id is not to be confused with the id below, which is the file id (or filename)
curr_dir <- paste0(base_dir, "Data/", curr_r[1, "id"], "/")
print(paste0("Current directory is ", curr_dir))

if(!dir.exists(curr_dir)){
  print("The current directory doesn't exist, meaning the data has not yet been downloaded for this sample so nothing can be run yet")
}else{
  curr_file_id <- curr_r[1, "filename"]
  curr_bam_file <- paste0(curr_dir, curr_file_id)
  curr_bam_bai_file <- paste0(strsplit(curr_bam_file, split = ".bam")[[1]][1], ".bai")
  
  #Shouldn't be any reads that designate "other" but pass this option just to confirm
  curr_fastq_0 <- paste0(curr_dir, "r0.fastq.gz")
  
  #Specify output names of the fastq files
  curr_fastq_1 <- paste0(curr_dir, "r1.fastq.gz")
  curr_fastq_2 <- paste0(curr_dir, "r2.fastq.gz")
  
  #Check for existence of the bam bai file because this is downloaded after the bam file
    #This ensures the process doesn't try to start while the bam file is being downloaded
    #Also, can't just check for existence of the fastq files because they are created at the beginning of the process
    #and added to so they could still exist even if the process failed
  if(file.exists(curr_bam_bai_file)){
    bam2fqcmd <- paste0("module load samtools;", "samtools bam2fq -N ", curr_bam_file, " -0 ", curr_fastq_0, " -1 ", curr_fastq_1, " -2 ", curr_fastq_2)
    print(paste0("bam2fq command to run is ", bam2fqcmd))
    
    bam2fq <- system(bam2fqcmd)
    VerifyCommandCompleted(bam2fq)
    
    #If the above code worked successfully, remove the current bam and bai files
    system(paste0("rm -f ", curr_bam_file))
    system(paste0("rm -f ", curr_bam_bai_file))
  }

  #Specify output names of the sorted fastq files
  sorted_fastq_1 <- paste0(curr_dir, "r1sorted.fastq.gz")
  sorted_fastq_2 <- paste0(curr_dir, "r2sorted.fastq.gz")
  
  #Now, these fastq's upon inspection didn't appear to be guaranteed to be in the same order-  sort them
    #This very helpful code was written by biostars user dariober, see https://www.biostars.org/p/15011/#103041 
    #though it may have come from this post first https://edwards.sdsu.edu/research/sorting-fastq-files-by-their-sequence-identifiers/
  
  sort_cmd1 <- paste0("zcat ", curr_fastq_1, " | paste - - - - | sort -V -k1,1 -S 3G | tr '\t' '\n' | gzip > ", sorted_fastq_1)
  sort_cmd2 <- paste0("zcat ", curr_fastq_2, " | paste - - - - | sort -V -k1,1 -S 3G | tr '\t' '\n' | gzip > ", sorted_fastq_2)
  print(paste0("sort command to run for r1 is ", sort_cmd1))
  print(paste0("sort command to run for r2 is ", sort_cmd2))
  
  if(file.exists(curr_fastq_1)){
    s_cmd1 <- system(sort_cmd1)
    VerifyCommandCompleted(s_cmd1)
    
    system(paste0("rm -f ", curr_fastq_1))
  }
  
  
  if(file.exists(curr_fastq_2)){
    s_cmd2 <- system(sort_cmd2)
    VerifyCommandCompleted(s_cmd2)
    
    system(paste0("rm -f ", curr_fastq_2))
  }
  
  #Now, once the sorted fastq files have been generated, run Salmon
  #Can reuse the same index that was used before - corresponding to GENCODE v27
  #See the GEUV1 Code for code to create this if it doesn't exist
  index <- paste0("/pine/scr/s/k/skvanbur/GEUV1/SalmonBootSamps/", "transcripts_index")
  
  salmon_output_dir <- paste0(base_dir, "SalmonBootSamps/", curr_barcode, "/")
  if(!dir.exists(salmon_output_dir)){dir.create(salmon_output_dir, recursive = T)}
  
  #If this file doesn't exist, Salmon has not been run or did not complete properly so repeat it
  if(!file.exists(paste0(salmon_output_dir, "aux_info/", "meta_info.json"))){
    ninf <- 100
    p2 <- "~/bin/Salmon0.11.3/bin/salmon quant -p 1 -i "
    #Add --dumpEq option to output equivalence class information
    #p3 <-  "-l A -1 "
    p3 <-  " --dumpEq -l A -1 "
    p4 <- " -2 " 
    p5 <- " -o "
    p6 <- paste0(" --seqBias --gcBias --numBootstraps ", ninf)
    
    
    cmd_to_run <- paste0(p2, index, p3, sorted_fastq_1, p4, sorted_fastq_2, p5, salmon_output_dir, p6)
    salmon_cmd <- system(cmd_to_run)
    VerifyCommandCompleted(salmon_cmd)
  }
  
}
