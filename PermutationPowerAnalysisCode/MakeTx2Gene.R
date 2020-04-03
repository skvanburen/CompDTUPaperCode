#Need to have run "SaveSalmonDataAsRData.R" for this to work
#Code to create the tx2gene file that matches transcripts to genes

onCluster <- FALSE

if(version$nickname=="Planting of a Tree"){
  .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
}else{
  if(onCluster==TRUE){source("~/.Rprofile")}
}

def_wd_clust <- "~/res/GEUV1Data/"
def_wd_local <- "/Users/Scott/Documents/Dissertation Data/GEUV1Data/"
func_loc_clust <- "~/code/CompFunctions.R"
func_loc_local <- "/Users/Scott/Documents/Dissertation/Code/CompFunctions.R"
txdb_loc_local <- "/Users/Scott/Documents/Dissertation Data/Gencodev27/gencode.v27.annotation.gtf.gz"
txdb_loc_clust <- "~/gencode.v27.annotation.gtf.gz"

if(onCluster==TRUE){
  setwd(def_wd_clust)
  func_loc <- func_loc_clust
  txdb_loc <- txdb_loc_clust
}else{
  setwd(def_wd_local)
  func_loc <- func_loc_local
  txdb_loc <- txdb_loc_local
}

source(func_loc)
#Build dataframe matching transcripts to genes using GENCODE if it doesn't exist
maketx2gene(txdb_loc = txdb_loc)
