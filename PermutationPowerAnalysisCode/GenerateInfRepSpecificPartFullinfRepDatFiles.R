array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if(version$nickname=="Planting of a Tree"){
  .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
}else{
  if(onCluster==TRUE){source("~/.Rprofile")}
}

curr_part_num <- array_val
func_loc <-  "~/code/CompFunctions.R"
source(func_loc)

GibbsSamps <- TRUE

if(GibbsSamps==TRUE){
  load_dir <- "/pine/scr/s/k/skvanbur/GEUV1/GibbsFullinfRepDat/"
}else{
  load_dir <- "/pine/scr/s/k/skvanbur/GEUV1/BootFullinfRepDat/"
}

nparts <- 100
ninfreps <- 100


SaveInfRepSpecificPartFullinfRepDatFiles(curr_part_num = curr_part_num, load_dir = load_dir, nparts = nparts, ninfreps = ninfreps)

print(gc())