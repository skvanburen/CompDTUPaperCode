#CorrectLowExpression is not needed for these analyses because that is done on the proportion scale before conversion to ilr scale 
  #and the data here is generated on the ilr scale directly

#array_val <- 1
if(version$nickname=="Planting of a Tree"){
  .libPaths("/nas/longleaf/home/skvanbur/bin/R3.6.0")
}else{
  source("~/.Rprofile")
}
library(gtools)
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

#Set directories
def_wd <- "~/res/GEUV1Data/RSimulationCode/"
save_dir <- "/pine/scr/s/k/skvanbur/GEUV1/"
setwd(def_wd)

#Set to true if you are running results for the NSamp values of 26, 50, or 80 (ie the additional values) and false if you are
  #running them for the original values of 10, 100, or 250
#Specify only one of RunAdditionalNSamps of RunAdditionalNSamps equal to TRUE
RunAdditionalNSamps <- FALSE

#Set to true if you are running results for the effect size values 1.375 and 1.75
  #This will run the results for these effect sizes for all sample sizes, including the additional ones
  #ie will run for nsamp of 10, 100, 250, and 26 and 50
#Specify at most one of RunAdditionalNSamps of RunAdditionalNSamps equal to TRUE
RunAdditionalEffectSizes <- TRUE

#Split each simulation up into 10 parts to speed up computation
  #ie for each simulation scenario (each row of RSimVals) run 10 jobs with 1000 runs each to give a total of
  #10000 repetitions for each simulation scenario
  #Since a slurm array specification value can't be larger than 40000 need to split running based on the values of NSamp
  #Specifically, if RunAdditionalNSamps is TRUE the results run will be for NSamp 26, 50, 80 and if FALSE will be 10,100,250
if(RunAdditionalNSamps==TRUE){
  sim_number <- ceiling(array_val/100) + 357
  part_of_sim_number <- array_val%%100
}else if(RunAdditionalEffectSizes==TRUE){
  sim_number <- ceiling(array_val/100) + 714
  part_of_sim_number <- array_val%%100
}else{
  sim_number <- ceiling(array_val/100)
  part_of_sim_number <- array_val%%100
}


Sys.sleep(sample(1:50, size = 1))
load("RSimVals.RData")


func_loc <- "RSimulationCodeFunctions.R"
source(func_loc)


curr_vals <- RSimVals[sim_number,]

fil <- "FilesToRunRSimulationCode3ResponsesDRIMSeqFiltered.RData"


nsamp <- curr_vals$nsamp
if(nsamp==250){
  stop("Don't run results for nsamp values of 250 for now to save computation time")
}
ncondlevels <- curr_vals$ncondlevels
GenWithMeasError <- curr_vals$GenWithMeasError
MeasErrorMultFactor <- curr_vals$MeasErrorMultFactor
BetweenCovMultFactor <- curr_vals$BetweenCovMultFactor
condEffSizeStr <- curr_vals$condEffSizeStr

startseed_sim <- curr_vals$StartSeed
#nreps_to_use_list <- c(1,10,50,100,250)
nreps_to_use_list <- c(1,10,50,100,250)

WithinSubjCovSampType <- "Permute"
Wishdf <- NA

#Nruns is set to 1000 here to get 10000 total runs per simulation scenario (10 parts * 1000 runs)
nruns <- 1000

load(fil)


#seedstouse <- seq(startseed, startseed + nruns - 1, 1)
startseed <- startseed_sim + ((part_of_sim_number-1) * nruns)
seedstouse <- seq(startseed, startseed + nruns - 1, 1)
startTime <- proc.time()
assign(paste0("SimRes", sim_number, "Part", part_of_sim_number), apply(as.matrix(1:nruns, ncol = 1), 1, RSimulationCode, fil = fil, nsamp = nsamp, 
                                          nreps_to_use_list = nreps_to_use_list, WithinSubjCovSampType = WithinSubjCovSampType, 
                                          Wishdf = Wishdf, seedstouse = seedstouse, ncondlevels = ncondlevels,  
                                          condEffSizeStr = condEffSizeStr, MeasErrorMultFactor = MeasErrorMultFactor, 
                                          BetweenCovMultFactor = BetweenCovMultFactor, GenWithMeasError = GenWithMeasError,
                                          useRobustCovEst = FALSE, GibbsCovsToUse = GibbsCovsToUse, res1 = res1, res1null = res1null, SigmaTildeAlt = SigmaTildeAlt,
                                          SigmaTildeNull = SigmaTildeNull, Y = Y, cond = cond, gene_id = gene_id))

curr_SimRes <- get(paste0("SimRes", sim_number, "Part", part_of_sim_number))
assign(paste0("SimResFinal", sim_number, "Part", part_of_sim_number), data.frame(data.table::rbindlist(curr_SimRes)))
comptime <- proc.time() - startTime
print(comptime)

#Code to save results on cluster
direc <- paste0(save_dir, "RSimulationCodeRes/", "SimNum", sim_number, "/")

if(!dir.exists(direc)){
  dir.create(direc, recursive = TRUE)
}

restosave <- paste0("SimResFinal", sim_number, "Part", part_of_sim_number)
print(restosave)

#Sys.sleep(sample(1:500, size = 1))
save(list = c(restosave), file = paste0(direc, "SimResFinal", sim_number, "Part", part_of_sim_number, ".RData"))

print(gc())




#This code will help find the ones that did not successfully run so you can resubmit them

if(!is.numeric(array_val)){
  load(paste0("~/res/GEUV1Data/RSimulationCode/", "RSimVals.RData"))
  array_vals_to_resubmit <- c()
  nparts_run_per_sim <- 10
  #j in 1:357
  for(j in 1:204){
    # if(j %in% seq(3, 357, 3)){
    #   next
    # }

    #sim_number <- j
    #sim_number <- j + 357
    sim_number <- j + 714
    if(RSimVals$nsamp[sim_number]==250){next}
    direc_to_look <- paste0("/pine/scr/s/k/skvanbur/GEUV1/RSimulationCodeRes/SimNum", sim_number, "/")
    for(h in 1:nparts_run_per_sim){
      part_of_sim_number <- h
      curr_fil <- paste0(direc_to_look, "SimResFinal", sim_number, "Part", part_of_sim_number, ".RData")

      if(!file.exists(curr_fil)){
        print(paste0("File does not exist for Simulation Number ", sim_number, " and Part ", part_of_sim_number))
        a_val <- ((j - 1) * 100) + part_of_sim_number
        array_vals_to_resubmit <- c(array_vals_to_resubmit, a_val)
      }
    }
  }

  for(jj in 1:length(array_vals_to_resubmit)){
    system(paste0("cat /pine/scr/s/k/skvanbur/GEUV1/RSimsLogs/RSimulationCode_", array_vals_to_resubmit[jj], ".Rout"))
    print(array_vals_to_resubmit[jj])
  }
  # print(paste0("Need to resubmit array values:"))
  # print(array_vals_to_resubmit)


  for(p in 1:length(array_vals_to_resubmit)){
    curr_v <- array_vals_to_resubmit[p]
    #system(paste0("cd ~/res/GEUV1Data/RSimulationCode/"))
    #system(paste0("module load r/3.6.0"))
    system(paste0("module load r/3.6.0; sbatch --array=", curr_v, " ~/res/GEUV1Data/RSimulationCode/RSimulationCode.sh"))
  }


}
