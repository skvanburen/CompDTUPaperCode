#Power Table code March 2020

library(xtable)
library(gtools)
library(tikzDevice)
library(data.table)

save_dir <- "/Users/Scott/Documents/Dissertation/Paper1/Tables"


#Neat little function that will load an .RData object and save its contents to whatever as the results of the output
#For example can do d <- loadRData("file.RData") and whatever is loaded is saved as d
loadRData <- function(fileName, objNameToGet = NULL){
  #loads an RData file, and returns it
  load(fileName)
  #print(ls()[ls() != "fileName"])
  if(is.null(objNameToGet)){
    rm(objNameToGet)
    #print(ls()[ls() != "fileName"])
    return(get(ls()[ls() != "fileName"]))
  }else{
    return(get(objNameToGet))
  }
  
}

fr <- function(x, n = getOption("digits"), InsertLessThanFormatting = TRUE) {
  if("data.table" %in% class(x)){
    leny <- nrow(x)
  }else{
    leny <- length(x)
  }
  y <- rep(0, leny)
  for (i in 1:leny) {
    curr_val <- x[i]
    if(!is.numeric(curr_val)){
      curr_val <- as.numeric(curr_val)
    }
    if (curr_val >=0 & curr_val< 10^-(n) & !is.na(curr_val) & InsertLessThanFormatting==TRUE){
      y[i] <- paste("$<$", ".", paste(rep("0", n-1), collapse=""), "1", sep = "")
    }
    if (is.numeric(curr_val) & !is.na(curr_val) &  (curr_val >= 10^-(n) | curr_val < 0 | InsertLessThanFormatting==FALSE)) {
      y[i] <- format(round(x[i], digits = n), nsmall = n, scientific = FALSE)
    }
    # Only truly missing values should ever show up as 0 in y created above, so change them back to missing
    #Note that this is not the same as having a value of 0 in the input vector (x), as these will correctly
    #Just be formatted to show as <0.0001, etc
    if(y[i] == "0") {
      y[i] <- " "
    }
  }
  return(y)
}
#setwd("/Users/Scott/Documents/Dissertation Data/GEUV1Data/TestNewGibbsApproachMVNRes1ResponseDRIMSeqFiltered")
#
#setwd("/Users/Scott/Documents/Dissertation Data/GEUV1Data/TestNewGibbsApproachMVNRes2ResponsesDRIMSeqFiltered")
#setwd("/Users/Scott/Documents/Dissertation Data/GEUV1Data/TestNewGibbsApproachMVNRes3ResponsesDRIMSeqFiltered")

#load("/Users/Scott/Documents/Dissertation/res/GEUV1Data/TestNewGibbsApproachMVNsimvals.RData")
#setwd("/Users/Scott/Documents/Dissertation Data/GEUV1Data/TestNewGibbsApproachMVNRes")


load("/Users/Scott/Documents/Dissertation/Paper1/SimulationCode/RSimulationCode/RSimVals.RData")

save_dir <- "/Users/Scott/Documents/Dissertation Data/GEUV1Data/RSimulationCodeRes/"
setwd(save_dir)

h <- nrow(RSimVals)
#h <- 357

ResultsSavedinParts <- TRUE

if(ResultsSavedinParts==TRUE){
  npartsrun <- 10
  for(i in 1:h){
    print(paste0("i is ", i, " out of ", h))
    curr_sim <- i
    res <- vector(mode = "list", length = npartsrun)
    
    if(!dir.exists(paste0(save_dir, "SimNum", curr_sim, "/"))){
      next
    }
    
    for(j in 1:npartsrun){
      res[[j]] <- loadRData(paste0(save_dir, "SimNum", curr_sim, "/", "SimResFinal", curr_sim, "Part", j, ".RData"), objNameToGet = paste0("SimResFinal", curr_sim, "Part", j))
    }
    
    res2 <- data.frame(data.table::rbindlist(res))
    assign(paste0("SimResFinal", curr_sim), res2)
    
    save(list = paste0("SimResFinal", curr_sim), file = paste0(save_dir, "SimRes", curr_sim, ".RData"))
  }
}


for(i in 1:h){
  if(file.exists(paste0("SimRes", i, ".RData"))){
    load(file = paste0("SimRes", i, ".RData"))
  }else{
    next
  }
  
}

for (i in 1:h){
  print(paste0("i is ", i, " out of ", h))
  if(file.exists(paste0("SimRes", i, ".RData"))){
    curr_power <- data.frame("A")
    curr_res <- get(paste0("SimResFinal", i))
    curr_cols <- colnames(curr_res)
    for(j in 1:length(curr_cols)){
      curr_p <- mean(curr_res[,curr_cols[j]] < 0.05, na.rm = TRUE)
      power_col <- paste0("power_", curr_cols[j])
      curr_power[,power_col] <- curr_p
    }
    curr_power$X.A. <- NULL
    assign(paste0("Power", i), curr_power)
  }else{
    #If the file does not exist, assign the power to -1 such that you know the results don't exist
    #Load SimResFinal1 just to get the column names of the methods
    curr_power <- data.frame("A")
    curr_res <- get(paste0("SimResFinal", 1))
    curr_cols <- colnames(curr_res)
    for(j in 1:length(curr_cols)){
      curr_p <- -1
      power_col <- paste0("power_", curr_cols[j])
      curr_power[,power_col] <- curr_p
    }
    curr_power$X.A. <- NULL
    assign(paste0("Power", i), curr_power)
  }
  
}

#Get the sim number for current parameter values
gPV <- function(n,cESS, G, H, B, method, nreps = NA){
  UniqString <- paste0("nsamp", n, "condEffSizeStr", cESS, 
                       "GenWithMeasError", G, "MeasErrorMultFactor", H, 
                       "BetweenCovMultFactor", B)
  curr_S <- as.numeric(rownames(RSimVals)[RSimVals$UniqString==UniqString])
  if(length(curr_S)==1){
    curr_SimNum <- curr_S
  }else if(length(curr_S)>1){
    print(paste0("The string ", UniqString, " was matched to multiple rows in the table.  This means the same experiment was repeated multiple times."))
    print(paste0("The first value will be returned as the correct values to take.  Ensure this is accurate.  The full list of values that match this string are given below"))
    print(curr_S)
    curr_SimNum <- curr_S[1]
  }else if(length(curr_S)==0){
    stop("No simulation was found corresponding to the values input.  Check that the inputs have been specified correctly, that all files have been downloaded, or that all results have been run")
  }else{
    stop("There was an error when extracting the simulation number corresponding to the current values")
  }
  
  #print(paste0("Current Sim is ", curr_SimNum))
  CurrentResults <- get(paste0("Power", curr_SimNum))
  
  if(method!="CompDTU" & is.null(nreps)){
    stop("Nreps must be specified for methods other than CompDTU")
  }
  
  if(method=="CompDTU"){
    pval_to_use <- "power_pval_CompDTU"
  }else if(method=="CompDTUme"){
    pval_to_use <- paste0("power_pval_CompDTUme", nreps, "Rep")
  }else if(method=="CompMICombinePvals"){
    pval_to_use <- paste0("power_pval_CompMICombinePvals", nreps, "Rep")
  }else if(method=="CompMICombineCoefs"){
    pval_to_use <- paste0("power_pval_CompMICombineCoefs", nreps, "Rep")
  }else{
    stop("Incorrect method specification")
  }
  
  return(CurrentResults[,pval_to_use])
}


PowerTableFn <- function(nsampcurr, nreps_curr, curr_methd, curr_H, curr_B, IncBetVarTable = FALSE){
  if(IncBetVarTable==TRUE){
    powerT <- data.table("NSamp" = nsampcurr, "NReps" = nreps_curr,  "Method" = curr_methd, "B" = curr_B, "1 (Null)" = gPV(nsampcurr, "c(1,1)", TRUE, curr_H, curr_B, curr_methd, nreps_curr), "1.10" = gPV(nsampcurr, "c(1,1.10)", TRUE, curr_H, curr_B, curr_methd, nreps_curr), "1.25" = gPV(nsampcurr, "c(1,1.25)", TRUE, curr_H, curr_B, curr_methd, nreps_curr), "1.375" = gPV(nsampcurr, "c(1,1.375)", TRUE, curr_H, curr_B, curr_methd, nreps_curr), "1.50" = gPV(nsampcurr, "c(1,1.50)", TRUE, curr_H, curr_B, curr_methd, nreps_curr), "1.75" = gPV(nsampcurr, "c(1,1.75)", TRUE, curr_H, curr_B, curr_methd, nreps_curr), "2.00" = gPV(nsampcurr, "c(1,2.00)", TRUE, curr_H, curr_B, curr_methd, nreps_curr), stringsAsFactors = FALSE)
  }else{
    powerT <- data.table("NSamp" = nsampcurr, "NReps" = nreps_curr,  "Method" = curr_methd, "H" = curr_H, "1 (Null)" = gPV(nsampcurr, "c(1,1)", TRUE, curr_H, curr_B, curr_methd, nreps_curr), "1.10" = gPV(nsampcurr, "c(1,1.10)", TRUE, curr_H, curr_B, curr_methd, nreps_curr), "1.25" = gPV(nsampcurr, "c(1,1.25)", TRUE, curr_H, curr_B, curr_methd, nreps_curr), "1.375" = gPV(nsampcurr, "c(1,1.375)", TRUE, curr_H, curr_B, curr_methd, nreps_curr), "1.50" = gPV(nsampcurr, "c(1,1.50)", TRUE, curr_H, curr_B, curr_methd, nreps_curr), "1.75" = gPV(nsampcurr, "c(1,1.75)", TRUE, curr_H, curr_B, curr_methd, nreps_curr), "2.00" = gPV(nsampcurr, "c(1,2.00)", TRUE, curr_H, curr_B, curr_methd, nreps_curr), stringsAsFactors = FALSE)
  }
  return(powerT)
}

InsertBlankLine <- function(table, IncBetVarTable = FALSE){
  if(IncBetVarTable==TRUE){
    return(data.table("NSamp" = "", "NReps" = "",  "Method" = "", "B" = "", "1 (Null)" = "", "1.10" = "", "1.25" = "", "1.375" = "", "1.50" = "", "1.75" = "", "2.00" = "", stringsAsFactors = FALSE))
  }else{
    return(data.table("NSamp" = "", "NReps" = "",  "Method" = "", "H" = "", "1 (Null)" = "", "1.10" = "", "1.25" = "", "1.375" = "", "1.50" = "", "1.75" = "", "2.00" = "", stringsAsFactors = FALSE))
  }
}


#Main text power table
nsampcurr <- 10
curr_B <- 1
curr_H <- 0.01
MainTextPowerTableT <- PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B)
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, InsertBlankLine(MainTextPowerTableT), fill = TRUE)

curr_H <- 1
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, InsertBlankLine(MainTextPowerTableT), fill = TRUE)

curr_H <- 2
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, InsertBlankLine(MainTextPowerTableT), fill = TRUE)
MainTextPowerTableT <- rbind(MainTextPowerTableT, InsertBlankLine(MainTextPowerTableT), fill = TRUE)

# curr_H <- 4
# MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
# MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
# MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
# MainTextPowerTableT <- rbind(MainTextPowerTableT, InsertBlankLine(MainTextPowerTableT), fill = TRUE)


nsampcurr <- 26

curr_H <- 0.01
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, InsertBlankLine(MainTextPowerTableT), fill = TRUE)

curr_H <- 1
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, InsertBlankLine(MainTextPowerTableT), fill = TRUE)

curr_H <- 2
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, InsertBlankLine(MainTextPowerTableT), fill = TRUE)
MainTextPowerTableT <- rbind(MainTextPowerTableT, InsertBlankLine(MainTextPowerTableT), fill = TRUE)

# curr_H <- 4
# MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
# MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
# MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
# MainTextPowerTableT <- rbind(MainTextPowerTableT, InsertBlankLine(MainTextPowerTableT), fill = TRUE)



nsampcurr <- 100
curr_H <- 0.01
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, InsertBlankLine(MainTextPowerTableT), fill = TRUE)

curr_H <- 1
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, InsertBlankLine(MainTextPowerTableT), fill = TRUE)

curr_H <- 2
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
#MainTextPowerTableT <- rbind(MainTextPowerTableT, InsertBlankLine(MainTextPowerTableT), fill = TRUE)

# curr_H <- 4
# MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
# MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
# MainTextPowerTableT <- rbind(MainTextPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
#MainTextPowerTableT <- rbind(MainTextPowerTableT, InsertBlankLine(MainTextPowerTableT), fill = TRUE)

RSimsPowerTableUpdatedOrderingT <- data.frame(MainTextPowerTableT)
RSimsPowerTableUpdatedOrderingT$H <- fr(as.numeric(RSimsPowerTableUpdatedOrderingT$H), 2)
for(j in 5:ncol(RSimsPowerTableUpdatedOrderingT)){
  RSimsPowerTableUpdatedOrderingT[,j] <- fr(as.numeric(RSimsPowerTableUpdatedOrderingT[,j]),3)
}
RSimsPowerTableUpdatedOrdering <- data.table(RSimsPowerTableUpdatedOrderingT)
colnames(RSimsPowerTableUpdatedOrdering) <- colnames(MainTextPowerTableT)
save(RSimsPowerTableUpdatedOrdering, file = paste0(save_dir, "/", "RSimsPowerTableUpdatedOrdering.RData"))
rm(curr_H)
rm(curr_B)
rm(nsampcurr)


#Supplementary Table for CompMI methods
curr_H <- 1
curr_B <- 1

nsampcurr <- 10
CompMISuppPowerTableT <- PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompMICombinePvals", curr_H = curr_H, curr_B = curr_B)
CompMISuppPowerTableT <- rbind(CompMISuppPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompMICombinePvals", curr_H = curr_H, curr_B = curr_B))

#nsampcurr <- 26
#CompMISuppPowerTableT <- rbind(CompMISuppPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompMICombinePvals", curr_H = curr_H, curr_B = curr_B))
#CompMISuppPowerTableT <- rbind(CompMISuppPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompMICombinePvals", curr_H = curr_H, curr_B = curr_B))

nsampcurr <- 50
CompMISuppPowerTableT <- rbind(CompMISuppPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompMICombinePvals", curr_H = curr_H, curr_B = curr_B))
CompMISuppPowerTableT <- rbind(CompMISuppPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompMICombinePvals", curr_H = curr_H, curr_B = curr_B))

nsampcurr <- 100
CompMISuppPowerTableT <- rbind(CompMISuppPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompMICombinePvals", curr_H = curr_H, curr_B = curr_B))
CompMISuppPowerTableT <- rbind(CompMISuppPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompMICombinePvals", curr_H = curr_H, curr_B = curr_B))


CompMISuppPowerTableT <- rbind(CompMISuppPowerTableT, InsertBlankLine(CompMISuppPowerTableT), fill = TRUE)


nsampcurr <- 10
CompMISuppPowerTableT <- rbind(CompMISuppPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompMICombineCoefs", curr_H = curr_H, curr_B = curr_B))
CompMISuppPowerTableT <- rbind(CompMISuppPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompMICombineCoefs", curr_H = curr_H, curr_B = curr_B))

#nsampcurr <- 26
#CompMISuppPowerTableT <- rbind(CompMISuppPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompMICombineCoefs", curr_H = curr_H, curr_B = curr_B))
#CompMISuppPowerTableT <- rbind(CompMISuppPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompMICombineCoefs", curr_H = curr_H, curr_B = curr_B))

nsampcurr <- 50
CompMISuppPowerTableT <- rbind(CompMISuppPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompMICombineCoefs", curr_H = curr_H, curr_B = curr_B))
CompMISuppPowerTableT <- rbind(CompMISuppPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompMICombineCoefs", curr_H = curr_H, curr_B = curr_B))

nsampcurr <- 100
CompMISuppPowerTableT <- rbind(CompMISuppPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompMICombineCoefs", curr_H = curr_H, curr_B = curr_B))
CompMISuppPowerTableT <- rbind(CompMISuppPowerTableT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompMICombineCoefs", curr_H = curr_H, curr_B = curr_B))


RSimsPowerTableSuppMIT <- data.frame(CompMISuppPowerTableT)
RSimsPowerTableSuppMIT$H <- fr(as.numeric(RSimsPowerTableSuppMIT$H), 2)
#RSimsPowerTableSuppMIT$H <- NULL
for(j in 5:ncol(RSimsPowerTableSuppMIT)){
  RSimsPowerTableSuppMIT[,j] <- fr(as.numeric(RSimsPowerTableSuppMIT[,j]),3)
}
RSimsPowerTableSupp <- data.table(RSimsPowerTableSuppMIT)
colnames(RSimsPowerTableSupp) <- colnames(CompMISuppPowerTableT)
RSimsPowerTableSupp$H <- NULL
save(RSimsPowerTableSupp, file = paste0(save_dir, "/", "RSimsPowerTableSuppMI.RData"))
rm(curr_H)
rm(curr_B)
rm(nsampcurr)



#First table of additional simulation results for varied measurement error
curr_B <- 1
nsampcurr <- 10

curr_H <- 0
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B)
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT), fill = TRUE)

curr_H <- 0.01
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT), fill = TRUE)

curr_H <- 0.50
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT), fill = TRUE)

curr_H <- 1
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT), fill = TRUE)

curr_H <- 2
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT), fill = TRUE)

curr_H <- 4
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT), fill = TRUE)
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT), fill = TRUE)


# curr_H <- 10
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT), fill = TRUE)



# curr_H <- 20
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT), fill = TRUE)





curr_B <- 1
nsampcurr <- 26

curr_H <- 0
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT), fill = TRUE)

curr_H <- 0.01
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT), fill = TRUE)

curr_H <- 0.50
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT), fill = TRUE)

curr_H <- 1
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT), fill = TRUE)

curr_H <- 2
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT), fill = TRUE)

curr_H <- 4
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
#RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT), fill = TRUE)


# curr_H <- 10
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
#RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT), fill = TRUE)



# curr_H <- 20
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
# #RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT), fill = TRUE)



RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT2 <- data.frame(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT)
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT2$H <- fr(as.numeric(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT2$H), 2, InsertLessThanFormatting = F)
#RSimsPowerTableSuppMIT$H <- NULL
for(j in 5:ncol(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT2)){
  RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT2[,j] <- fr(as.numeric(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT2[,j]),3)
}
RSimsPowerTableUpdatedOrderingSuppIncMeasError <- data.table(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT2)
colnames(RSimsPowerTableUpdatedOrderingSuppIncMeasError) <- colnames(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorT)
save(RSimsPowerTableUpdatedOrderingSuppIncMeasError, file = paste0(save_dir, "/", "RSimsPowerTableUpdatedOrderingSuppIncMeasError.RData"))
rm(curr_H)
rm(curr_B)
rm(nsampcurr)



#Second table of additional simulation results for varied measurement error (this time with 50 and 100 samples)
curr_B <- 1
nsampcurr <- 50

curr_H <- 0
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B)
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT), fill = TRUE)

curr_H <- 0.01
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT), fill = TRUE)

curr_H <- 0.50
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT), fill = TRUE)

curr_H <- 1
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT), fill = TRUE)

curr_H <- 2
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT), fill = TRUE)

curr_H <- 4
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT), fill = TRUE)
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT), fill = TRUE)


# curr_H <- 10
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT), fill = TRUE)



# curr_H <- 20
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT), fill = TRUE)





curr_B <- 1
nsampcurr <- 100

curr_H <- 0
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT), fill = TRUE)

curr_H <- 0.01
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT), fill = TRUE)

curr_H <- 0.50
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT), fill = TRUE)

curr_H <- 1
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT), fill = TRUE)

curr_H <- 2
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT), fill = TRUE)

curr_H <- 4
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
#RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT), fill = TRUE)


# curr_H <- 10
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
#RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT), fill = TRUE)



# curr_H <- 20
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B))
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
# RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B))
# #RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT <- rbind(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT, InsertBlankLine(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT), fill = TRUE)



RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT2 <- data.frame(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT)
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT2$H <- fr(as.numeric(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT2$H), 2, InsertLessThanFormatting = F)
#RSimsPowerTableSuppMIT$H <- NULL
for(j in 5:ncol(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT2)){
  RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT2[,j] <- fr(as.numeric(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT2[,j]),3)
}
RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSamp <- data.table(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT2)
colnames(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSamp) <- colnames(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSampT)
save(RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSamp, file = paste0(save_dir, "/", "RSimsPowerTableUpdatedOrderingSuppIncMeasErrorAddNSamp.RData"))
rm(curr_H)
rm(curr_B)
rm(nsampcurr)


curr_B <- 1
nsampcurr <- 26
curr_H <- 1
#Now, tables for increased between-subject variance
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE)
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, InsertBlankLine(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, IncBetVarTable = TRUE), fill = TRUE)

curr_B <- 1.50
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, InsertBlankLine(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, IncBetVarTable = TRUE), fill = TRUE)

curr_B <- 2.00
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, InsertBlankLine(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, IncBetVarTable = TRUE), fill = TRUE)

curr_B <- 4.00
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, InsertBlankLine(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, IncBetVarTable = TRUE), fill = TRUE)
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, InsertBlankLine(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, IncBetVarTable = TRUE), fill = TRUE)


# curr_B <- 10.00
# RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
# RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
# RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
# RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, InsertBlankLine(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, IncBetVarTable = TRUE), fill = TRUE)

# curr_B <- 20.00
# RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
# RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
# RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
# RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, InsertBlankLine(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, IncBetVarTable = TRUE), fill = TRUE)


curr_B <- 1
nsampcurr <- 100
#Now, tables for increased between-subject variance
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, InsertBlankLine(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, IncBetVarTable = TRUE), fill = TRUE)

curr_B <- 1.50
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, InsertBlankLine(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, IncBetVarTable = TRUE), fill = TRUE)

curr_B <- 2.00
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, InsertBlankLine(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, IncBetVarTable = TRUE), fill = TRUE)

curr_B <- 4.00
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
#RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, InsertBlankLine(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, IncBetVarTable = TRUE), fill = TRUE)

# curr_B <- 10.00
# RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
# RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
# RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
#RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, InsertBlankLine(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, IncBetVarTable = TRUE), fill = TRUE)

# curr_B <- 20.00
# RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = NA, curr_methd = "CompDTU", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
# RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 50, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
# RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, PowerTableFn(nsampcurr = nsampcurr, nreps_curr = 100, curr_methd = "CompDTUme", curr_H = curr_H, curr_B = curr_B, IncBetVarTable = TRUE))
# #RSimsPowerTableUpdatedOrderingIncBetSubjVarT <- rbind(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, InsertBlankLine(RSimsPowerTableUpdatedOrderingIncBetSubjVarT, IncBetVarTable = TRUE), fill = TRUE)

RSimsPowerTableUpdatedOrderingIncBetSubjVarT2 <- data.frame(RSimsPowerTableUpdatedOrderingIncBetSubjVarT)
RSimsPowerTableUpdatedOrderingIncBetSubjVarT2$B <- fr(as.numeric(RSimsPowerTableUpdatedOrderingIncBetSubjVarT2$B), 2, InsertLessThanFormatting = F)
#RSimsPowerTableSuppMIT$H <- NULL
for(j in 5:ncol(RSimsPowerTableUpdatedOrderingIncBetSubjVarT2)){
  RSimsPowerTableUpdatedOrderingIncBetSubjVarT2[,j] <- fr(as.numeric(RSimsPowerTableUpdatedOrderingIncBetSubjVarT2[,j]),3)
}
RSimsPowerTableSuppIncBetSubjCovDiagonalsUpdatedOrdering <- data.table(RSimsPowerTableUpdatedOrderingIncBetSubjVarT2)
colnames(RSimsPowerTableSuppIncBetSubjCovDiagonalsUpdatedOrdering) <- colnames(RSimsPowerTableUpdatedOrderingIncBetSubjVarT)
save(RSimsPowerTableSuppIncBetSubjCovDiagonalsUpdatedOrdering, file = paste0(save_dir, "/", "RSimsPowerTableSuppIncBetSubjCovDiagonalsUpdatedOrdering.RData"))
rm(curr_H)
rm(curr_B)
rm(nsampcurr)
