#q is df of condition variable
calcPillaiPval <- function(SigmaTildeNull, SigmaTildeAlt, res1 = NULL, q, nsamp, GibbsCovs = NULL, df_residual = NA, fullReturn = FALSE, returnTestStatOnly = FALSE){
  
  
  if(is.null(res1) & is.na(df_residual)){
    stop("df. residual must be specified to calcPillaiPval or else an lm object must specified in res1 to extract df.residual from")
  }
  
  Etilde <- nsamp * SigmaTildeAlt
  Htilde <- nsamp * (SigmaTildeNull - SigmaTildeAlt)
  
  vals <- tryCatch(diag((Htilde %*%solve(Etilde + Htilde))), error = function(x){})
  if(is.null(vals)){
    return(NA)
  }
  pill_stat <- sum(vals)
  
  
  #See the Multivariate ANOVA Testing pdf document (from the SAS help file) for the necessary formulas
  
  #v <- nsamp - ncol(model.matrix(res1))
  
  if(!is.na(df_residual)){
    v <- df_residual
  }else{
    v <- res1$df.residual
  }
  #v is the error/residual df- also extract from the r anova fit
  
  
  #p is the number of eigenvales (ie rank of Etilde + Htilde)
  p <- length(vals)
  
  s <- min(p,q)
  
  m <- 0.5 * (abs(p - q) - 1)
  n <- 0.5 * (v - p - 1)
  
  #Formulas come from the SAS help file Multivariate ANOVA Testing
  piece1 <- 2*n + s + 1
  piece2 <- 2*m + s + 1
  fstat_pillai <- (piece1/piece2) * (pill_stat/(s - pill_stat))
  
  if(fstat_pillai < 0){
    return(NA)
  }
  
  
  if(returnTestStatOnly==TRUE){
    return(fstat_pillai)
  }
  numdf_pillai <- s * piece2
  denomdf_pillai <- s * piece1
  pval_pillai <- 1-pf(fstat_pillai, df1 = numdf_pillai, df2 = denomdf_pillai)
  
  if(fullReturn==TRUE){
    return(list(pval_pillai = pval_pillai, pill_stat = pill_stat, df_residual = v, Etilde = Etilde, Htilde = Htilde, numdf_pillai = numdf_pillai, denomdf_pillai = denomdf_pillai))
  }else{
    return(pval_pillai)
  }
  
}

#CorrectLowExpression is not needed for these analyses because that is done on the proportion scale before conversion to ilr scale 
  #and the data here is generated on the ilr scale directly
RSimulationCode <- function(x, fil, nsamp, nreps_to_use_list, WithinSubjCovSampType = "Permute",
                            Wishdf = NULL, seedstouse, ncondlevels = NULL, condEffSizeStr = NA,
                            MeasErrorMultFactor = NA, BetweenCovMultFactor = NA, GenWithMeasError = NA,
                            useRobustCovEst = FALSE, GibbsCovsToUse, res1, res1null, SigmaTildeAlt,
                            SigmaTildeNull, Y, cond, gene_id){
  #load(fil)

  #Set seed to current simulation
  curr_seed <- seedstouse[x]
  
  #1 needs to be in the nreps_to_use list because 1 Rep corresponds to CompDTU and CompDTUme needs these results to run
    #In particular 1 should be the first nreps element specified
  if(!(1 %in% nreps_to_use_list)){
    nreps_to_use_list <- c(1, nreps_to_use_list)
  }
  
  if(nreps_to_use_list[1]!=1){
    nreps_to_use_list <- c(1, nreps_to_use_list)
  }
  
  library(matrixcalc)
  if(x %% 10 == 0){
    print(paste0("Current simuation number is ", x))
  }
  
  #Load initial values based on the data for the observed gene that is the basis of rht simulation values
  XNull <- model.matrix(res1null)
  X <- model.matrix(res1)
  
  nsamporig <- nrow(X)
  
  betNull <- coef(res1null)
  MeanNull <- XNull%*%betNull
  MeanVNull <- as.vector(t(MeanNull))
  
  
  bet <- coef(res1)
  MeanAlt <- X%*%bet
  MeanVAlt <- as.vector(t(MeanAlt))
  
  k <- ncol(Y)
  q <- ncondlevels - 1
  
  
  #Set the number of conditions the simulated data are to be generated from if it is not specified to 2
  if(is.na(ncondlevels)==TRUE){
    ncondlevels <- nlevels(cond)
  }
  
  if(ncondlevels > nlevels(cond)){
    stop(paste0("You have specified more unique conditions than this data can support.  The maximum number of conditions can be ", nlevels(cond)))
  }
  
  #Set effect size vectors for the conditions
  condEffSize <- eval(parse(text = condEffSizeStr))
  
  if(length(condEffSize)!=ncondlevels){
    stop(paste0("The effectsizes much be specified for each level of the condition, in this case there need to be ", ncondlevels, " levels"))
  }
  
  #Generate a dataset with the maximum desired number of replicates that can be subsetted to obtain fewer replicates
  nreps_sim <- max(nreps_to_use_list)
  SimData <- GenerateRSimDataset(curr_seed = curr_seed, fil = fil, nsamp = nsamp, nreps_sim = nreps_sim, 
                                  GenWithMeasError = GenWithMeasError, condEffSize = condEffSize,
                                  condEffSizeStr = condEffSizeStr, ncondlevels = 2, WithinSubjCovSampType = "Permute",
                                  MeasErrorMultFactor = MeasErrorMultFactor, BetweenCovMultFactor = BetweenCovMultFactor,
                                  Wishdf = Wishdf)
  
  #Full simulated dataset
  #browser()
  datSim <- SimData$datSim
  
  #key with one row per sample (with Sample identifier and condition assignment)
  key_sim <- SimData$key_sim
  
  #key with one row per replicate (with Sample identifier and condition assignment)
  key_simrep <- SimData$key_simrep
  
  #Number of samples to set for each condition
  nsampspercond <- ceiling(nsamp/ncondlevels)
  
  #Obtain key for only those samples that are used in this analysis
    #This allows the code to accommodate a case that only a subset of samples were used, 
    #such as would occur is nsamp is not evenly divisible by ncondlevels
  samps_to_use <- list()
  condsim <- key_sim$condsim
  for(j in 1:length(levels(condsim))){
    curr_cond <- levels(condsim)[j]
    sub_key_sim <- subset(key_sim, condsim==curr_cond)
    
    #Now, if the number of samples to use per condition (say C) is less than the total number, take the first C from the total
    samps_to_use_curr_cond <- sub_key_sim[1:nsampspercond,]
    samps_to_use[[j]] <- samps_to_use_curr_cond
  }
  
  #key corresponding to the samples being used for the current analysis
    #These do not need to be ordered, as the key will be compared to the already ordered key_simrep at the 4th line of the loop below
  key_samps_to_use <- data.frame(data.table::rbindlist(samps_to_use))
  
  #Repeat analysis for each nrep of interest
  for(j in 1:length(nreps_to_use_list)){
    curr_nrep <- nreps_to_use_list[j]
    #print(paste0("curr_nrep is ", curr_nrep))
    #Subset the full key to only correspond to Samples and replicates in use by the current analysis
    assign(paste0("key_simrep_sub", curr_nrep, "Rep"), subset(key_simrep, (key_simrep$repNum %in% 1:curr_nrep & key_simrep$Sample %in% key_samps_to_use$Sample)))
    curr_key_simrep_sub <- get(paste0("key_simrep_sub", curr_nrep, "Rep"))
    
    #Current condition factor for the current analysis
    assign(paste0("condsimrep", curr_nrep, "Rep"), curr_key_simrep_sub$condsim)
    curr_condsimrep <- get(paste0("condsimrep", curr_nrep, "Rep"))
    
    #Subset of simulated data datSim corresponding to the currently used data
    assign(paste0("datSimToUse", curr_nrep, "Rep"), datSim[(key_simrep$Sample %in% key_samps_to_use$Sample & key_simrep$repNum %in% 1:curr_nrep),])
    curr_dat <- get(paste0("datSimToUse", curr_nrep, "Rep"))
    #browser()
    #Covariance of simulated data- not currently used
    #assign(paste0("SimDataObsCov", curr_nrep, "Rep"), cov(curr_dat))
    
    #Fit CompDTU model
    assign(paste0("model1null", curr_nrep, "Rep"), lm(curr_dat ~ 1))
    assign(paste0("model1alt", curr_nrep, "Rep"), lm(curr_dat ~ curr_condsimrep))
    curr_nullmodel <- get(paste0("model1null", curr_nrep, "Rep"))
    curr_altmodel <- get(paste0("model1alt", curr_nrep, "Rep"))
    
    assign(paste0("nobs_null", curr_nrep, "Rep"), nobs(curr_nullmodel))
    assign(paste0("nobs_alt", curr_nrep, "Rep"), nobs(curr_altmodel))
    curr_nobs_null <- get(paste0("nobs_null", curr_nrep, "Rep"))
    curr_nobs_alt <- get(paste0("nobs_alt", curr_nrep, "Rep"))
    
    #MANOVA cov matrix - see slide 78 of Helwig's MANOVA notes from Minnesota
    #This is the exact MLE from MANOVA in case of no subject specific covariance
    #Can calculate both of these quantities whether data is generated under null or alternative,
    #though which one you want to compare everything to will depend on how data is generated
    #Denom is nobs(model1null) or nobs(model1alt) not nsamp because if nreps > 1 this needs to be the total
    #of unique observations input to the lm fit
    assign(paste0("simMatMANOVANull", curr_nrep, "Rep"), cov(as.matrix(curr_nullmodel$residuals, ncol = k)) * ((curr_nobs_null-1)/curr_nobs_null))
    assign(paste0("simMatMANOVAAlt", curr_nrep, "Rep"), cov(as.matrix(curr_altmodel$residuals, ncol = k)) * ((curr_nobs_alt-1)/curr_nobs_alt))
    curr_simMatMANOVANull <- get(paste0("simMatMANOVANull", curr_nrep, "Rep"))
    curr_simMatMANOVAAlt <- get(paste0("simMatMANOVAAlt", curr_nrep, "Rep"))
    #print("Step 1")
    if(curr_nrep==1){
      pval_CompDTU <- calcPillaiPval(SigmaTildeNull = simMatMANOVANull1Rep, SigmaTildeAlt = simMatMANOVAAlt1Rep, 
                                     res1 = model1alt1Rep, q = q, nsamp = nobs(model1alt1Rep))
    }else{
      #Generate within sample covariance of simulated data for each sample
      WithinSubjCovofSimData <- vector(mode = "list", length = nsamp)
      MeanSampledDataT <- vector(mode = "list", length = nsamp)
      MedianSampledDataT <- vector(mode = "list", length = nsamp)
      rnames <- rownames(curr_dat)
      rnames2 <- strsplit(rnames, "Rep")
      rnames3 <- lapply(rnames2, function(x){x[1]})
      rnames4 <- as.character(rnames3)
      #subj_spec_dat <- list()
      #print("Step 2")
      for(i in 1:nsamp){
        curr_samp <- paste0("Sample", i)
        curr_rows <- rnames4==curr_samp
        curr_cov <- NULL
        sub_dat_curr_samp <- curr_dat[curr_rows,]
        #subj_spec_dat[[i]] <- sub_dat_curr_samp
        #browser()
        if(useRobustCovEst==TRUE){
          RobCovEstT <- tryCatch(robust::covRob(matrix(sub_dat_curr_samp, ncol = ncol(datSim))), error = function(x){})
          if(is.null(RobCovEstT)){
            curr_cov <- cov(sub_dat_curr_samp)
          }
          
          if(is.null(curr_cov)){
            RobCovEst <- RobCovEstT$cov
            if(ncol(RobCovEst)!=ncol(datSim)){
              curr_cov <- cov(sub_dat_curr_samp)
            }else{
              curr_cov <- RobCovEst
              #print(paste0("The robust covariance estimator is used for Sample ", i))
            }
            
          }
        }else{
          curr_cov <- cov(sub_dat_curr_samp)
        }
        curr_mean <- colMeans(sub_dat_curr_samp)
        curr_median <- matrixStats::colMedians(sub_dat_curr_samp)
        MeanSampledDataT[[i]] <- data.frame(matrix(curr_mean, nrow = 1, ncol = ncol(sub_dat_curr_samp)))
        MedianSampledDataT[[i]] <- data.frame(matrix(curr_median, nrow = 1, ncol = ncol(sub_dat_curr_samp)))
        WithinSubjCovofSimData[[i]] <- curr_cov
      }
      
      #browser()
      #save(WithinSubjCovofSimData, subj_spec_dat, file = "~/tempCov.RData")
      #stop()
      
      #First, calculate CompDTUme results
      #print("Step 3")
      #mean of within sample covariance
      curr_mean.withinhat <-  Reduce("+", WithinSubjCovofSimData)/nsamp
      
      #Updated covariance terms
      curr_UpdatedCovAlt <- curr_simMatMANOVAAlt - curr_mean.withinhat
      curr_UpdatedCovNull <- curr_simMatMANOVANull - curr_mean.withinhat

      assign(paste0("pval_CompDTUme", curr_nrep, "Rep"), calcPillaiPval(SigmaTildeNull = curr_UpdatedCovNull, SigmaTildeAlt = curr_UpdatedCovAlt,
                                                                        res1 = NULL, q = q, nsamp = nsamp, df_residual = nsamp - ncol(model.matrix(model1alt1Rep))))
      
      
      #print("Step 4")
      #Now, calculate CompDTU results for values that are the mean/median of the sampled replicated
      MeanSampledData <- as.matrix(data.frame(data.table::rbindlist(MeanSampledDataT)))
      MedianSampledData <- as.matrix(data.frame(data.table::rbindlist(MedianSampledDataT)))
      rownames(MeanSampledData) <- paste0("Sample", 1:nsamp)
      rownames(MedianSampledData) <- paste0("Sample", 1:nsamp)
      colnames(MeanSampledData) <- NULL
      colnames(MedianSampledData) <- NULL
      
      #Previous results calculated the two covariance matricies and the Pillai statistic by hand, but
        #can just use anova function now since these matricies will not be needed by future calculations
      anova_MeanSampled <- tryCatch(anova(lm(MeanSampledData ~ condsim), test = "Pillai"))
      
      if(is.null(anova_MeanSampled)){
        assign(paste0("pval_CompAbFromMeanInfRep", curr_nrep, "Rep"), NA)
      }else{
        assign(paste0("pval_CompAbFromMeanInfRep", curr_nrep, "Rep"), anova_MeanSampled["condsim", "Pr(>F)"])
      }
      
      #print("Step 6")
      rm(MeanSampledData, MedianSampledData, anova_MeanSampled, anova_MedianSampled)
      
    }#end overall else loop that is run if curr_nrep is not equal to 1
  }
  

      objtoretT <-c(ls(pattern = "pval_CompDTU"), ls(pattern = "pval_CompAbFromMeanInfRep"))
    
    objtoret <- data.frame(mget(objtoretT))
    return(objtoret)
    
}



#CorrectLowExpression is not needed for these analyses because those corrections are done on the proportion scale before conversion to ilr scale 
  #and the data here is generated on the ilr scale directly
GenerateRSimDataset <- function(curr_seed, fil, nsamp, nreps_sim, GenWithMeasError, condEffSize, condEffSizeStr, ncondlevels = 2, 
                                 WithinSubjCovSampType = "Permute", MeasErrorMultFactor = 1, BetweenCovMultFactor = 1,
                                 Wishdf = NULL, save_dir){
  set.seed(curr_seed)
  load(fil)
  
  #Load initial values based on the data for the observed gene
  XNull <- model.matrix(res1null)
  X <- model.matrix(res1)
  
  nsamporig <- nrow(X)
  
  betNull <- coef(res1null)
  MeanNull <- XNull%*%betNull
  MeanVNull <- as.vector(t(MeanNull))
  
  
  bet <- coef(res1)
  MeanAlt <- X%*%bet
  MeanVAlt <- as.vector(t(MeanAlt))
  
  k <- ncol(Y)
  q <- ncondlevels - 1
  
  #Assign the level to be the first 2 levels from the E-GEUV-1 data condition
    #These are just labels, could also be "Cond1", "Cond2", etc
  levelstouse <- levels(cond)[1:ncondlevels]
  ###########################
  #Create condition variable
  #Condition needs to be created in this way to ensure an even number of samples in each condition with more than 2 condition levels
  #I tried simpler approaches, but they didn't always guarantee an even number of samples per condition, especially when
  #using a small number of samples
  
  
  #Number of samples per condition to generate
  nsampspercond <- ceiling(nsamp/ncondlevels)
  samps_vec <- 1:nsamp
  samps_curr_cond <- vector(mode = "list", length = ncondlevels)
  samps_already_assigned <- vector(mode = "list", length = ncondlevels)

  for(j in 1:ncondlevels){
    #For the last group, assign the remaining samples that have not yet been assigned
    if(j==ncondlevels){
      samps_to_assign <- samps_vec[!(samps_vec) %in% samps_already_assigned]
      if(length(samps_to_assign)!=nsampspercond){
        stop("Something is wrong with the assigned samples")
      }
      samps_curr_cond[[j]] <- samps_to_assign
      rm(samps_to_assign)
      next
    }

    #Sample numbers to assign to the current condition
    samps_curr_cond[[j]] <- sample(samps_vec[!(samps_vec %in% samps_already_assigned)], nsampspercond, replace = FALSE)

    #Remove samples already assigned above to make sure they do not get assigned to multiple groups
    samps_already_assigned <- c(samps_already_assigned, samps_curr_cond[[j]])
  }

  #Create condition vector, with level values that are assigned to be the same as from the input cond
    #ie CEU, FIN, etc from the E-GEUV-1 data
  condsim <- factor(character(length = nsamp), levels = levelstouse)
  for(j in 1:ncondlevels){
    condsim[samps_curr_cond[[j]]] <- levelstouse[j]
  }

  SampNames <- paste0("Sample", 1:nsamp)
  key_sim <- data.frame(condsim, SampNames, stringsAsFactors = F)
  colnames(key_sim) <- c("condsim", "Sample")
  
  # End create condition variable
  ###########################
  
  MeanAlt2 <- unique(MeanAlt)
  MeanAlt3 <- MeanAlt2
  rownames(MeanAlt3) <- as.character(strsplit(rownames(MeanAlt2), split = "TPM"))
  MeanAlt4 <- MeanAlt3
  rownames(MeanAlt4) <- cond[names(cond) %in% rownames(MeanAlt3)]
  
  #Specify Mean Values to use in data simulation
  #Under the alternative these correspond to the specific condition, under the null they are all the same
  s <- length(levelstouse)
  
  #Sample within-sample covariance matrices
  if(WithinSubjCovSampType=="Permute"){
    #Sample covariance matrices with replacement from ones that were in observed data
    PermCov <- vector(mode = "list", length = nsamp)
    m <- length(GibbsCovsToUse)
    PermMatIndicies <- sample(1:m, nsamp, replace = TRUE, prob = rep(1/m, m))
    for(i in 1:nsamp){
      PermCov[[i]] <- GibbsCovsToUse[[PermMatIndicies[i]]]
    }
    SubjSpecCov <- PermCov
    
    #Increase off diagonals to result in increased measurement error assuming MeasErrorMultFactor is specified and not = to 1
    if(!is.na(MeasErrorMultFactor)){
      for(i in 1:nsamp){
        curr_mat <- SubjSpecCov[[i]] * MeasErrorMultFactor
        SubjSpecCov[[i]] <- curr_mat
        
        #Older approach that only multiplies diagonal by the MeasErrorMultFactor, not the whole matrix
        # for(j in 1:k){
        #   curr_mat[j,j] <- curr_mat[j,j] * MeasErrorMultFactor
        # }
      }
    }
    
    
    rm(m)
    
    #Possible method to sampele covariance matrices from the Wishart distribution
  }else if(WithinSubjCovSampType=="Sample"){
    #Use k as the df parameter to let the sampled matricies vary as much as possible
    #or, you can specify another value as one of the input arguments
    #Regardless, use param values of Wishdf, GibbsCovsMean/Wishdf such that the mean of the sampled matricies is GibbsCovsMean
    
    #Calculate mean and sd of covariance matrix of inferential replicates across all samples
    if(!is.null(GibbsCovsToUse)){
      if(dim(GibbsCovsToUse[[1]])[1]==1){
        tempdat <- data.frame()
        for(i in 1:length(GibbsCovsToUse)){
          tempdat[i,1] <- GibbsCovsToUse[[i]]
        }
        GibbsCovsMean <- mean(tempdat[,1])
        GibbsCovsSD <- sd(tempdat[,1])
      }else{
        GibbsCovsMean <- apply(simplify2array(GibbsCovsToUse), 1:2, mean)
        GibbsCovsSD <- apply(simplify2array(GibbsCovsToUse), 1:2, sd)
      }
    }
    
    if(is.na(Wishdf)==TRUE){
      Wishdf <- k
      #Wishdf <- nsamp/2
    }
    
    #Increase off diagonals to result in increased measurement error assuming MeasErrorMultFactor is specified and not 1
    #For this wishart simuation try increasing the diagonals of the matrix the cov matricies are generated from, not the diagonals
    #of the generated matricies themselves because this should help you have more control over the amount of measurement error
    if(!is.na(MeasErrorMultFactor)){
      for(j in 1:k){
        GibbsCovsMean[j,j] <- GibbsCovsMean[j,j] * MeasErrorMultFactor
      }
    }
    
    rCovs <- rWishart(nsamp, Wishdf, GibbsCovsMean/Wishdf)
    rCovsList <- lapply(seq(dim(rCovs)[3]), function(x) rCovs[ , , x])
    SubjSpecCov <-rCovsList
    
  }else{
    stop("WithinSubjCovSampType specification is invalid")
  }
  
  
  
  #Generate matrix of mean values to use in the simulation
  MVals <- matrix(0, nrow = nsamp, ncol = k)
  
  #Make sure at least one of the elements in condEffSize is 1 to have a null group that has is unchanged from the values
    #of the observed data
    #Usually this will be the first element, but it wouldn't have to be
    #If condEffSize is 1 for all j, all mean values (MVals) are the same, corresponding to the null hypothesis
  for (j in 1:s){
    currv <- matrix(MeanNull[1,] * condEffSize[j], nrow = 1, ncol = length(MeanNull[1,]))
    ind <- which(condsim==levelstouse[j])
    for(z in 1:length(ind)){
      MVals[ind[z],] <- currv
    }
  }
  
  #Simulate all runs under the same between-sample covariance, using the values from CompDTU under the null for the observed data
  SigmaTildeToUse <- SigmaTildeNull
  
  #Increase diagonals of between subject covariance matricies to inc variance of between-subject cov matricies
  if(!is.na(BetweenCovMultFactor)){
    for(j in 1:k){
      SigmaTildeToUse[j,j] <- SigmaTildeToUse[j,j] * BetweenCovMultFactor
    }
  }
  datSim <- matrix(NA, nrow = nsamp*nreps_sim, ncol = k)
  if(GenWithMeasError==TRUE){
    for(i in 1:nsamp){
      #Sample between-subject residual
      between.sub.e = mvtnorm::rmvnorm(n = 1, sigma = SigmaTildeToUse) # one btw subject residual per subject
      between.sub.e.to.use <- matrix(NA, nrow = nreps_sim, ncol = ncol(between.sub.e))
      
      #Repeat between-subject residual for each replicate for that individual
      for(j in 1:nreps_sim){
        between.sub.e.to.use[j,] <- between.sub.e
      }
      start <- (nreps_sim*(i-1) + 1)
      end <- nreps_sim*i
      datSim[start:end,] <- mvtnorm::rmvnorm(nreps_sim, mean = MVals[i,], sigma = SubjSpecCov[[i]]) + between.sub.e.to.use
    }
  }else if(GenWithMeasError==FALSE){
    #Sample data without any within-subject error
    for(i in 1:nsamp){
        start <- (nreps_sim*(i-1) + 1)
        end <- nreps_sim*i
        datSim[start:end,] <- mvtnorm::rmvnorm(nreps_sim, mean = MVals[i,], sigma = SigmaTildeToUse)
      
    }
  }
  
  #Create key corresponding to each row of the simulated data (datSim) that has each element of key_sim repeated nreps times at its current position
    #This is done by repeating everything nreps times and importantly sorting by Sample
    #Use mixedsort in gtools such that Sample100 is after Sample99, etc
  key_simrep_temp <- data.frame(rep(key_sim$condsim, nreps_sim), rep(key_sim$Sample, nreps_sim), stringsAsFactors = F)
  colnames(key_simrep_temp) <- c("condsim", "Sample")
  key_simrep <- key_simrep_temp[gtools::mixedorder(key_simrep_temp$Sample),]
  #A way to do the above without gtools, though this is much slower than the gtools version
  #key_simrep2 <- key_simrep[rep(seq_len(nrow(key_simrep)), each = nreps_sim), ]
  rownames(key_simrep) <- 1:nrow(key_simrep)
  
  key_simrep$repNum <- rep(1:nreps_sim, nsamp)
  
  rownames(datSim) <- paste0(key_simrep$Sample, "RepNum", key_simrep$repNum)
  
  return(list(datSim = datSim, key_sim = key_sim, key_simrep = key_simrep))
}