CompDTUMethodsPowerAnalysis <- function(genes, curr_cond, onCluster,
                                      AllGroupCombinations, curr_change, CalculateCompDTUmeRes, nsamp, nboot,
                                      ngrpcombos = 252, BootSamps = TRUE, samps = NULL, GEUV1Data = FALSE,
                                      txsfiltered = NULL, 
                                      CalcCompDTUResForAbFromInfReps = FALSE,
                                      UpdatedAbDatasetsAbFromInfRepsDir = NULL,
                                      DatMeanBoot = NULL){
  
  print(paste0("nsamp is ", nsamp))
  print(paste0("nboot/ngibbs is ", nboot))
  print(paste0("Length of samps is ", length(samps)))
  print(paste0("First 12 samps are printed below"))
  print(head(samps, 12))
  res <- vector(mode = "list", length = length(genes))
  if(is.null(samps)){
    samps <- paste0("Sample", 1:nsamp)
  }
  
  if(length(samps)!=nsamp){
    stop("The length of samps doesn't match nsamp.  Something must be specified incorrectly.")
  }
  for(l in 1:length(genes)){
    print(paste0("Current gene number is ", l, " out of ", length(genes)))
      curr_gene <- genes[l]
      #print(paste0("l is ", l))
      
      if(is.null(txsfiltered)){
        #Non inferential replicate data
        dObs <- unclass(ilr(CorrectLowExpression(NonBootData[[curr_gene]][paste0(samps, "TPM"),])))

        #Bootstrap data replicate data
        
        bootrownames <- vector(mode = "character", length = nsamp*nboot)
        if(BootSamps==TRUE){
          str1 <- "Boot"
        }else{
          str1 <- "Gibbs"
        }
        for(q in 1:nboot){
          start <- 1 + (q-1)*nsamp
          end <- start + (nsamp - 1)
          curr_bo <- q
          bootrownames[start:end] <- paste0(samps, "TPM", str1, curr_bo)
        }
        
        dboot <- unclass(ilr(CorrectLowExpression(BootData[[curr_gene]][bootrownames,])))
        
      }else{
        #Non inferential replicate data
        cu <- NonBootData[[curr_gene]]
        cuu <- cu[,colnames(cu) %in% txsfiltered, drop = FALSE]
        cuuu <- cu[paste0(samps, "TPM"),]
        
        if(ncol(cuuu)==1){
          res[[l]] <- list(pval_CompDTU = NA)
          next
        }    
    
        dObs <- unclass(ilr(CorrectLowExpression(cuuu)))
        
        
        if(l==1){
          print("Head of the non bootstrap dataset is given below.  The data should be ordered Sample4TPM, Sample7TPM, etc")
          print(head(rownames(dObs), 12))
        }
        
        
        #Bootstrap data replicate data
        cu2 <- BootData[[curr_gene]]
        cuu2 <- cu2[,colnames(cu2) %in% txsfiltered, drop = FALSE]
        
        bootrownames <- vector(mode = "character", length = nsamp*nboot)
        if(BootSamps==TRUE){
          str1 <- "Boot"
        }else{
          str1 <- "Gibbs"
        }
        for(q in 1:nboot){
          start <- 1 + (q-1)*nsamp
          end <- start + (nsamp - 1)
          curr_bo <- q
          bootrownames[start:end] <- paste0(samps, "TPM", str1, curr_bo)
        }
        
        cuuu2 <- cuu2[bootrownames,]
        if(ncol(cuuu2)==1){
          res[[i]] <- list(pval_CompDTU = NA)
          next
        }  
        
        dboot <- unclass(ilr(CorrectLowExpression(cuuu2)))
        
      }
      
      

    
    
    if(l==1){
      print("Head of the bootstrap/gibbs sample dataset is given below.  The data should be ordered Sample4TPMBoot (or  Sample4TPMGibbs), Sample7TPMBoot1 (or  Sample7TPMGibbs), etc")
      print(head(rownames(dboot), 12))
    }
    
    if(nrow(dObs)!=nsamp | nrow(dboot)!=(nsamp*nboot)){
      stop("The datasets don't have the right dimensions, something is wrong.")
    }
    
    nc <- ncol(dObs)
    
    CompDTUPowerFn <- function(data, curr_cond, returnTestStat = FALSE){
      res_Obs <- lm(data ~ curr_cond)
      anova_res_Obs <- tryCatch(anova(res_Obs), error = function(x){})
      
      if(is.null(anova_res_Obs)){
        test_stat_CompDTU <- NA
      }else{
        if(ncol(data)==1){
          test_stat_CompDTU <- anova_res_Obs["curr_cond", "F value"]
          df1 <- anova_res_Obs["curr_cond", "Df"]
          df2 <- anova_res_Obs["Residuals", "Df"]
          
          pval_CompDTU <- anova_res_Obs["curr_cond", "Pr(>F)"]
        }else{
          test_stat_CompDTU <- anova_res_Obs["curr_cond", "approx F"]
          df1 <- anova_res_Obs["curr_cond", "num Df"]
          df2 <- anova_res_Obs["curr_cond", "den Df"]
          
          pval_CompDTU <- anova_res_Obs["curr_cond", "Pr(>F)"]
        }
      }
      
      if(returnTestStat==TRUE){
        return(test_stat_CompDTU)
      }else{
        return(pval_CompDTU)
      }
      
    }
    
    pval_CompDTU <- CompDTUPowerFn(dObs, curr_cond)
    
    test_stat_CompDTU <- CompDTUPowerFn(dObs, curr_cond, returnTestStat = TRUE)
    
    if(CalcCompDTUResForAbFromInfReps==TRUE){
      
      RunResForAbFromInfReps <- function(curr_gene, DataToUse, samps, curr_cond){
        if(is.null(txsfiltered)){
          dObs <- unclass(ilr(CorrectLowExpression(DataToUse[[curr_gene]][paste0(samps, "TPM"),])))
        }else{
          cu <- DataToUse[[curr_gene]]
          cuu <- cu[,colnames(cu) %in% txsfiltered, drop = FALSE]
          cuuu <- cu[paste0(samps, "TPM"),]
          
          if(ncol(cuuu)==1){
            res[[l]] <- list(pval_CompDTU = NA)
            next
          }    
          
          dObs <- unclass(ilr(CorrectLowExpression(cuuu)))
        }
        
        res_Obs <- lm(dObs ~ curr_cond)
        anova_res_Obs <- tryCatch(anova(res_Obs), error = function(x){})
        if(is.null(anova_res_Obs)){
          pval_CompDTU <- NA
        }else{
          pval_CompDTU <- anova_res_Obs["curr_cond", "Pr(>F)"]
        }
        
        return(pval_CompDTU)
      }
      
      pval_CompDTUFromMeanBoot <- RunResForAbFromInfReps(curr_gene = curr_gene, DataToUse = DatMeanBoot, samps = samps, curr_cond = curr_cond)
      
      
      res[[l]] <- data.frame(gene_id = curr_gene, 
                             pval_CompDTU = pval_CompDTU, 
                             pval_CompDTUFromMeanBoot = pval_CompDTUFromMeanBoot)
      next
      
    }
    
    
    
    if(CalculateCompDTUmeRes==TRUE){
      #curr_WithinCov <- list()
      #To get the ordered row names
      rows_mixedsort <- mixedsort(rownames(dboot))
      CompDTUmePowerFn <- function(data, curr_cond, rows_mixedsort, nsamp, nboot, NonBootData, returnTestStatistics = FALSE){
        #First, fit nonbootstrap model to the non bootstrap data just to extract the necessary df value for the pvalue calculation
        res_Obs <- lm(NonBootData ~ curr_cond)
        
        curr_WithinCov <- vector(mode = "list", length = nsamp)
        #ct <- proc.time()
        
        curr_e4 <- data[rows_mixedsort, , drop = FALSE]
        for(j in 1:nsamp){
          start <- nboot*(j-1) +1
          end <- nboot*j
          #print(paste0("start is ", start, " and end is ", end))
          curr_e5 <- curr_e4[start:end, , drop = FALSE]
          #print(rownames(curr_e5)[1:5])
          curr_WithinCov[[j]] <- cov(curr_e5)
        }
        
        #proc.time() - ct
        
        MeanWithinCov <- Reduce("+", curr_WithinCov)/nsamp
        
        boot_alt_res <- lm(data ~ rep(curr_cond, nboot))
        boot_null_res <- lm(data ~ 1)
        
        
        
        TotalCovAlt <- crossprod(data - boot_alt_res$fitted.values)/nrow(data)
        TotalCovNull <- crossprod(data - boot_null_res$fitted.values)/nrow(data)
        
        
        UpdatedCovAlt <- TotalCovAlt - MeanWithinCov
        UpdatedCovNull <- TotalCovNull - MeanWithinCov
        
        #qq is the df of the condition variable, which is 1 here since there are 2 levels
        qq <- 1
        
        if(returnTestStatistics==TRUE){
          Obs_TestStatCompDTUme <- calcPillaiPval(SigmaTildeNull = UpdatedCovNull, SigmaTildeAlt = UpdatedCovAlt,
                                                  res1 = res_Obs, q = qq, nsamp = nsamp, fullReturn = FALSE, returnTestStatOnly = TRUE)
          
          Obs_TestStatCompDTUmeNome <- calcPillaiPval(SigmaTildeNull = TotalCovNull, SigmaTildeAlt = TotalCovAlt,
                                                      res1 = res_Obs, q = qq, nsamp = nsamp, fullReturn = FALSE, returnTestStatOnly = TRUE)
          
          return(list(Obs_TestStatCompDTUme = Obs_TestStatCompDTUme, Obs_TestStatCompDTUmeNome = Obs_TestStatCompDTUmeNome))
        }else{
          pval_CompDTUme <- calcPillaiPval(SigmaTildeNull = UpdatedCovNull, SigmaTildeAlt = UpdatedCovAlt,
                                           res1 = res_Obs, q = qq, nsamp = nsamp, fullReturn = FALSE)
          
          pval_CompDTUmeNome <- calcPillaiPval(SigmaTildeNull = TotalCovNull, SigmaTildeAlt = TotalCovAlt,
                                               res1 = res_Obs, q = qq, nsamp = nsamp, fullReturn = FALSE)
          
          return(list(pval_CompDTUme = pval_CompDTUme, pval_CompDTUmeNome = pval_CompDTUmeNome))
        }

      }
      
      CompDTUmeRes <- CompDTUmePowerFn(dboot, curr_cond, rows_mixedsort, nsamp, nboot, dObs)
      
      
      pval_CompDTUme <- CompDTUmeRes$pval_CompDTUme
      pval_CompDTUmeNome <- CompDTUmeRes$pval_CompDTUmeNome
    
      
      res[[l]] <- list(pval_CompDTU = pval_CompDTU, 
                         pval_CompDTUme = pval_CompDTUme, 
                         pval_CompDTUmeNome = pval_CompDTUmeNome)
      
    }#end loop for CalculateCompDTUmeRes
    

  }# end l loop indexing over genes
  
  res2 <- data.table::rbindlist(res, fill = TRUE)
  res3 <- data.frame(res2)
  rownames(res3) <- genes
  return(res3)
}
