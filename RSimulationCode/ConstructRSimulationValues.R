#Construct data.frame of simulation values to use for R simulations

library(data.table)

def_wd <- "/Users/Scott/Documents/Dissertation/Paper1/SimulationCode/RSimulationCode"
setwd(def_wd)

RSimVals <-                 data.table(nsamp = 10, ncondlevels = 2, condEffSizeStr = "c(1,1)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1)
RSimVals <- rbind(RSimVals, data.table(nsamp = 100, ncondlevels = 2, condEffSizeStr = "c(1,1)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1))
RSimVals <- rbind(RSimVals, data.table(nsamp = 250, ncondlevels = 2, condEffSizeStr = "c(1,1)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1))

RSimVals <-  rbind(RSimVals, data.table(nsamp = 10, ncondlevels = 2, condEffSizeStr = "c(1,1.025)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1))
RSimVals <- rbind(RSimVals, data.table(nsamp = 100, ncondlevels = 2, condEffSizeStr = "c(1,1.025)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1))
RSimVals <- rbind(RSimVals, data.table(nsamp = 250, ncondlevels = 2, condEffSizeStr = "c(1,1.025)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1))

RSimVals <-  rbind(RSimVals, data.table(nsamp = 10, ncondlevels = 2, condEffSizeStr = "c(1,1.05)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1))
RSimVals <- rbind(RSimVals, data.table(nsamp = 100, ncondlevels = 2, condEffSizeStr = "c(1,1.05)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1))
RSimVals <- rbind(RSimVals, data.table(nsamp = 250, ncondlevels = 2, condEffSizeStr = "c(1,1.05)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1))

RSimVals <-  rbind(RSimVals, data.table(nsamp = 10, ncondlevels = 2, condEffSizeStr = "c(1,1.10)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1))
RSimVals <- rbind(RSimVals, data.table(nsamp = 100, ncondlevels = 2, condEffSizeStr = "c(1,1.10)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1))
RSimVals <- rbind(RSimVals, data.table(nsamp = 250, ncondlevels = 2, condEffSizeStr = "c(1,1.10)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1))

RSimVals <-  rbind(RSimVals, data.table(nsamp = 10, ncondlevels = 2, condEffSizeStr = "c(1,1.25)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1))
RSimVals <- rbind(RSimVals, data.table(nsamp = 100, ncondlevels = 2, condEffSizeStr = "c(1,1.25)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1))
RSimVals <- rbind(RSimVals, data.table(nsamp = 250, ncondlevels = 2, condEffSizeStr = "c(1,1.25)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1))

RSimVals <-  rbind(RSimVals, data.table(nsamp = 10, ncondlevels = 2, condEffSizeStr = "c(1,1.50)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1))
RSimVals <- rbind(RSimVals, data.table(nsamp = 100, ncondlevels = 2, condEffSizeStr = "c(1,1.50)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1))
RSimVals <- rbind(RSimVals, data.table(nsamp = 250, ncondlevels = 2, condEffSizeStr = "c(1,1.50)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1))

RSimVals <-  rbind(RSimVals, data.table(nsamp = 10, ncondlevels = 2, condEffSizeStr = "c(1,2.00)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1))
RSimVals <- rbind(RSimVals, data.table(nsamp = 100, ncondlevels = 2, condEffSizeStr = "c(1,2.00)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1))
RSimVals <- rbind(RSimVals, data.table(nsamp = 250, ncondlevels = 2, condEffSizeStr = "c(1,2.00)", GenWithMeasError = TRUE, MeasErrorMultFactor = 1, BetweenCovMultFactor = 1))

curr_nrow <- nrow(RSimVals)
nsimblocks <- 1

RSimVals <- rbind(RSimVals, RSimVals[((1):(curr_nrow)),])
RSimVals[((nsimblocks*curr_nrow +1):((nsimblocks + 1)*curr_nrow)), "MeasErrorMultFactor"] <- 1.50
nsimblocks <- nsimblocks + 1

RSimVals <- rbind(RSimVals, RSimVals[((1):(curr_nrow)),])
RSimVals[((nsimblocks*curr_nrow +1):((nsimblocks + 1)*curr_nrow)), "MeasErrorMultFactor"] <- 2
nsimblocks <- nsimblocks + 1

RSimVals <- rbind(RSimVals, RSimVals[((1):(curr_nrow)),])
RSimVals[((nsimblocks*curr_nrow +1):((nsimblocks + 1)*curr_nrow)), "MeasErrorMultFactor"] <- 4
nsimblocks <- nsimblocks + 1

RSimVals <- rbind(RSimVals, RSimVals[((1):(curr_nrow)),])
RSimVals[((nsimblocks*curr_nrow +1):((nsimblocks + 1)*curr_nrow)), "MeasErrorMultFactor"] <- 10
nsimblocks <- nsimblocks + 1

RSimVals <- rbind(RSimVals, RSimVals[((1):(curr_nrow)),])
RSimVals[((nsimblocks*curr_nrow +1):((nsimblocks + 1)*curr_nrow)), "MeasErrorMultFactor"] <- 20
nsimblocks <- nsimblocks + 1


RSimVals <- rbind(RSimVals, RSimVals[((1):(curr_nrow)),])
RSimVals[((nsimblocks*curr_nrow +1):((nsimblocks + 1)*curr_nrow)), "BetweenCovMultFactor"] <- 1.50
nsimblocks <- nsimblocks + 1


RSimVals <- rbind(RSimVals, RSimVals[((1):(curr_nrow)),])
RSimVals[((nsimblocks*curr_nrow +1):((nsimblocks + 1)*curr_nrow)), "BetweenCovMultFactor"] <- 2
nsimblocks <- nsimblocks + 1

RSimVals <- rbind(RSimVals, RSimVals[((1):(curr_nrow)),])
RSimVals[((nsimblocks*curr_nrow +1):((nsimblocks + 1)*curr_nrow)), "BetweenCovMultFactor"] <- 4
nsimblocks <- nsimblocks + 1

RSimVals <- rbind(RSimVals, RSimVals[((1):(curr_nrow)),])
RSimVals[((nsimblocks*curr_nrow +1):((nsimblocks + 1)*curr_nrow)), "BetweenCovMultFactor"] <- 10
nsimblocks <- nsimblocks + 1

RSimVals <- rbind(RSimVals, RSimVals[((1):(curr_nrow)),])
RSimVals[((nsimblocks*curr_nrow +1):((nsimblocks + 1)*curr_nrow)), "BetweenCovMultFactor"] <- 20
nsimblocks <- nsimblocks + 1

RSimVals <- rbind(RSimVals, RSimVals[((1):(curr_nrow)),])
RSimVals[((nsimblocks*curr_nrow +1):((nsimblocks + 1)*curr_nrow)), "MeasErrorMultFactor"] <- 0.50
nsimblocks <- nsimblocks + 1

RSimVals <- rbind(RSimVals, RSimVals[((1):(curr_nrow)),])
RSimVals[((nsimblocks*curr_nrow +1):((nsimblocks + 1)*curr_nrow)), "MeasErrorMultFactor"] <- 0.25
nsimblocks <- nsimblocks + 1

RSimVals <- rbind(RSimVals, RSimVals[((1):(curr_nrow)),])
RSimVals[((nsimblocks*curr_nrow +1):((nsimblocks + 1)*curr_nrow)), "MeasErrorMultFactor"] <- 0.10
nsimblocks <- nsimblocks + 1

RSimVals <- rbind(RSimVals, RSimVals[((1):(curr_nrow)),])
RSimVals[((nsimblocks*curr_nrow +1):((nsimblocks + 1)*curr_nrow)), "MeasErrorMultFactor"] <- 0.05
nsimblocks <- nsimblocks + 1

RSimVals <- rbind(RSimVals, RSimVals[((1):(curr_nrow)),])
RSimVals[((nsimblocks*curr_nrow +1):((nsimblocks + 1)*curr_nrow)), "MeasErrorMultFactor"] <- 0.01
nsimblocks <- nsimblocks + 1

RSimVals <- rbind(RSimVals, RSimVals[((1):(curr_nrow)),])
RSimVals[((nsimblocks*curr_nrow +1):((nsimblocks + 1)*curr_nrow)), "MeasErrorMultFactor"] <- 0
nsimblocks <- nsimblocks + 1

#Add simulations for 26 and 50 samples 
  #Add these new nsamp values below below the other simulation values to make sure the results saved for lower number correspond
  #properly to how they were run

#Need to split the RSimsVals because the cluster is unable to accept array values greater than 40000
sub_all_sims <- RSimVals[seq(1, nrow(RSimVals), by = 3),]
sub_all_sims2 <- sub_all_sims
sub_all_sims3 <- sub_all_sims

sub_all_sims$nsamp <- 26
sub_all_sims2$nsamp <- 50
sub_all_sims3$nsamp <- 80

RSimVals <- rbind(RSimVals, sub_all_sims, sub_all_sims2, sub_all_sims3)


#Now, add results for effect sizes c(1,1.375) and c(1,1.75)
  #Want to calculate these results for every scenario mentioned above, so subset to other unique combinations of this
t2 <- subset(RSimVals, RSimVals$condEffSizeStr=="c(1,1)")
t3 <- t2
t2$condEffSizeStr <- "c(1,1.375)"
t3$condEffSizeStr <- "c(1,1.75)"

RSimVals <- rbind(RSimVals, t2, t3)
RSimVals$UniqString <- paste0("nsamp", RSimVals$nsamp, "condEffSizeStr", RSimVals$condEffSizeStr, 
                              "GenWithMeasError", RSimVals$GenWithMeasError, "MeasErrorMultFactor", RSimVals$MeasErrorMultFactor, 
                              "BetweenCovMultFactor", RSimVals$BetweenCovMultFactor)

RSimVals$SimVal <- 1:nrow(RSimVals)
RSimVals$StartSeed <- RSimVals$SimVal * 1e6


save(RSimVals, file = "RSimVals.RData", version = 2)


#curr_nrow <- nrow(RSimVals)
#RSimVals <- rbind(RSimVals, RSimVals)
#RSimVals[((1):(curr_nrow)), "GenWithMeasError"] <- TRUE
