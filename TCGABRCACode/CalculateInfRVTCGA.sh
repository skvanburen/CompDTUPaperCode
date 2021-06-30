#!/bin/bash

#  CalculateInfRVGEUV1Data.sh
#  
#
#  Created by Scott Van Buren on 11/7/19.
#
#SBATCH --time=71:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=29g

srun R CMD BATCH --vanilla CalculateInfRVTCGA.R CalculateInfRVTCGA.Rout
