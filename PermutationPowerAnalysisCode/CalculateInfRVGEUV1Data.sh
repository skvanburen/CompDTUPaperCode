#!/bin/bash

#  CalculateInfRVGEUV1Data.sh
#  
#
#  Created by Scott Van Buren on 11/7/19.
#
#SBATCH --time=47:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16g

srun R CMD BATCH --vanilla CalculateInfRVGEUV1Data.R CalculateInfRVGEUV1Data_$SLURM_ARRAY_TASK_ID.Rout
