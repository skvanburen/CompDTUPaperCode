#!/bin/bash

#  JobArrayTrial.sh
#
#
#  Created by Scott Van Buren on 11/8/17.
#
### Keep the memory used higher if generating bootstrap samples, reduce to rerun to generate equivalence class files
### Also keep time higher if generating bootstrap samples, generally takes around 3-6 hours (with sleep statments to reduce file loading errors) if generating bootstrap samples and maybe 1-2 of not

#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=15g

sleep $(( RANDOM % 1000))

#mkdir -p ~/res/GEUV1Data/SalmonLogs/
#srun R CMD BATCH --vanilla ~/res/GEUV1Data/RunSalmonGEUV1Data.R ~/res/GEUV1Data/SalmonLogs/SalmonRunLog_$SLURM_ARRAY_TASK_ID.Rout

mkdir -p ~/res/GEUV1Data/SalmonBootSampsLogs/
srun R CMD BATCH --vanilla ~/res/GEUV1Data/RunSalmonGEUV1Data.R ~/res/GEUV1Data/SalmonBootSampsLogs/SalmonRunLog_$SLURM_ARRAY_TASK_ID.Rout




