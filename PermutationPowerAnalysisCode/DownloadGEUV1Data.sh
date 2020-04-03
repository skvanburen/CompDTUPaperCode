#!/bin/bash

#  JobArrayTrial.sh
#  
#
#  Created by Scott Van Buren on 11/8/17.
#

#SBATCH --time=5:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64m



srun R CMD BATCH --vanilla DownloadGEUV1Data.R DownloadLogs/DownloadGEUV1Data_$SLURM_ARRAY_TASK_ID.Rout

