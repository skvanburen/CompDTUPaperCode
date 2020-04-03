#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16g


mkdir -p SaveGibbsAsRLogs
srun R CMD BATCH --vanilla SaveGibbsDataAsRData.R SaveGibbsAsRLogs/SaveGibbsDataAsRData_$SLURM_ARRAY_TASK_ID.Rout
