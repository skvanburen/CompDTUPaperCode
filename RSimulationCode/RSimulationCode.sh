#!/bin/bash

###Usually takes < 5 hours, but some values can take much longer such that it is best to give it a much longer time limit
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1g
#SBATCH --output=/pine/scr/s/k/skvanbur/GEUV1/RSimsslurmLogs/slurmLog_%A_%a.out


sleep $(( RANDOM % 100))
mkdir -p /pine/scr/s/k/skvanbur/GEUV1/RSimsLogs
###mkdir -p /pine/scr/s/k/skvanbur/GEUV1/RSimsLogsAdditionalNSamps
srun R CMD BATCH ~/res/GEUV1Data/RSimulationCode/RSimulationCode.R  /pine/scr/s/k/skvanbur/GEUV1/RSimsLogs/RSimulationCode_$SLURM_ARRAY_TASK_ID.Rout
