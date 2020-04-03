#!/bin/bash

#
#
#  Created by Scott Van Buren on 6/6/18.
#

###Longer time is used if all abDatasets/cntDatasets for the Gibbs/Boot samples are created
##### SBATCH --time=72:00:00
#SBATCH --time=59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8g
#SBATCH --output=/pine/scr/s/k/skvanbur/GEUV1/GenerateInfRepSpecificPartFullinfRepDatFilesslurmLogs/slurmLog_%A_%a.out

sleep $(( RANDOM % 100))

mkdir -p GenerateInfRepSpecificPartFullinfRepDatFilesLogs
srun R CMD BATCH --vanilla GenerateInfRepSpecificPartFullinfRepDatFiles.R GenerateInfRepSpecificPartFullinfRepDatFilesLogs/GenerateFullInfRepDatasets_$SLURM_ARRAY_TASK_ID.Rout
