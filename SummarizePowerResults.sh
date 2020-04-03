#!/bin/bash

#  SummarizePowerResults.sh
#  
#
#  Created by Scott Van Buren on 7/10/18.
#
## # # Takes more on the order of 3-6 hours for DRIMSeqFiltering=TRUE with less than 10GB of RAM also
## In the past this time has been more on the order of 16-20 hours depending on the number of genes kept in the analysis.  And memory needed to be more like 22-24 GB also
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16g
## SBATCH --output=/pine/scr/s/k/skvanbur/GEUV1/DRIMSeqPowerLogs
## SBATCH --output=/dev/null
## SBATCH --output=/pine/scr/s/k/skvanbur/GEUV1/CompositionalPowerLogs2/slurmLog_%A_%a.out

sleep $(( RANDOM % 100))
#mkdir -p /pine/scr/s/k/skvanbur/GEUV1/CompositionalPowerLogs2/

#baseloc=/pine/scr/s/k/skvanbur/GEUV1/DRIMSeqPowerLogs
srun R CMD BATCH --vanilla SummarizePowerResultsNewAug2019.R SummarizePowerResultsNewAug2019_$SLURM_ARRAY_TASK_ID.Rout
