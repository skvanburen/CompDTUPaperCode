#!/bin/bash

#  GEUV1DRIMSeqPower.sh
#  
#
#  Created by Scott Van Buren on 7/3/18.
#

### The 20 sample analysis should take less than an hour (so give it mabe 4-6 to be safe)
    ### with maybe 15-16 GB of RAM required
    ### and the 100 sample case should take longer so give is 72 hours just to be safe
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=24g
## SBATCH --output=/pine/scr/s/k/skvanbur/GEUV1/DRIMSeqPowerLogs
## SBATCH --output=/dev/null
## SBATCH --output=/pine/scr/s/k/skvanbur/GEUV1/CompositionalPowerLogs2/slurmLog_%A_%a.out

sleep $(( RANDOM % 100))
mkdir -p DRIMSeqPowerLogs/
#baseloc=/pine/scr/s/k/skvanbur/GEUV1/DRIMSeqPowerLogs
srun R CMD BATCH --vanilla GEUV1DRIMSeqPower.R DRIMSeqPowerLogs/GEUV1DRIMSeqPower_$SLURM_ARRAY_TASK_ID.Rout
