#!/bin/sh

#  RATsPower.sh
#  
#
#  Created by Scott Van Buren on 4/4/19.
#

####Memory needs to be at 24g for the 100 sample case but can be set to 12g for the 20 sample case
#####Computation time should be set at 48 hours for the 100 sample case and can be set to 12 hours
##### for the 20 sample case but it should actually compute much faster
#SBATCH --time=71:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=28g

mkdir -p RATsPowerLogs

srun R CMD BATCH --vanilla RATsPower.R RATsPowerLogs/RATsPower_$SLURM_ARRAY_TASK_ID.Rout
