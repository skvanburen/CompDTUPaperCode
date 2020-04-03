#!/bin/bash

#
#
#  Created by Scott Van Buren on 6/6/18.
#

##### SBATCH --time=72:00:00
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=29g

#SBATCH --output=/pine/scr/s/k/skvanbur/GEUV1/updateilrMeansCovsAndabDatasetsslurmLogs/slurmLog_%A_%a.out
sleep $(( RANDOM % 1000))

mkdir -p updateilrMeansCovsAndabDatasetsLogs
srun R CMD BATCH --vanilla updateilrMeansCovsAndabDatasets.R updateilrMeansCovsAndabDatasetsLogs/updateilrMeansCovsAndabDatasetsLogs_$SLURM_ARRAY_TASK_ID.Rout
