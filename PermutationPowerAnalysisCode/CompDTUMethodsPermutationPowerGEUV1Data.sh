#!/bin/bash

#  CompDTUMethodsPermutationPower.sh
#  
#
#  Created by Scott Van Buren on 10/11/19.
#  


#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2g
#SBATCH --output=/pine/scr/s/k/skvanbur/GEUV1/ExamineImputeandPermuteApproachslurmLogs/slurmLog_%a.out

sleep $(( RANDOM % 100))


mkdir -p /pine/scr/s/k/skvanbur/GEUV1/ExamineImputeandPermuteApproachLogs
srun R CMD BATCH --vanilla ~/res/GEUV1Data/ExamineImputeandPermuteApproachGEUV1Data.R /pine/scr/s/k/skvanbur/GEUV1/ExamineImputeandPermuteApproachLogs/ExamineImputeandPermuteApproach_$SLURM_ARRAY_TASK_ID.Rout
