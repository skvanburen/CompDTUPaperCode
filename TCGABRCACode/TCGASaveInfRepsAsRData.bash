#!/bin/bash



#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16g
#SBATCH --output=/nas/longleaf/home/skvanbur/res/TCGABRCAAnalysis/TCGASaveInfRepsAsRDataLogs/slurmLog_%a.out
#SBATCH --error=/nas/longleaf/home/skvanbur/res/TCGABRCAAnalysis/TCGASaveInfRepsAsRDataLogs/slurmLog_%a.err

sleep $(( RANDOM % 1000))

mkdir -p ~/res/TCGABRCAAnalysis/TCGASaveInfRepsAsRDataRLogs/
srun R CMD BATCH --vanilla ~/res/TCGABRCAAnalysis/TCGASaveInfRepsAsRData.R ~/res/TCGABRCAAnalysis/TCGASaveInfRepsAsRDataRLogs/TCGASaveInfRepsAsRData_$SLURM_ARRAY_TASK_ID.Rout
