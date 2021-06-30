#!/bin/bash



#SBATCH --time=71:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=25g
#SBATCH --output=/nas/longleaf/home/skvanbur/res/TCGABRCAAnalysis/TCGASaveNecessaryDatasetsForCompDTURegLogs/slurmLog_%a.out
#SBATCH --error=/nas/longleaf/home/skvanbur/res/TCGABRCAAnalysis/TCGASaveNecessaryDatasetsForCompDTURegLogs/slurmLog_%a.err

sleep $(( RANDOM % 1000))

mkdir -p ~/res/TCGABRCAAnalysis/TCGASaveNecessaryDatasetsForCompDTURegRLogs/
srun R CMD BATCH --vanilla ~/res/TCGABRCAAnalysis/TCGASaveNecessaryDatasetsForCompDTUReg.R ~/res/TCGABRCAAnalysis/TCGASaveNecessaryDatasetsForCompDTURegRLogs/TCGASaveNecessaryDatasetsForCompDTUReg_$SLURM_ARRAY_TASK_ID.Rout
