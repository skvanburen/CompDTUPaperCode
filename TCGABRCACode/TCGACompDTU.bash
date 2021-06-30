#!/bin/bash



#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32g
#SBATCH --partition=general_big
#SBATCH --nodelist=c0405
#SBATCH --output=/nas/longleaf/home/skvanbur/res/TCGABRCAAnalysis/TCGACompDTULogs/slurmLog_%a.out
#SBATCH --error=/nas/longleaf/home/skvanbur/res/TCGABRCAAnalysis/TCGACompDTULogs/slurmLog_%a.err


mkdir -p ~/res/TCGABRCAAnalysis/TCGACompDTURLogs/
srun ~/bin/R CMD BATCH --vanilla ~/res/TCGABRCAAnalysis/TCGACompDTU.R ~/res/TCGABRCAAnalysis/TCGACompDTURLogs/TCGACompDTU_$SLURM_ARRAY_TASK_ID.Rout
