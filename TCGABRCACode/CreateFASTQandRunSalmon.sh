#!/bin/bash

#SBATCH --time=47:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16g
#SBATCH --output=/nas/longleaf/home/skvanbur/res/TCGABRCAAnalysis/CreateFASTQandRunSalmonLogs/slurmLog_%a.out
#SBATCH --error=/nas/longleaf/home/skvanbur/res/TCGABRCAAnalysis/CreateFASTQandRunSalmonLogs/slurmLog_%a.err

sleep $(( RANDOM % 1000))

#mkdir -p ~/res/TCGABRCAAnalysis/SalmonBootSampsLogs/
mkdir -p ~/res/TCGABRCAAnalysis/CreateFASTQandRunSalmonRLogs/
srun R CMD BATCH --vanilla ~/res/TCGABRCAAnalysis/CreateFASTQandRunSalmon.R ~/res/TCGABRCAAnalysis/CreateFASTQandRunSalmonRLogs/CreateFASTQandRunSalmon_$SLURM_ARRAY_TASK_ID.Rout




