#!/bin/bash

#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50g

srun R CMD BATCH --vanilla SumToGene.R SumToGene_$SLURM_ARRAY_TASK_ID.Rout
