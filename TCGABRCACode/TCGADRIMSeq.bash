#!/bin/bash



#SBATCH --time=71:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32g
#SBATCH --partition=general_big
#SBATCH --nodelist=c0405
#SBATCH --output=/nas/longleaf/home/skvanbur/res/TCGABRCAAnalysis/TCGADRIMSeqLogs/slurmLog_%a.out
#SBATCH --error=/nas/longleaf/home/skvanbur/res/TCGABRCAAnalysis/TCGADRIMSeqLogs/slurmLog_%a.err


mkdir -p ~/res/TCGABRCAAnalysis/TCGADRIMSeqRLogs/
srun R CMD BATCH --vanilla ~/res/TCGABRCAAnalysis/TCGADRIMSeq.R ~/res/TCGABRCAAnalysis/TCGADRIMSeqRLogs/TCGADRIMSeq_$SLURM_ARRAY_TASK_ID.Rout
