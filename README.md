This repository contains code to reproduce results from "Differential Transcript Usage Analysis Incorporating Quantification Uncertainty via Compositional Measurement Error Regression Modeling."  The paper is currently available as a preprint at https://www.biorxiv.org/content/10.1101/2020.05.22.111450v1.  Additionally, the described methods in the preprint can be run using the CompDTUReg package available here https://github.com/skvanburen/CompDTUReg.

Last updated May 27, 2020

Commands 1-19 reproduce the permutation power analysis, while the instructions 20-22 will reproduce the Simulation based power analyses.  


First, a few commands that may be useful for reproducing results on a slurm based cluster:

Open an interactive session of R/bash:<br/>
````module load r/3.6.0````<br/>
````srun -p interact -N 1 -n 1 --mem=32g --time=8:00:00 --pty R --no-init-file````<br/>

````srun -p interact -N 1 -n 1 --time=8:00:00 --mem=4g --pty bash````<br/>


Now, the main instructions to reproduce the results from the paper:

1. First, need to download the E-GEUV-1 data using the following lines, which end up running the R file DownloadGEUV1Data.R .  This array job has 924 parts because each job downloads one paired-end file and there are 462 samples (462*2=924)<br/>



````module load r/3.6.0````<br/>
````sbatch --array=1-924 DownloadGEUV1Data.sh````<br/>



2. Then, generate the Salmon index that will be used in step 3 to quantify the data. This index is generated from gencodeV27 (for the reference transcripts only).  The necessary GENCODE files can be downloaded from https://www.gencodegenes.org/human/release_27.html

````sbatch -N 1 -n 1 --mem=32g --time=2:00:00 -o /pine/scr/s/k/skvanbur/GEUV1/Salmon/indexlog.txt --wrap "~/bin/Salmon0.11.3/bin/salmon index -p 1 -t gencode.v27.transcripts.fa.gz -i /pine/scr/s/k/skvanbur/GEUV1/SalmonBootSamps/transcripts_index --gencode"````




3. Next, run the Salmon quantifications themselves on all 462 samples using the following lines.  This code runs the R file RunSalmonGEUV1Data.R
 


````module load r/3.6.0````<br/>
````sbatch --array=1-462 RunSalmonGEUV1Data.sh````<br/>



4. Next, create the the tx2gene data frame in r that matches transcripts to genes.  This will be needed by various downstream code. This code will not take too long to run, so can be run interactively or like this:

````module load r/3.6.0````<br/>
````sbatch -N 1 -n 1 --mem=8g --time=2:00:00 -o /pine/scr/s/k/skvanbur/GEUV1/MakeTx2Gene.out --wrap “Rscript MakeTx2Gene.R"````<br/> 




5. Next, save the Salmon output (not including bootstrap/Gibbs samples for now) and the key file containing the population membership of each sample to a .RData file for later use.  This file should not take long to run so can probably run interactively (using the command at the top of the file) or like this:

````module load r/3.6.0````<br/>
````sbatch -N 1 -n 1 --mem=32g --time=2:00:00 -o /pine/scr/s/k/skvanbur/GEUV1/SaveSalmonDataAsRData.out --wrap “Rscript GEUV1SaveSalmonDataAsRData.R"````<br/>



6. Next, run “SumToGene.R” to summarize results to the gene level, both for counts and TPM values.  This file saves the data in two structures.  The abGene and cntGene files respectively contain transcripts on the rows and sample-specific expression values on the columns as well as average expression across gene/conditions, information on the highest expressed transcript across samples (major transcript), etc.  The abGene file contains TPM measurements, and the cntGene file contains counts and effective length measurements for each transcript/sample.  The abDatasets and cntDatasets files contain lists of data frames, with each data frame being reformatted data for each gene.  Specifically, in these data frames samples are on the rows and the columns contain transcript-level expression for a given sample.  The list of transcripts in these files will be all transcripts in the gene for any files saved with “NoOtherGroups” in the name and  transcripts collapsed into other groups for files names “abDatasets.RData” or “cntDatasets.RData”.  These other groups combine and transcript into an “Other” category if the transcript has less than 5% of the overall gene-level expression.  We had originally used this as out filtering approach but switched to using DRIMSeq’s filtering procedure for all results the paper.  DRIMseq’s filtering procedure is implemented in the file “GEUV1FilterGenesForAnalysis.R”  below.  Make sure to set “OnlyCalculateAbundanceFromInfReps”  in the R file to FALSE initially. Additionally, the file can be submitted as an array job using the second command below to save the results for abundance data that is generated by taking the mean or median of bootstrap samples (using the option infRepStat in tximport) once it has been saved for the usual abundance estimates using the command above it.  To do this, make sure to set OnlyCalculateAbundanceFromInfReps to TRUE in the R file.



````module load r/3.6.0````<br/>
````sbatch -N 1 -n 1 --mem 16g --time=2:00:00 -o ~/res/GEUV1Data/SumToGene.txt --wrap "Rscript SumToGene.R"````<br/>



````module load r/3.6.0````<br/>
````sbatch --array=1-4 SumToGene.sh````<br/>


7. Now, filter the output of “SumToGene.R” from above using DRIMSeq’s filtering procedure.  First, run the top set of commands so save the filtered results for the regular abundance estimates.  Then, if desired, this job can be submitted as an array job using the second command below to save the filtered datasets for abundance data that is generated by taking the mean or median of bootstrap samples (using the option infRepStat in tximport) once it has been saved for the usual abundance estimates using the command above it.  This needs to be first completed for the usual abundance estimates because the list of transcripts that passes from the usual estimates is used as the set that passes for the abundance estimates calculated using the mean/median of bootstrap/gibbs samples to make comparisons between the different datasets as smooth as possible.



````module load r/3.6.0````<br/>
````sbatch -N 1 -n 1 --mem 16g --time=4:00:00 -o ~/res/GEUV1Data/GEUV1FilterGenesForAnalysis.txt --wrap "Rscript GEUV1FilterGenesForAnalysis.R"````<br/>



````module load r/3.6.0````<br/>
````sbatch --array=1-4 GEUV1FilterGenesForAnalysis.sh````<br/>



8. Now, save the inferential replicates as separate .RData files for each biological sample.  Results need to be saved separately for each biological sample because otherwise the files get too large.  While the file is called “SaveGibbsDataAsRData.R”, it saves data for either bootstrap or Gibbs inferential replicates.



````module load r/3.6.0````<br/>
````sbatch --array=1-462 SaveGibbsDataAsRData.sh````<br/>


9. Now, save full inferential replicate data files.  These files reorganize the files saved by “SaveGibbsDataAsRData.R” to contain all inferential replicate values for all samples for all transcripts from a specific subset of genes.  To get this subset of genes, the full gene list that remains after filtering using DRIMSeq’s procedure is completed in the file “GEUV1FilterGenesForAnalysis.R” is split up into 100 parts such that all 100 files taken together contain information for all genes that pass filtering. Results must be split up this way because otherwise the files are too large.



````module load r/3.6.0````<br/>
````sbatch --array=1-100 GenerateFullInfRepDatasets.sh````<br/>



10. Now, once the previous step has completed also save files that split results to be from a specific inferential replicate number for a specific chunk of the genes.  This will speed up code that modifies the datasets for the power analysis.



````module load r/3.6.0````<br/>
````sbatch --array=1-100 GenerateInfRepSpecificPartFullinfRepDatFiles.sh````<br/>



11. Now, generate the files that are needed to run the power analysis on the data, including the group assignments and a necessary key file, etc with the file “GenerateValuesForPowerAnalysis1.R”.  This file will not take long to run so can be run interactively or like this:

````module load r/3.6.0````<br/>
````sbatch -N 1 -n 1 --mem=8g --time=2:00:00 -o /pine/scr/s/k/skvanbur/GEUV1/GenerateValuesForPowerAnalysis1.out --wrap “Rscript GenerateValuesForPowerAnalysis1.R”````<br/> 

12. Now, update the datasets in accordance with the procedure described in the paper.  Each array value will run one of 100 parts (chunks) of genes.  Additionally, different change values (comparable to effect sizes) are run using different values.  Specifically, for the 100 total samples power analysis, array values 1-100 correspond to a change value of 4, 101-200 to a value of 1, 201-300 to a value of 2, and 301-400 to a value of 8.  Additionally, adding 1000 to the values generates the datasets for the 20 total sample power analysis (such that now array values 1001-1100 correspond to a change value of 4, 1101-1200 to a value of 1, 1201-1300 to a value of 2, and 1301-1400 to a value of 8). The code is structured in such a way that if the code does not complete, the job can be resubmitted and it will leave off at whichever file it was at before instead of starting from the beginning and overwriting files that are already there.



````module load r/3.6.0````<br/>
````sbatch --array=1-300,1001-1400 updateilrMeansCovsAndabDatasets.sh````<br/>



13. Now, run the power analyses for CompDTU and CompDTUme.  These results are executed and output by a file called “CompDTUMethodsPermutationPowerGEUV1Data.R”.  This code runs the CompDTU and CompDTUme methods in addition to the CompMI approaches as well.  Setting TwentySamplesTotalAnalysis at the top of the file to be TRUE will run the results for the 20 sample total analysis, and setting it to FALSE will run the results for the 100 total sample analysis.   Lastly, the multiple imputation based results in the supplement are run from this file by setting RunCompMIResults to be TRUE.  Setting this to TRUE will cause it to run the MI methods only and not run the CompDTU and CompDTUme results.  As a result, the MI values will run very quickly (usually less than 10-20 minutes for each array value). If you are only interested in the CompDTU and CompDTUme results you can and should set the option calcCompDTUandCompDTUmeOnly=TRUE in the function call.



````module load r/3.6.0````<br/>
````sbatch --array=1-30000 CompDTUMethodsPermutationPowerGEUV1Data.sh````<br/>



14.  Now, run the DRIMSeq power analyses.  The default results presented in the paper add 1 to every count to greatly stabilize the DRIMSeq power results, though this option can be changed and theoretically shouldn’t be needed, especially if using the add_uniform=TRUE option.  Still, we found in many cases adding 1 to every count was needed to stabilize the results.



````module load r/3.6.0````<br/>
````sbatch --array=1-300 GEUV1DRIMSeqPower.sh````<br/>



15. Now, run the RATs power analyses:



````module load r/3.6.0````<br/>
````sbatch --array=1-300 RATsPower.sh````<br/>

16. Now, calculate the inferential variability of the genes on the ilr scale, which is used to group genes and results post hoc (but not by the method to determine significance directly).  These files need to be saved for the SummarizePowerResults.R results to work (otherwise, you can set useInfRV to FALSE in that file).  Set array_val equal to 2 for bootstrap samples (the default used in the paper) and 1 for Gibbs samples.



````module load r/3.6.0````<br/>
````sbatch --array=2 CalculateInfRVGEUV1Data.sh````<br/>



17. Now, run the file that calculates the overlap values for each gene.  This is done via the file “CalcTransOverlap.R”, which runs quickly such that it can be run interactively.




18. Now, run the file that summarizes the power results for the various methods and calculate sens, spec, etc for all genes and various subsets of genes that are of interest.  The array values control the current method that is being run as well as the current change value being evaluated (ie values < 100 are for change value of 1, 101-200 are a change value of 2, etc).


````module load r/3.6.0````<br/>
````sbatch --array=1-9,101-109,201-209,301-309,1001-1009,1101-1109,1201-1209 SummarizePowerResults.sh````<br/>

19. Once the results in the previous step have been completed, download the resulting files (from the “GEUV1PowerResBootDRIMSeqFiltering” folder since use of bootstraps and DRIMSeq’s filtering procedure is the default) to your computer and generate the various ROC and other curves of interest using the file “GenerateROCCurves.R”.  Within this file you can change the change value and plotType being plotted and the methods being used/their colors, etc.  Additionally, note that changing SuppPlots to TRUE will generate plots including methods presented in the supplement (ie including the RATs methods) and setting to FALSE will only include methods included in the main paper.


Now, instructions to reproduce the simulation power analysis.  The original version of this code is saved within the RSimulationCode directory within this folder such that no additional syncing of files from other locations is necessary.

20. To run the simulation code, first run the file “ConstructRSimulationValues.R”.  This file specifies the parameter values for each simulation and saves these values in the data file RSimVals.RData.  This file can be run locally or interactively and takes very little time to run.

21. Once the file above has completed, run the file RSimulationCode.sh (which in turn calls RSimulationCode.R) within the directory.   This needs to be submitted as an array job, where each simulation setup is split up into parts to enable the code to run in parallel across different jobs.  Each job is currently split up such that each “part” runs 1000 simulation runs, and each setup is currently run 10 times to give 10,000 total runs per setup.  Thus, each simulation scenario has 10 parts of results.  For example array values 1-10 correspond to simulation scenario 1 parts 1-10, values 101-110 correspond to scenario 2 parts 1-10, 201-210 correspond to scenario 3 parts 1-10, etc.  Thus values ending on 11-99 are not currently used and would be reserved to submit more parts for given scenarios.  The total number of simulation scenarios in the RSimVals.RData file is currently 714, which would require array value numbers greater than the maximum of about 40000 to be submitted.  To fix this, change the option RunAdditionalNSamps to TRUE and submit using array values corresponding to simulations 1-357 such that now values 1-10 correspond to simulation scenario 358 parts 1-10, values 101-110 correspond to scenario 359 parts 1-10, 201-210 correspond to scenario 360 parts 1-10, etc.  To submit these jobs to the cluster use code like the following, making sure to modify options such as RunAdditionalNSamps and and CalculateImputePermutePvals as necessary. 

````cd ~/res/GEUV1Data/RSimulationCode````<br/>
````module load r/3.6.0````<br/>
````for i in `seq 1 357`;````<br/>
````do````<br/>
````let k=(100*($i-1))+1````<br/> 
````let l=(100*($i-1))+10````<br/>
````sbatch --array=$k-$l RSimulationCode.sh````<br/> 
````done````<br/>


22. The computation time table results can be reproduced with the scripts contained within the ComputationTimeTable subdirectory.  Lastly, power table results are produced using the file “PowerTableCodeMar2020.R”.


