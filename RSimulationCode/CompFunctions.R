#Functions for Analysis corresponding to the CompDTU paper
#Last updated 4-3-20


#maketx2gene requires GenomicFeatures from Bioconductor the annotation to be saved as a .gtf.gz file and specified as an input
maketx2gene <- function(txdb_loc, save_loc = NULL, savetsv = FALSE){
  library(GenomicFeatures)
  temp1 <- makeTxDbFromGFF(txdb_loc)
  temp2 <- DataFrame(transcripts(temp1, columns = c("tx_id", "tx_name", "gene_id")))
  tx2genetemp <- data.frame(as.character(temp2$tx_name),as.character(temp2$gene_id), stringsAsFactors = FALSE)
  colnames(tx2genetemp) <- c("tx_id", "gene_id")

  #Get number of transcripts per gene
  numtranspergene <- data.frame(table(tx2genetemp$gene_id), stringsAsFactors = FALSE)
  colnames(numtranspergene) <- c("gene_id", "NTrans")

  tx2gene <- merge(tx2genetemp, numtranspergene, by = "gene_id")


  #Save transcript to gene file for potential use later (in a possibly specified save_loc location)
  if(is.null(save_loc)){
    save(tx2gene, file = "tx2gene.RData")
  }else{
    curr_wd <- getwd()
    setwd(save_loc)
    save(tx2gene, file = "tx2gene.RData")
    setwd(curr_wd)
  }
  
  if(savetsv==TRUE){
    tx2gene2 <- tx2gene[,c("tx_id", "gene_id")]
    rownames(tx2gene2) <- NULL
    colnames(tx2gene2) <- NULL
    
    write.table(tx2gene2, file='tx2gene.tsv', quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
  }

}


#Outer function called by SumToGene.R that saves data output by SaveSalmonDatatoRData (QuantSalmon) into formats needed to run analyses
sumToGene <- function(QuantSalmon, key, tx2gene, clust, countsFromAbundance, GenAllGroupCombos = FALSE, abFromInfRepFunc = NULL, GibbsSamps = NULL){
  nsamp <- ncol(QuantSalmon$abundance)
  Grouptemp <- key$Condition
  Group <- relevel(as.factor(Grouptemp), ref = 1)

  countsFromAbundance <- QuantSalmon$countsFromAbundance

  abundance <- data.frame(QuantSalmon$abundance)
  #abundance$tx_id <- rownames(abundance)

  counts <- data.frame(QuantSalmon$counts)
  #counts$tx_id <- rownames(counts)

  lengths <- data.frame(QuantSalmon$length)
  #lengths$tx_id <- rownames(lengths)

  ObsData <- sumToGeneHelper(abundance = abundance, counts = counts, lengths = lengths, tx2gene = tx2gene,
                             Group = Group, clust = clust, nsamp = nsamp, key = key,
                             GenAllGroupCombos = GenAllGroupCombos, useExistingMajorTrans = FALSE)

  #Use abundances (ie TPMs) instead of the counts to generate the compositions for compositional analysis
  abGene <- ObsData$abGene
  cntGene <- ObsData$cntGene

  abDatasets <- ObsData$abDatasets
  cntDatasets <- ObsData$cntDatasets

  AllGroupCombinations <- ObsData$AllGroupCombinations

  fullgenenames <- ObsData$fullgenenames

  nonfullrankab <- ObsData$nonfullrankab
  nonfullrankabnames <- ObsData$nonfullrankabnames
  nonfullrankcnt <- ObsData$nonfullrankcnt
  nonfullrankcntnames <- ObsData$nonfullrankcntnames

  abDatasetsCompTime <- ObsData$abDatasetsCompTime
  cntDatasetsCompTime <- ObsData$cntDatasetsCompTime
  abGenecntGeneCompTime <- ObsData$abGenecntGeneCompTime

  ncomb <- nrow(AllGroupCombinations)
  
  if(is.null(abFromInfRepFunc)){
    fil_mod <- ""
  }else if(abFromInfRepFunc=="median" & GibbsSamps==TRUE){
    fil_mod <- "AbRowMedianInfRepGibbs"
  }else if(abFromInfRepFunc=="mean" & GibbsSamps==TRUE){
    fil_mod <- "AbRowMeanInfRepGibbs"
  }else if(abFromInfRepFunc=="median" & GibbsSamps==FALSE){
    fil_mod <- "AbRowMedianInfRepBoot"
  }else if(abFromInfRepFunc=="mean" & GibbsSamps==FALSE){
    fil_mod <- "AbRowMeanInfRepBoot"
  }else{
    stop("Check specification of abFromInfRepFunc function")
  }

  if(countsFromAbundance =="no"){
    cntGenefil <- paste0("cntGene", fil_mod, ".RData")
    cntDatasetsfil <- paste0("cntDatasets", fil_mod, ".RData")
    cntDatasetsNoOtherfil <- paste0("cntDatasetsNoOtherGroups", fil_mod, ".RData")
  }else if(countsFromAbundance =="scaledTPM"){
    cntGenefil <- paste0("cntGenecntsScaledTPM", fil_mod, ".RData")
    cntDatasetsfil <- paste0("cntDatasetscntsScaledTPM", fil_mod, ".RData")
    cntDatasetsNoOtherfil <- paste0("cntDatasetsNoOtherGroupscntsScaledTPM", fil_mod, ".RData")
  }
  #Save abundance (abGene) and counts (cntGene) with nsamp information for use downstream
  save(abGene, nsamp, key, abGenecntGeneCompTime, file = paste0("abGene", fil_mod, ".RData"))
  save(cntGene, nsamp, countsFromAbundance, key, abGenecntGeneCompTime, file = cntGenefil)

  save(abDatasets, nsamp, key, fullgenenames, Group, nonfullrankab, nonfullrankabnames, abDatasetsCompTime, file = paste0("abDatasets", fil_mod, ".RData"))

  save(cntDatasets, nsamp, key, fullgenenames, countsFromAbundance, Group, nonfullrankcnt,
       nonfullrankcntnames, cntDatasetsCompTime, file = cntDatasetsfil)
  save(AllGroupCombinations, ncomb, file = "AllGroupCombinations.RData")

  rm(ObsData)
  rm(abGene)
  rm(cntGene)

  #rm(abDatasets)
  rm(cntDatasets)
  rm(AllGroupCombinations)


  #Now, generate the abDatasets without using OtherGroups and save results
  ObsDataNoOtherGroups <- sumToGeneHelper(abundance = abundance, counts = counts, lengths = lengths, tx2gene = tx2gene,
                                          Group = Group, clust = clust, nsamp = nsamp, key = key,
                                          abCompDatasets = abDatasets, useOtherGroups = FALSE, useExistingMajorTrans = TRUE)

  abDatasetsNoOtherGroups <- ObsDataNoOtherGroups$abDatasets
  cntDatasetsNoOtherGroups <- ObsDataNoOtherGroups$cntDatasets

  abDatasetsNoOtherGroupsCompTime <- ObsDataNoOtherGroups$abDatasetsCompTime
  cntDatasetsNoOtherGroupsCompTime <- ObsDataNoOtherGroups$cntDatasetsCompTime

  #Do not need to resave cntGene or abGene when turning other groups off because those do
  #not use other groups at all and won't change


  fullgenenames <- ObsDataNoOtherGroups$fullgenenames

  nonfullrankabNoOtherGroups <- ObsDataNoOtherGroups$nonfullrankab
  nonfullrankabnamesNoOtherGroups <- ObsDataNoOtherGroups$nonfullrankabnames
  nonfullrankcntNoOtherGroups <- ObsDataNoOtherGroups$nonfullrankcnt
  nonfullrankcntnamesNoOtherGroups <- ObsDataNoOtherGroups$nonfullrankcntnames

  save(abDatasetsNoOtherGroups, nsamp, key, fullgenenames, Group, nonfullrankabNoOtherGroups, nonfullrankabnamesNoOtherGroups, abDatasetsNoOtherGroupsCompTime, file = paste0("abDatasetsNoOtherGroups", fil_mod, ".RData"))
  save(cntDatasetsNoOtherGroups, nsamp, key, fullgenenames, countsFromAbundance, Group, nonfullrankcntNoOtherGroups, nonfullrankcntnamesNoOtherGroups, cntDatasetsNoOtherGroupsCompTime, file = cntDatasetsNoOtherfil)
}

sumToGeneHelper <- function(abundance, counts, lengths, tx2gene, Group, clust, nsamp, key = NULL, abCompDatasets = NULL,
                      useExistingOtherGroups = FALSE,  useOtherGroups = TRUE, useExistingMajorTrans = TRUE,
                      useExistingGenes = FALSE, GenAllGroupCombos = FALSE){

  if(sum(colnames(abundance) != paste0("Sample", 1:nsamp)) != 0){
    stop("Columns must be ordered Sample1, Sample2, etc and tx_ids must just be row names and not a column")
  }

  # if(is.null(clust)){
  #   clust <- makeCluster(1)
  # }
  #Get data into initial format needed

  ST1 <- proc.time()
  initialData <- prepareData(abundance = abundance, counts = counts, lengths = lengths, tx2gene = tx2gene,
                             nsamp = nsamp, key = key)

  abGenecntGeneCompTimeP1 <- proc.time() - ST1
  abGeneTempF <- initialData$abGeneTempF
  cntGeneTempF <- initialData$cntGeneTempF


  ##################################################################################################################
  #Generate list of data frames to be used in compositional analysis
  ##################################################################################################################

  #Cant use if the gene only has 1 trans (since no isoform switching/differential splicing to test for otherwise)
    #or if gene expression across all samples is 0 for all samples
  #So remove genes with only 1 total transcript or ones with no expression in any trans/sample combination
  CompAbGene <- subset(abGeneTempF, (abGeneTempF$NTrans!=1 & abGeneTempF$SumTGE!=0))
  CompCntGene <- subset(cntGeneTempF, (cntGeneTempF$NTrans!=1 & cntGeneTempF$SumTGE!=0))

  #Gene names, make sure to sort this so the right transcript go with the right genes down stream
  #Genenames should be the same using abundance or count data, confirm this with the line below (should be 0)
  if(sum(sort(unique(CompAbGene$gene_id))!=sort(unique(CompCntGene$gene_id))) !=0){
    stop("Something is wrong the the genenames, they should be the same between count and TPM data but are not")
  }
  
  if(useExistingGenes==TRUE){
    fullgenenames <- names(abCompDatasets)
  }else{
    fullgenenames <- sort(unique(CompAbGene$gene_id))
  }

  genestouse <- fullgenenames

  if(useExistingOtherGroups==TRUE | useExistingMajorTrans==TRUE){
    abD <- abCompDatasets
  }else{
    abD <- NULL
  }

  ST2 <- proc.time()
  # abDatasets <- parLapply(clust, genestouse, generateData, dat = CompAbGene,
  #                         nsamp = length(Group), abundance = TRUE, abData = CompAbGene, abCompDatasets = abD,
  #                         useExistingOtherGroups = useExistingOtherGroups, useOtherGroups = useOtherGroups,
  #                         useExistingMajorTrans = useExistingMajorTrans)
  
  abDatasets <- lapply(genestouse, generateData, dat = CompAbGene,
                          nsamp = length(Group), abundance = TRUE, abData = CompAbGene, abCompDatasets = abD,
                          useExistingOtherGroups = useExistingOtherGroups, useOtherGroups = useOtherGroups,
                          useExistingMajorTrans = useExistingMajorTrans)
  # library(plyr)
  #
  # abDatasets <- laply(genestouse, generateData, dat = CompAbGene,
  #                         nsamp = length(Group), abundance = TRUE, abData = CompAbGene, abCompDatasets = abD,
  #                         useExistingOtherGroups = useExistingOtherGroups, useOtherGroups = useOtherGroups,
  #                         useExistingMajorTrans = useExistingMajorTrans, .progress = "text", .inform = TRUE)

  names(abDatasets) <- genestouse
  abDatasetsCompTime <- proc.time() - ST2

  #Use the other groups from abDatasets created just above
  ST3 <- proc.time()
  # cntDatasets <- parLapply(clust, genestouse, generateData, dat = CompCntGene,
  #                          nsamp = length(Group), abundance = FALSE, abData = CompAbGene,
  #                          abCompDatasets = abDatasets, useExistingOtherGroups = TRUE,
  #                          useOtherGroups = useOtherGroups, useExistingMajorTrans = TRUE)
  
  
  cntDatasets <- lapply(genestouse, generateData, dat = CompCntGene,
                           nsamp = length(Group), abundance = FALSE, abData = CompAbGene,
                           abCompDatasets = abDatasets, useExistingOtherGroups = TRUE,
                           useOtherGroups = useOtherGroups, useExistingMajorTrans = TRUE)

  # cntDatasets <- laply(genestouse, generateData, dat = CompCntGene,
  #                          nsamp = length(Group), abundance = FALSE, abData = CompAbGene,
  #                          abCompDatasets = abDatasets, useExistingOtherGroups = TRUE,
  #                          useOtherGroups = useOtherGroups, useExistingMajorTrans = TRUE,
  #                          .progress = "text", .inform = TRUE)
  names(cntDatasets) <- genestouse
  cntDatasetsCompTime <- proc.time() - ST3
  #Use this laply loop to help with debugging (will say which gene the code failed at)




  #Add the major transcript to the abundance and count dataframes abGene/cntGene
  #Do this before filtering out genes to get abDatasets and cntDatasets since you could always
    #filter those out later if needed
    #Major trans is the trans witin a gene that has highest average TPM measurement across al samples, even for the count data
    #This is to ensure that the major trans for a gene is the same between the TPM and count measurements
    #This becomes relevant for the power analyses, when the counts of the major transcripts aer modified
  ST4 <- proc.time()
  Temp1 <- addMajorTrans(genestouse = genestouse, abGeneTempF = abGeneTempF, cntGeneTempF = cntGeneTempF, abDatasets = abDatasets)
  abGene <- Temp1$abGene
  cntGene <- Temp1$cntGene
  abGenecntGeneCompTimeP2 <- proc.time() - ST4

  abGenecntGeneCompTime <- abGenecntGeneCompTimeP1 + abGenecntGeneCompTimeP2

  #Generate all possible group arrangements to be able to easily rerun code on each group arrangement
    #This is only possible if there aren't too many samples and is only working for the 10 SQCCData samples for now
    #So by default this is turned of
  if(GenAllGroupCombos==TRUE){
    library(gtools)
    nsamp <- length(Group)
    numCond1 <- sum(Group==levels(Group)[1])

    #Choose elements corresponding to level 1
    Combs <- combinations(nsamp, numCond1)
    ncomb <- nrow(Combs)
    #Construct all complete Group combinations, each row is a possible arrangement of group
    AllGroupCombinations <- matrix(NA, nrow = ncomb, ncol = nsamp)
    for (i in 1:ncomb){
      for (j in Combs[i,]){
        AllGroupCombinations[i,j] <- "A"
      }
      for(k in 1:nsamp)
        if(is.na(AllGroupCombinations[i,k])){
          AllGroupCombinations[i,k] <- "B"
        }
    }
  }else{
    AllGroupCombinations <- NULL
  }

  #Keep track of how many genes within abDatasets/cntDatasets are not full rank/ which genes specifically are not full rank
  nonfullrankab <- 0
  nonfullrankabnames <- c()

  nonfullrankcnt <- 0
  nonfullrankcntnames <- c()

  for (i in 1:length(genestouse)){
    if(is.null(attr(abDatasets[[i]], "FullRank"))){
      next
    }
    if(attr(abDatasets[[i]], "FullRank")==FALSE){
      nonfullrankab <- nonfullrankab + 1
      nonfullrankabnames <- c(nonfullrankabnames, names(abDatasets)[i])
    }

    #If there is an error, I think it is at this line
    if(attr(cntDatasets[[i]]$Counts, "FullRank")==FALSE){
      nonfullrankcnt <- nonfullrankcnt + 1
      nonfullrankcntnames <- c(nonfullrankcntnames, names(cntDatasets)[i])
    }
  }


  #Quick check to make sure the majortrans are the same between abundances or counts
  ndiff <- 0
  diff_i <- c()
  for (i in 1:length(genestouse)){
    if(is.null(attr(abDatasets[[i]], "MajorTrans"))){
      next
    }
    if(attr(abDatasets[[i]], "MajorTrans")!=attr(cntDatasets[[i]]$Counts, "MajorTrans")){
      ndiff <- ndiff + 1
      diff_i <- c(diff_i, i)

    }
  }
  if(ndiff!=0){
    stop("Something is wrong with the major trans, it should be the same between TPMs and counts but is not.")
  }

  return(list(abGene = abGene, cntGene = cntGene, abDatasets = abDatasets,
              cntDatasets = cntDatasets, AllGroupCombinations = AllGroupCombinations,
              fullgenenames = fullgenenames, nonfullrankab = nonfullrankab,
              nonfullrankabnames = nonfullrankabnames, nonfullrankcnt = nonfullrankcnt,
              nonfullrankcntnames = nonfullrankcntnames, abDatasetsCompTime = abDatasetsCompTime,
              cntDatasetsCompTime = cntDatasetsCompTime, abGenecntGeneCompTime = abGenecntGeneCompTime))

}


#Prepare raw count/TPM data for use with generateData and later functions
#abundance, counts, and lengths should be data.frames with rows being transcript level values and
#columns corresponding to each sample expression (names should only be Sample1, Sample2, etc not Sample1TPM, etc)
#Samps argument gives sample names
#Only give columns tx_id, Sample1Cnt(or TPM) (or just Sample1, Sample2), etc - nothing else including gene_id
#If key is not specified, the sum of the total gene expression (sumTGE) cannot be calculated for each condition (neither can meanTGE)
prepareData <- function(abundance, counts, lengths, tx2gene, nsamp, key = NULL, infReps = "none", samps = NULL) {
  abundance <- data.frame(abundance)

  if(is.null(abundance$tx_id)){
    abundance$tx_id <- rownames(abundance)
  }

  counts <- data.frame(counts)
  if(is.null(counts$tx_id)){
    counts$tx_id <- rownames(counts)
  }


  lengths <- data.frame(lengths)
  if(is.null(lengths$tx_id)){
    lengths$tx_id <- rownames(lengths)
  }

  if(is.null(abundance$NTrans)){
    abGeneTemp1 <- merge(tx2gene, abundance, by = "tx_id")
  }else{
    abGeneTemp1 <- abundance
  }


  if(is.null(counts$NTrans)){
    cntGeneTemp1 <- merge(tx2gene, counts, by = "tx_id")
  }else{
    cntGeneTemp1 <- counts
  }

  abGeneTemp2 <- as.data.frame(abGeneTemp1[order(abGeneTemp1$gene_id, abGeneTemp1$tx_id),])
  cntGeneTemp2 <- as.data.frame(cntGeneTemp1[order(cntGeneTemp1$gene_id, cntGeneTemp1$tx_id),])

  #t3 <- subset(abGeneTemp2, abGeneTemp2$tx_id=="ENST00000161559.10")
  if(infReps == "none"){
    rownames(abGeneTemp2) <- abGeneTemp2$tx_id
    rownames(cntGeneTemp2) <- cntGeneTemp2$tx_id
  }



  #column numbers of the TPM/count info for all samples
  if(is.null(samps)==FALSE){
    abcolnums <- which(colnames(abGeneTemp2) %in% samps)
    cntcolnums <- which(colnames(cntGeneTemp2) %in% samps)

    if(length(abcolnums)==0){
      abcolnums <- which(colnames(abGeneTemp2) %in% paste0(samps, "TPM"))
    }

    if(length(cntcolnums)==0){
      cntcolnums <- which(colnames(cntGeneTemp2) %in% paste0(samps, "Cnt"))
    }

    colnames(abGeneTemp2)[abcolnums] <- paste0(samps, "TPM")
    colnames(cntGeneTemp2)[cntcolnums] <- paste0(samps, "Cnt")
  }else{
    abcolnums <- which(colnames(abGeneTemp2) %in% paste0("Sample", 1:nsamp))
    cntcolnums <- which(colnames(cntGeneTemp2) %in% paste0("Sample", 1:nsamp))

    if(length(abcolnums)==0){
      abcolnums <- which(colnames(abGeneTemp2) %in% paste0("Sample", 1:nsamp, "TPM"))
    }

    if(length(cntcolnums)==0){
      cntcolnums <- which(colnames(cntGeneTemp2) %in% paste0("Sample", 1:nsamp, "Cnt"))
    }

    colnames(abGeneTemp2)[abcolnums] <- paste0("Sample", 1:nsamp, "TPM")
    colnames(cntGeneTemp2)[cntcolnums] <- paste0("Sample", 1:nsamp, "Cnt")
  }


  #Generate total sums for each sample/gene combo by
    #aggregating over each gene (ie total expression for that gene in that sample)
  abgenetotals <- aggregate(abGeneTemp2[,abcolnums], by = list(gene_id = abGeneTemp2$gene_id), FUN = "sum", drop = FALSE)
  cntgenetotals <- aggregate(cntGeneTemp2[,cntcolnums], by = list(gene_id = cntGeneTemp2$gene_id), FUN = "sum", drop = FALSE)

  #TGE short for total gene expression, ie total expression of that gene for that sample across all transcripts

  if(is.null(samps)==FALSE){
    colnames(abgenetotals) <- c("gene_id", paste0(samps, "TGE"))
    colnames(cntgenetotals) <- c("gene_id", paste0(samps, "TGE"))
  }else{
    colnames(abgenetotals) <- c("gene_id", paste0("Sample", 1:nsamp, "TGE"))
    colnames(cntgenetotals) <- c("gene_id", paste0("Sample", 1:nsamp, "TGE"))
  }


  #Calculate total gene expression across all samples because a gene will have to be dropped from analysis if its expression
  # is 0 for all samples - ie, will later filter out genes with SumTGE=0
  abgenetotals$SumTGE <- rowSums(abgenetotals[,-1], na.rm = TRUE)
  cntgenetotals$SumTGE <- rowSums(cntgenetotals[,-1], na.rm = TRUE)

  #Also calculate MeanTGE level (which is juse mean of SumTGE above)
  abgenetotals$MeanTGE <- abgenetotals$SumTGE/nsamp
  cntgenetotals$MeanTGE <- cntgenetotals$SumTGE/nsamp

  #Also add in information on condition specific TGE
  #If key is not specified, the sum of the total gene expression (sumTGE) cannot be calculated for each condition only (neither can meanTGE)
  if(!is.null(key)){
    ncond <- length(unique(key$Condition))
    Group <- relevel(as.factor(key$Condition), ref = 1)
    for(i in 1:ncond){
      assign(paste0("SampsCond", i), subset(key$Identifier, Group==levels(key$Condition)[i]))
      Samps <- get(paste0("SampsCond", i))
      cols <- paste0(Samps, "TGE")
      condSumAb <- rowSums(abgenetotals[,cols])
      condSumCnt <- rowSums(cntgenetotals[,cols])
      val1 <- paste0("SumTGECond", i)
      abgenetotals[,val1] <- condSumAb
      cntgenetotals[,val1] <- condSumCnt

      val2 <- paste0("MeanTGECond", i)
      abgenetotals[,val2] <- condSumAb/length(get(paste0("SampsCond", i)))
      cntgenetotals[,val2] <- condSumCnt/length(get(paste0("SampsCond", i)))

      cols2 <- paste0(Samps, "TPM")
      tempp <- abGeneTemp2[,cols2]
      val3 <- paste0("MeanTPMCond", i)
      abGeneTemp2[,val3] <- rowMeans(tempp)

      cols3 <- paste0(Samps, "Cnt")
      tempp <- cntGeneTemp2[,cols3]
      val4 <- paste0("MeanCntCond", i)
      cntGeneTemp2[,val4] <- rowMeans(tempp)

    }
  }


  #Merge gene/sample totals back in
  abGeneTemp3 <- merge(abGeneTemp2, abgenetotals, by = "gene_id")
  cntGeneTemp3 <- merge(cntGeneTemp2, cntgenetotals, by = "gene_id")


  #Create relative abundances (proportions) for each sample/gene combo
    #For use in the compositional regression analysis and possibly in later filtering
  if(is.null(samps)==FALSE){
    for (i in 1:length(samps)){
      abGeneTemp3[paste0(samps[i], "RTA")] <-  abGeneTemp3[paste0(samps[i], "TPM")]/abGeneTemp3[paste0(samps[i], "TGE")]
      cntGeneTemp3[paste0(samps[i], "RTA")] <-  cntGeneTemp3[paste0(samps[i], "Cnt")]/cntGeneTemp3[paste0(samps[i], "TGE")]
    }
  }else{
    for (i in 1:nsamp) {
      abGeneTemp3[paste0("Sample", i, "RTA")] <-  abGeneTemp3[paste0("Sample", i, "TPM")]/abGeneTemp3[paste0("Sample", i, "TGE")]
      cntGeneTemp3[paste0("Sample", i, "RTA")] <-  cntGeneTemp3[paste0("Sample", i, "Cnt")]/cntGeneTemp3[paste0("Sample", i, "TGE")]
    }
  }


  #extract which columns are the RTA columns
  if(is.null(samps)==FALSE){
    abRTAcols <- which(colnames(abGeneTemp3) %in% paste0(samps, "RTA"))
    cntRTAcols <- which(colnames(cntGeneTemp3) %in% paste0(samps, "RTA"))
  }else{
    abRTAcols <- which(colnames(abGeneTemp3) %in% paste0("Sample", 1:nsamp, "RTA"))
    cntRTAcols <- which(colnames(cntGeneTemp3) %in% paste0("Sample", 1:nsamp, "RTA"))
  }


  abGeneTemp3$meanRTA <- rowMeans(abGeneTemp3[,abRTAcols], na.rm = TRUE)
  cntGeneTemp3$meanRTA <- rowMeans(cntGeneTemp3[,cntRTAcols], na.rm = TRUE)


  #Merge in length information into count data file for possible use later
    #Keep in mind that the length is the effective length which varies per sample
  len <- lengths
  if(!("tx_id" %in% colnames(len))){
    if(is.null(samps)==FALSE){
      colnames(len) <- paste0(samps, "Len")
    }else{
      colnames(len) <- paste0("Sample", 1:nsamp, "Len")
    }

    len$tx_id <- rownames(len)
  }else{
    pos <- which(colnames(len)=="tx_id")
    if(is.null(samps)==FALSE){
      colnames(len)[-pos] <- paste0(samps, "Len")
    }else{
      colnames(len)[-pos] <- paste0("Sample", 1:nsamp, "Len")
    }

  }

  #len$tx_id <- sapply(strsplit(as.character(rownames(len)), "\\."), "[[", 1)
  cntGeneTemp4 <- merge(cntGeneTemp3, len, by = "tx_id")


  #Order/sort dataframe by gene then by trans within a gene
  abGeneTempF <- abGeneTemp3[order(abGeneTemp3$gene_id, abGeneTemp3$tx_id),]
  cntGeneTempF <- cntGeneTemp4[order(cntGeneTemp4$gene_id, cntGeneTemp4$tx_id),]

  if(infReps=="none"){
    rownames(abGeneTempF) <- abGeneTempF$tx_id
    rownames(cntGeneTempF) <- cntGeneTempF$tx_id
  }

  return(list(abGeneTempF = abGeneTempF, cntGeneTempF = cntGeneTempF))
}


#Generate the datasets with all information that will be needed to do downstream analyses
#Need to run this on abundance data first then use those results for count data
#Generally used with an apply loop and not called directly by user
#Agruments:
#x: genename
#dat: the observed counts/TPM data
#n: number of samples
#abundance; (T/F) is dat adundance (TPM) data or not
#abData: results from dat being abundance (abundance=T) (only non-NULL when running on count data)
#abCompDatasets: list of gene-wise dataframes output when dat is abundance (abundance=T)
#useExistingOtherGroups: (T/F) Should the "Other" groups from the observed abundance be used or new ones computed
#useOtherGroups: Should other groups be created (default Yes)
#infReps: One of "none", "Gibbs", and "Boot", corresponding to none, Gibbs samples, and Boot Samples respectively
generateData <- function(x, dat, nsamp, abundance, abData, abCompDatasets = NULL, useExistingOtherGroups, useOtherGroups = TRUE, 
                         useExistingMajorTrans = TRUE, infReps = "none", ninfreps = NA, samps = NULL, CompMI = FALSE){
  library(Matrix)
  #for (i in 1:length(fullgenenames)) {
  #x <- genestouse[1]
  temp <- subset(dat, dat$gene_id==x)

  #print(paste0("Currently Running generateData for gene", x))

  #If no rows in temp, no data exists for that gene
    #This can occur especially when using the Gibbs samples if no gibbs samples
    #exist for an entire gene that has "regular" observations for some reason
  if(nrow(temp)==0){
    return(NULL)
  }
  #temp <- subset(dat, dat$gene_id=="ENSG00000002919")

  if(infReps=="Gibbs" | infReps=="GibbsThin16" | infReps=="GibbsThin100"){
    rownames(temp) <- paste0(temp$tx_id, "Gibbs", temp$infRepNum)
    tnames <- unique(temp$tx_id)[unique(temp$tx_id) %in% attributes(abCompDatasets[[x]])$FullTrans]
  }else if(infReps=="Boot"){
    rownames(temp) <- paste0(temp$tx_id, "Boot", temp$infRepNum)
    tnames <- unique(temp$tx_id)[unique(temp$tx_id) %in% attributes(abCompDatasets[[x]])$FullTrans]
  }else if(infReps=="none"){
    tnames <- temp$tx_id
  }

  if(abundance==TRUE){
    if(is.null(samps)==FALSE){
      subcols <- which(colnames(temp) %in% paste0(samps, "TPM"))
    }else{
      subcols <- which(colnames(temp) %in% paste0("Sample", 1:nsamp, "TPM"))
    }

  }else{
    if(is.null(samps)==FALSE){
      subcols <- which(colnames(temp) %in% paste0(samps, "Cnt"))
    }else{
      subcols <- which(colnames(temp) %in% paste0("Sample", 1:nsamp, "Cnt"))
    }
  }

  GDinf <- function(x, data, subcols, infReps){
    data2 <- subset(data, data$infRepNum==x)
    data3 <- data2[,subcols, with = FALSE]

    data4 <- data.frame(data3)
    data5 <- t(data4)
    if(infReps=="Gibbs"| infReps=="GibbsThin16" | infReps=="GibbsThin100"){
      rownames(data5) <- paste0(rownames(data5), "Gibbs", x)
    }else if(infReps=="Boot"){
      rownames(data5) <- paste0(rownames(data5), "Boot", x)
    }
    colnames(data5) <- data2$tx_id
    data6 <- data.frame(data5)
    data6$id <- rownames(data6)
    return(data6)
  }

  if(infReps=="none"){

    sub <- temp[,subcols]
    rownames(sub) <- tnames
    d1 <- t(sub)

  }else{
    sub <- temp[,subcols, with = FALSE]
    rownames(sub) <- rownames(temp)
    d1t <- data.frame(rbindlist(apply(as.matrix(1:ninfreps), 1, GDinf, data = temp,
                                      subcols = subcols, infReps = infReps), use.names = TRUE))
    rownames(d1t) <- d1t$id
    d1t$id <- NULL

    #Use the same trans for Gibbs/Boot replicates datasets as is used in the observed datasets
      #Otherwise, the results may be hard to compare, as the set of transcripts could be different 
    d1 <- as.matrix(d1t[,colnames(d1t) %in% attributes(abCompDatasets[[x]])$FullTrans], nrow = nsamp)
    colnames(d1) <- colnames(d1t)[colnames(d1t) %in% attributes(abCompDatasets[[x]])$FullTrans]
    rownames(d1) <- rownames(d1t)
  }



  #If d1 is NULL or 0, gene has no valid data and its dataframe is set to NULL
  if(is.null(ncol(d1))){
    d3 <- NULL
    return(d3)
  }

  if(ncol(d1)==0){
    d3 <- NULL
    return(d3)
  }

  #Ignore transcripts with total count across all samples of 0
    #Just drop these transcripts now so they aren't included in with the Other category if it is used
    #These couldn't ever have a "switching" event so we don't really need to include them
  #Use the existing cols if you want to use the existing other groups or if this is constructing
    #the abDatasets for the Gibbs replicates because the columns/number of columns has to match
    #for the analysis to work right
  if(useExistingOtherGroups==TRUE){
    cols <- which(colnames(d1) %in% attributes(abCompDatasets[[x]])$FullTrans)
    #cols <- 1:ncol(d1)
    #names(cols) <- colnames(d1)
  }else{
    cols <- colSums(d1)>0
  }


  #Need to create the return object in this weird way because of the default r behavior or making a
  # matrix with only one column into a numeric vector
  # Need to include drop = FALSE option so dimensions are not dropped in the case of 1 column
    #As dropping the dimensions would break the R code downstream
  d2 <- data.frame(as.matrix(d1[,cols], nrow = nsamp))
  colnames(d2) <- tnames[cols]

  #Choose which transcripts contribute to the "Other" category based on TPM to ensure the Other category includes the same
  # transcripts for counts and abundance data - otherwise the makeup of this Other category would differ between
  # counts and TPM category
  # So need to run the results on the adundance data first and load in those results into the count results
  if(useOtherGroups==TRUE & useExistingOtherGroups==FALSE){
    #Combine all transcripts with RTA across all samples less than 5% into an "Other" category if there is more than 1
      #If there is exactly 1 trans with RTA < 5%, just drop it for now
    fullcolnames <- colnames(d2)

    #Col #s that have RTA < 5% after summing across all samples- not the same as RTA quantity, which is calculted per sample
    rarecols <- which(colSums(d2)/sum(d2) < 0.05)
    rarecolnames <- colnames(d2)[rarecols]

    #Need to generate like this and not just -rarecols because there could be no rarecols
    # and code would break in that case
    notrarecolnames <- colnames(d2)[!(colnames(d2) %in% rarecolnames)]
  }else if(useOtherGroups==TRUE & useExistingOtherGroups==TRUE){
    fullcolnames <- attr(abCompDatasets[[x]], "FullTrans")
    rarecolnames <- attr(abCompDatasets[[x]], "OtherTrans")
    #rarecols <- which(colnames(d2) %in% attr(abCompDatasets[[x]], "OtherTrans"))
    rarecols <- attr(abCompDatasets[[x]], "OtherTrans")
    notrarecolnames <- attr(abCompDatasets[[x]], "NotOtherTrans")
  }else if(useOtherGroups==FALSE & is.null(abCompDatasets)){
    fullcolnames <- colnames(d2)
    rarecolnames <- NULL
    rarecols <- NULL
    notrarecolnames <- colnames(d2)
  }else if(useOtherGroups==FALSE & !is.null(abCompDatasets)){
    fullcolnames <- attr(abCompDatasets[[x]], "FullTrans")
    rarecolnames <- NULL
    rarecols <- NULL
    notrarecolnames <- NULL
  }

  #If length of tnames is not the same as length(union(rarecolnames, notrarecolnames)),
    #there is a transcript that is not a member of the other groups
    #that is also not in the gibbs replicates.  This will cause a problem when trying
    #to calculate the covariances on the ilr scale because the number of columns
    #won't match
    #So, in this case set all values of this transcript to 0 to be able to get the dataset
    #To generate and have the proper dimensions to match the regular counts/TPMs
  if((infReps!="none" & length(tnames)!=length(fullcolnames)) | CompMI==TRUE){
    misstrans <- fullcolnames[!(fullcolnames %in% tnames)]
    if(length(misstrans)!=0){
      for(l in 1:length(misstrans)){
        if(is.na(misstrans[l])){next}
        print(paste0("A transcript that is missing in infRep data but not in regular had to be re-inserted with all 0s"))
        d2[,misstrans[l]] <- 0
      }
    }
  }

  # If d2 has no columns now, move on to the next gene
  if(ncol(d2)==0){
    return(NULL)
  }


  #If there is only 1 rare column, for now drop that trans instead of combining it with anything else
    #d3 is the count or abundance data frame that will be output for each gene

  #If no trans with RTA < 5%, just use d2 dataset from above
  if(length(rarecols)==0){
    d3 <- d2
  }

  #If exactly 1 trans with a trans with RTA < 5%, just drop that trans for now to avoid needing to combine with a higher one
    #Could maybe reevaluate this or perhaps add an option for how this could be handled?
  if(length(rarecols)==1){
    d3 <- d2[notrarecolnames]
  }

  #If multiple trans with RTA < 5%, combine them into the "Other" category
  if(length(rarecols) > 1){
    d2$Other <- as.numeric(rowSums(d2[,rarecols]))
    d3 <- d2[c(notrarecolnames, "Other")]
  }

  #As of now, don't do anything in this file if all trans are "rare" (ie all have RTA <5%)
    #just filter that out before doing any analysis on it- this should be very rare
    #This is a problem since isn't really possible to pick a transcript to choose to not be included in the other category
    #And having only an Other category means we can't do isoform switching analysis (and it looses interpretation also)

  #If there are no cols left, all transcripts are "rare" and the gene can't be used accurately as of now
    #So, just return it as null and go on to the next gene
  if(ncol(d3)==0){
    d3 <- NULL
    return(d3)
  }

  #Use the MajorTrans from the observed TPM values to be the major trans both for the counts
    #and for any data based on the Gibbs replicates
  #If there are multiple "major" transcripts (ie the meanRTA for two transcripts across all samples
    #are exactly tied for the max value) just assign the first transcript of those to be the major transcript
  if(useExistingMajorTrans==FALSE){
    #Now, return column number and name of the major transcript (transcript with highest
      # average TPM across) but don't allow the "Other" category to be the major trans
    if(length(notrarecolnames)==1){
      majtrans <- notrarecolnames
    }else{
      colm <- colMeans(d3[,notrarecolnames], na.rm = TRUE)
      majtrans <- names(which(colm==max(colm)))
    }

    if(length(majtrans) >1){
      majtrans <- majtrans[1]
    }

    attr(d3, "MajorTrans") <-   majtrans
  }else{
    #Use the major trans from TPM for the count data so they are the same
    attr(d3, "MajorTrans") <- attr(abCompDatasets[[x]], "MajorTrans")
  }

  if(length(attr(d3, "MajorTrans")) >1){stop("Something is wrong with the major transcript calculation, there are multiple of them when there should only be 1")}


  #Remember that a transcripts membership in the Other group is determined from observed TPMs always, not the counts
    #or from any specific Gibbs replicate that may be currently being used
  attr(d3, "OtherTrans") <- rarecolnames
  attr(d3, "NotOtherTrans") <- notrarecolnames
  attr(d3, "FullTrans") <- fullcolnames


  #calculating rank for the Gibbs/Boot datasets causes computational issues, and is never needed based on the way the analysis is done, so don't bother calculating it
  if(infReps=="none"){
    attr(d3, "FullRank") <- (rankMatrix(d3) == ncol(d3))
  }

  ##################################################################################################################
  #Now ,create length data files for use with count data in case they would ever be needed
  #Don't include these with abundance data (since theyve already been normalized)
  #And, don't include if using Gibbs/Boot Samps since that will just serve to take up space and these
  #can be loaded from the regular cntDatasets file
  ##################################################################################################################
  if(abundance==FALSE & infReps=="none"){
    #First, get the lengths

    subcols2 <- which(colnames(temp) %in% paste0("Sample", 1:nsamp, "Len"))
    sub2 <- temp[,subcols2]
    rownames(sub2) <- tnames
    ln <- data.frame(t(sub2))

    #Again, drop transcripts that are always 0- these are the ones in cols from above
    #Need to create the return object in this weird way because of the default r behavior or making a
    # matrix with only one column into a numeric vector
    # This would break code downstream that relies on dim(), etc so need to create data this way
    ln2 <- data.frame(as.matrix(ln[,fullcolnames], nrow = nsamp))

    #As before, if there are no trans with RTA < 5% just use ln result
    if(length(rarecols)==0){
      ln3 <- ln2
      colnames(ln3) <- notrarecolnames
    }

    #With exactly one trans with <5% RTA, drop for now- could change later
    if(length(rarecols)==1){
      ln3 <- data.frame(as.matrix(ln2[,notrarecolnames], nrow = nsamp))
      colnames(ln3) <- notrarecolnames

    }

    #For >=2 RTA, create other category, with the length for that category just
    #being the sum of the offsets for the trans that make up the other category
    if(length(rarecols) > 1){
      ln2$Other <- rowSums(ln2[,rarecols])
      ln3 <- ln2[,c(notrarecolnames, "Other")]
    }

    rownames(ln3) <- paste0("Sample", 1:nsamp, "Len")
    if(useOtherGroups==FALSE){
      colnames(ln3) <- colnames(d3)
    }
    attr(ln3, "OtherTrans") <- rarecolnames
    attr(ln3, "NotOtherTrans") <- notrarecolnames
    attr(ln3, "FullTrans") <- fullcolnames

    ###############################################################################################
    #This code below was creating an offset, which as of now is never used so is being skipped
    #Offset as of now is the log difference between
    #the count value and the the TPM value (ie offset = log(TPM/Count))- check on this
    #To get the offset for the "other" category, for now sum TPM for the categories used in other and
    #sum counts for the other categories and do Offset =log(sum(TPM)/sum(Count)), where
    #the sums are over all transcripts that make up the other category
    ###############################################################################################

    # tempoff <- subset(abData, abData$gene_id==x)
    # #tempoff <- subset(abData, abData$gene_id=="ENSG00000002919")
    #
    # tnames <- tempoff$tx_id
    #
    # subcols3 <- which(colnames(tempoff) %in% paste0("Sample", 1:nsamp, "TPM"))
    # suboff <- tempoff[,subcols3]
    # rownames(suboff) <- tnames
    #
    # off1 <- t(suboff)
    #
    # #Drop columns (transcripts) that always have 0 expression, as is done above with creating the count object
    # off2 <- as.matrix(off1[,fullcolnames], nrow = nsamp)
    # colnames(off2) <- fullcolnames
    # rownames(off2) <- paste0("Sample", 1:nsamp)
    #
    # denom <- as.matrix(d3[,notrarecolnames], nrow = nsamp)
    # colnames(denom) <- notrarecolnames
    # rownames(denom) <- paste0("Sample", 1:nsamp)
    #
    # #Offset for "regular" (not rare transcripts) is TPM value divided by count value
    # OffsetNotOther <- data.frame(as.matrix(off2[,notrarecolnames]/denom, nrow = nsamp))
    # colnames(OffsetNotOther) <- notrarecolnames
    # rownames(OffsetNotOther) <- paste0("Sample", 1:nsamp)
    #
    #
    # #Offset for "Other" transcripts is the sum of TPM for the other trans divided
    # #by the sum of counts for other trans
    # #Note that if there is exactly 1 rare transcript, it is just dropped above instead of trying to combine
    # #it with another one.  So, only create this "Other" group offset if there are 2 or more raretrans
    # if(length(rarecolnames) > 1){
    #   AbOther <- data.frame(as.matrix(rowSums(as.matrix(off2[,rarecols], nrow = nsamp)), nrow = nsamp, ncol = 1))
    #   colnames(AbOther) <- "Other"
    #   rownames(AbOther) <- paste0("Sample", 1:nsamp, "Off")
    #   OffsetOther <- AbOther/(as.matrix(d3$Other, nrow = nsamp))
    #   #OffsetFull <- log(cbind(OffsetNotOther, OffsetOther))
    #   OffsetFull <- cbind(OffsetNotOther, OffsetOther)
    # }else{
    #   #OffsetFull <- log(OffsetNotOther)
    #   OffsetFull <- OffsetNotOther
    # }
    #
    # if(ncol(d3) != ncol(OffsetFull)){
    #   stop(print(paste("ncol of Counts not equal to ncol of Offsets, somethings wrong with gene", x)))
    # }
    #
    # #Return a list with 1st element being the counts, second element being the lengths, third being Offset
    # ret <- list(d3, ln3, OffsetFull)
    # names(ret) <- c("Counts", "Lengths", "Offsets")


    ret <- list(d3, ln3)
    names(ret) <- c("Counts", "Lengths")


    #dataloop[[i]] <- ret
    return(ret)
  }else{
    return(d3)
  }

} #end generateData


#This function reads in the major transcript information from abDatasets and adds it into the temporary
  #abGene and cntGene files
addMajorTrans <- function(genestouse, abGeneTempF, cntGeneTempF, abDatasets, CompMI = FALSE){
  MajorTrans <- data.frame(genestouse, NA)
  colnames(MajorTrans) <- c("gene_id", "MajorTransName")
  MajorTrans$MajorTransName <- lapply(genestouse, function(x) {attr(abDatasets[[x]], "MajorTrans")})

  if(CompMI == TRUE){
    abGene <- merge(abGeneTempF, MajorTrans, by = "gene_id", all = FALSE)
    abGene$MajorTrans <- as.numeric(abGene$MajorTransName==abGene$tx_id)
    rownames(abGene) <- abGene$tx_id

    cntGene <- merge(cntGeneTempF, MajorTrans, by = "gene_id", all = FALSE)
    cntGene$MajorTrans <- as.numeric(cntGene$MajorTransName==cntGene$tx_id)
    rownames(cntGene) <- cntGene$tx_id
  }else{
    abGene <- merge(abGeneTempF, MajorTrans, by = "gene_id", all = TRUE)
    abGene$MajorTrans <- as.numeric(abGene$MajorTransName==abGene$tx_id)
    rownames(abGene) <- abGene$tx_id

    cntGene <- merge(cntGeneTempF, MajorTrans, by = "gene_id", all = TRUE)
    cntGene$MajorTrans <- as.numeric(cntGene$MajorTransName==cntGene$tx_id)
    rownames(cntGene) <- cntGene$tx_id
  }


  abGene <-  abGene[order(abGene$gene_id, abGene$tx_id),]
  cntGene <- cntGene[order(cntGene$gene_id, cntGene$tx_id),]
  return(list(abGene = abGene, cntGene = cntGene))
}



#This function filters like DRIMSeq and saves the results
  #If abFromInfRepFunc is specified, the set of transcripts that passes filtering is fixed to be the same as those that pass it 
  #from the regular abundnace data to keep the results as consistent as possible
  #In this case the abGene and cntGene files need to be pre filtered using the list of filtered transcript names before plugging into this function
  #Also, specify abFromInfRepFunc and GibbsSamps as non-null in this case too
DRIMSeqFilter <- function(abGene, cntGene, failedgibbs = NULL, abFromInfRepFunc = NULL, GibbsSamps = NULL, abGeneRegularAbValues = NULL){

  if(is.null(abFromInfRepFunc)){
    temp1 <- PreDRIMSeq(cntGene = cntGene, failedgibbssamps = failedgibbs, key = key)
    cnts <- temp1$cnts
    samp <- temp1$samp
    
    if(is.null(sampstouse)){
      samp2 <- samp
    }else{
      samp2 <- subset(samp, samp$sample_id %in% sampstouse)
      samp2$group <- relevel(factor(samp2$group), ref = 1)
    }
    
    DRIMData <- dmDSdata(counts = cnts, samples = samp2)
    
    DRIMData2 <- dmFilter(DRIMData,
                          min_samps_feature_expr=min_samps_feature_expr, min_feature_expr=min_feature_expr,
                          min_samps_feature_prop=min_samps_feature_prop, min_feature_prop=min_feature_prop,
                          min_samps_gene_expr=min_samps_gene_expr, min_gene_expr=min_gene_expr)
    
    
    DRIMSeqFilteredData <- counts(DRIMData2)
    transcriptstouse <- DRIMSeqFilteredData$feature_id
    fullgenenames_filtered <- unique(DRIMSeqFilteredData$gene_id)
  }else{
    print(paste0("Results are being run for abundances calculated using the inferential replicates"))
    transcriptstouse <- abGeneRegularAbValues$tx_id
    fullgenenames_filtered <- unique(abGeneRegularAbValues$gene_id)
  }
  

  counts <- cntGene[cntGene$tx_id %in% transcriptstouse, c("tx_id", paste0(key$Identifier, "Cnt"))]
  rownames(counts) <- counts$tx_id
  counts$tx_id <- NULL
  colnames(counts) <- key$Identifier
  abundance <- abGene[abGene$tx_id %in% transcriptstouse, c("tx_id", paste0(key$Identifier, "TPM"))]
  rownames(abundance) <- abundance$tx_id
  abundance$tx_id <- NULL
  colnames(abundance) <- key$Identifier
  lengths <- cntGene[cntGene$tx_id %in% transcriptstouse, c("tx_id", paste0(key$Identifier, "Len"))]
  rownames(lengths) <- lengths$tx_id
  lengths$tx_id <- NULL
  colnames(lengths) <- key$Identifier

  ST1 <- proc.time()
  initialData <- prepareData(abundance = abundance, counts = counts, lengths = lengths, tx2gene = tx2gene,
                             nsamp = nsamp, key = key)

  abGenecntGeneCompTimeP1 <- proc.time() - ST1
  abGeneTempF <- initialData$abGeneTempF
  cntGeneTempF <- initialData$cntGeneTempF

  #library(parallel)
  #clust <- makeCluster(1)
  #FilteredDat <- sumToGeneHelper(abundance = abundance, counts = counts, lengths = lengths, tx2gene = tx2gene, Group = Group, clust = NULL, nsamp = length(Group), key = key, useOtherGroups = FALSE, useExistingOtherGroups = FALSE, useExistingMajorTrans = FALSE)

  #abDatasets are only input to extract existing MajorTrans information from
  FilteredDat <- sumToGeneHelper(abundance = abundance, counts = counts, lengths = lengths, tx2gene = tx2gene, Group = Group, clust = NULL, nsamp = length(Group),
                                 key = key, useOtherGroups = FALSE, useExistingOtherGroups = FALSE, useExistingMajorTrans = FALSE, abCompDatasets = NULL)


  abDatasetsFiltered <- FilteredDat$abDatasets
  cntDatasetsFiltered <- FilteredDat$cntDatasets

  abDatasetsFilteredCompTime <- FilteredDat$abDatasetsCompTime
  cntDatasetsFilteredCompTime <- FilteredDat$cntDatasetsCompTime

  abGeneFiltered <- FilteredDat$abGene
  cntGeneFiltered <- FilteredDat$cntGene

  NTransFilabGene <- data.frame(table(abGeneFiltered$gene_id))
  colnames(NTransFilabGene) <- c("gene_id", "NTransFiltered")

  NTransFilcntGene <- data.frame(table(cntGeneFiltered$gene_id))
  colnames(NTransFilcntGene) <- c("gene_id", "NTransFiltered")

  abGeneFiltered <- merge(abGeneFiltered, NTransFilabGene, by = "gene_id")
  cntGeneFiltered <- merge(cntGeneFiltered, NTransFilcntGene, by = "gene_id")

  rownames(abGeneFiltered) <- abGeneFiltered$tx_id
  rownames(cntGeneFiltered) <- cntGeneFiltered$tx_id
  
  
  
  
  if(is.null(abFromInfRepFunc)){
    fil_mod <- ""
  }else if(abFromInfRepFunc=="median" & GibbsSamps==TRUE){
    fil_mod <- "AbRowMedianInfRepGibbs"
  }else if(abFromInfRepFunc=="mean" & GibbsSamps==TRUE){
    fil_mod <- "AbRowMeanInfRepGibbs"
  }else if(abFromInfRepFunc=="median" & GibbsSamps==FALSE){
    fil_mod <- "AbRowMedianInfRepBoot"
  }else if(abFromInfRepFunc=="mean" & GibbsSamps==FALSE){
    fil_mod <- "AbRowMeanInfRepBoot"
  }else{
    stop("Check specification of abFromInfRepFunc function")
  }
  
  

  if(countsFromAbundance=="scaledTPM"){
    save(abGeneFiltered, nsamp, key, file = paste0("abGeneFiltered", fil_mod, ".RData"))
    save(cntGeneFiltered, nsamp, countsFromAbundance, key, file = paste0("cntGenecntsScaledTPMFiltered", fil_mod, ".RData"))

    save(abDatasetsFiltered, nsamp, key, fullgenenames_filtered, Group, abDatasetsFilteredCompTime, file = paste0("abDatasetsNoOtherGroupsFiltered", fil_mod, ".RData"))
    save(cntDatasetsFiltered, nsamp, key, fullgenenames_filtered, countsFromAbundance, Group, cntDatasetsFilteredCompTime, file = paste0("cntDatasetscntsScaledTPMNoOtherGroupsFiltered", fil_mod, ".RData"))
  }else if(countsFromAbundance=="no"){
    save(abGeneFiltered, nsamp, key, abGeneFilteredCompTime, file = paste0("abGeneFiltered", fil_mod, ".RData"))
    save(cntGeneFiltered, nsamp, countsFromAbundance, key, cntGeneFilteredCompTime, file = paste0("cntGeneFiltered", fil_mod, ".RData"))

    save(abDatasetsFiltered, nsamp, key, fullgenenames_filtered, Group, abDatasetsFilteredCompTime, file = paste0("abDatasetsNoOtherGroupsFiltered", fil_mod, ".RData"))
    save(cntDatasetsFiltered, nsamp, key, fullgenenames_filtered, countsFromAbundance, Group, cntDatasetsFilteredCompTime, file = paste0("cntDatasetsNoOtherGroupsFiltered", fil_mod, ".RData"))
  }
}




#Now, filter the datasets using the values where the abundances are calculated as the mean/median of the bootstrap/gibbs samples
returnFilMod <- function(abFromInfRepFunc, GibbsSamps){
  if(is.null(abFromInfRepFunc)){
    fil_mod <- ""
  }else if(abFromInfRepFunc=="median" & GibbsSamps==TRUE){
    fil_mod <- "AbRowMedianInfRepGibbs"
  }else if(abFromInfRepFunc=="mean" & GibbsSamps==TRUE){
    fil_mod <- "AbRowMeanInfRepGibbs"
  }else if(abFromInfRepFunc=="median" & GibbsSamps==FALSE){
    fil_mod <- "AbRowMedianInfRepBoot"
  }else if(abFromInfRepFunc=="mean" & GibbsSamps==FALSE){
    fil_mod <- "AbRowMeanInfRepBoot"
  }else{
    stop("Check specification of abFromInfRepFunc function")
  }
  return(fil_mod)
}

#Run the DRIMSeq filtering approach on the datasets with 
runDRIMSeqFilterOnAbFromInfReps <- function(abFromInfRepFuncT, GibbsSampsT, countsFromAbundance = "scaledTPM", abGeneRegularAbValues){
  abFromInfRepFunc <- abFromInfRepFuncT
  GibbsSamps <- GibbsSampsT
  fil_mod <- returnFilMod(abFromInfRepFunc = abFromInfRepFunc, GibbsSamps = GibbsSamps)
  
  abGeneT <- loadRData(paste0("abGene", fil_mod, ".RData"), objNameToGet = "abGene")
  
  if(countsFromAbundance=="scaledTPM"){
    cntGeneT <- loadRData(paste0("cntGenecntsScaledTPM", fil_mod, ".RData"), objNameToGet = "cntGene")
  }else if(countsFromAbundance=="no"){
    cntGeneT <- loadRData(paste0("cntGene", fil_mod, ".RData"), objNameToGet = "cntGene")
  }
  
  abGeneFiltered <- loadRData(paste0(dir1, "abGeneFiltered.RData"), objNameToGet = "abGeneFiltered")
  
  abGene <- subset(abGeneT, abGeneT$tx_id %in% abGeneFiltered$tx_id)
  cntGene <- subset(cntGeneT, cntGeneT$tx_id %in% abGeneFiltered$tx_id)
  
  DRIMSeqFilter(abGene = abGene, cntGene = cntGene, abFromInfRepFunc = abFromInfRepFunc, GibbsSamps = GibbsSamps, abGeneRegularAbValues = abGeneRegularAbValues)
}



#This function transforms from a completed cntDatasets object to an updated cnts file
  #This is useful for a getting a cnts file with genes/transcripts groups (including other groups) that exactly match
  #cntDatasets.  This is useful in giving the same values all around to then be run in DRIMSeq to be able to more directly
  #compare Compositional Methods to DRIMSeq without the confounder of different genes/trans being used or the fact that DRIMSeq
  #doesn't create "Other" Groups
  #Note that the count values that will be used eventually are contained in the cntGene file except for the "Other" group which
  #come fron cntDatasets.  This is ok even for the power analyses where the counts in major trans are changed because the changing
  #is done within cntGene such that the new counts would be successfully carried through this way and because the OtherGroups
  #are never changed in any of the power analyses such that it is ok to still pull OtherGroups from the cntDatasets files

convertcntGeneToOtherGroupsForDRIMSeq <- function(cntGene, cntDatasets, samps = NULL, infReps = "none"){
 genestouse <- names(cntDatasets)
 print(samps)
 if(is.null(samps)){
   sampscolnames <- paste0("Sample", 1:nsamp, "Cnt")
 }else{
   sampscolnames <- paste0(samps, "Cnt")
   nsamp <- length(samps)
 }


 colstokeep <- c("gene_id", "tx_id", sampscolnames , "NTrans", "SumTGE", "MajorTransName", "MajorTrans")


 newcntsfun <- function(x, colstokeep, cntGene, cntDatasets, infReps){
   curr_gene <- x
   if(is.null(cntDatasets[[curr_gene]])){
     return(NULL)
   }
   if(infReps=="none"){
     cntDatasetsSub <- cntDatasets[[curr_gene]]$Counts
   }else{
     cntDatasetsSub <- cntDatasets[[curr_gene]]
   }
   OtherTrans <- attr(cntDatasetsSub, "OtherTrans")
   NotOtherTrans <- attr(cntDatasetsSub, "NotOtherTrans")
   FullTrans <- attr(cntDatasetsSub, "FullTrans")
   cntgenesub <- cntGene[FullTrans, colstokeep]

   if(("Other" %in% colnames(cntDatasetsSub))==FALSE){
     cntgenesub2 <- cntgenesub[colnames(cntDatasetsSub),]
     return(cntgenesub2)
   }
   #cntgenesub <- subset(cntGene[,colstokeep], cntGene$gene_id==curr_gene)


   if(length(OtherTrans)==0){
     return(cntgenesub)
   }

   #Update the final result to combine other category if it is present
   cntgenesub2 <- cntgenesub[c(NotOtherTrans, OtherTrans[1]),]
   cntgenesub2[OtherTrans[1],"tx_id"] <- paste0("Other", curr_gene)
   rownames(cntgenesub2) <- cntgenesub2$tx_id

   for(i in 1:nsamp){
     curr_samp <- samps[i]
     if(infReps=="none"){
       curr_nam <- paste0(curr_samp, "Cnt")
       othernam <- paste0("Other", curr_gene)
     }else if(infReps=="Boot"){
       curr_nam <- paste0(curr_samp, "Cnt", "Boot", 1:ninfreps)
       othernam <- paste0("Other", curr_gene, "Boot", 1:ninfreps)
     }else if(infReps=="Gibbs"){
       curr_nam <- paste0(curr_samp, "Cnt", "Gibbs", 1:ninfreps)
       othernam <- paste0("Other", curr_gene, "Gibbs", 1:ninfreps)
     }

     if(infReps=="none"){
       #Update the count for the "Other Group" to be the value for the "Other" group from the correspondng cntDataset
       #This is fine even if the data has been changed for a power analysis because the OtherGroup values are never changed
       cntgenesub2[othernam, curr_nam] <- cntDatasetsSub[curr_nam, "Other"]
     }else{
       for(j in 1:length(curr_nam)){
         cntgenesub2[othernam[j], curr_nam[j]] <- cntDatasetsSub[curr_nam[j], "Other"]
       }
     }

   }
   rownames(cntgenesub2) <- 1:nrow(cntgenesub2)
   return(cntgenesub2)
 }

 CntsTemp <- lapply(genestouse, newcntsfun, colstokeep, cntGene, cntDatasets, infReps = infReps)
 print(colnames(CntsTemp[[1]]))
 print(colnames(CntsTemp[[2]]))
 #restemp <- laply(genestouse[6000:20000], newcntsfun, colstokeep, cntGene, cntDatasets, .inform = T)
 Cnts <- data.frame(rbindlist(CntsTemp))
 return(Cnts)
 }


#Inputs
#cnts dataframe of counts with columns with names Sample1Cnt, Sample2Cnt, ... and a column for tx_id
#Additionally, if this does nt have length columns (Sample1Len, Sample2Len, etc) specify len argument
#nsamp: number of biological samples
#len: An optional argument with length information.  dataframe with (Sample1Len, Sample2Len, etc) and a tx_id column at the end
#samps: An optional vector containing the sample names.  Need this if sample names are not just paste0("Sample", 1:nsamp) without any missing.
cntsToTPM <- function(cnts, nsamp, len = NULL, samps = NULL, returnScaledTPMCounts = FALSE){
  #print("head of samps within cntsToTPM is ")
  #print(head(samps))
  if(is.null(len)==TRUE){
    cnts2 <- cnts
  }else{
    #Confirm that the transcripts are in the same order in the new TPM and abGene
    if(sum(rownames(len)!=rownames(cnts)) != 0){
      stop("Transcript row names are not aligned, dataframes are not in the same order")
    }
    len2 <- data.frame(len)
    if(is.null(len2$tx_id)){
      len2$tx_id <- rownames(len2)
    }
    if(is.null(cnts$tx_id)){
      cnts$tx_id <- rownames(cnts)
    }
    cnts2 <- merge(cnts, len2, by = "tx_id")
    rownames(cnts2) <- cnts2$tx_id
    #lengthCols <- which(colnames(len) %in% paste0("Sample", 1:nsamp, "Len"))
    #t3 <- len[,lengthCols]
  }
  if(is.null(samps)==TRUE){
    countCols <- which(colnames(cnts2) %in% paste0("Sample", 1:nsamp, "Cnt"))
    lengthCols <- which(colnames(cnts2) %in% paste0("Sample", 1:nsamp, "Len"))
  }else{
    countCols <- which(colnames(cnts2) %in% paste0(samps, "Cnt"))
    lengthCols <- which(colnames(cnts2) %in% paste0(samps, "Len"))
  }

  #print("head of countCols in cntsToTPM is ")
  #print(head(countCols))

  #print("head of lengthCols in cntsToTPM is ")
  #print(head(lengthCols))

  t2 <- cnts2[,countCols, drop = FALSE]
  t3 <- cnts2[,lengthCols, drop = FALSE]

  #print("head of colnames of t2 (countscols) in cntsToTPM is ")
  #print(head(colnames(t2)))
  #print("head of colnames of t3 (lengthcols) in cntsToTPM is ")
  #print(head(colnames(t3)))

  z <- t2/t3 # Will do element wise division of counts and (effective) lengths as desired

  #print(paste0("length of samps within cntstoTPM is ", length(samps)))
  #print(paste0("nsamp within cntstoTPM is ", nsamp))
  if(returnScaledTPMCounts==TRUE){
    if(is.null(samps)==TRUE){
      ScaledTPMCounts <- as.data.frame(apply(as.matrix(1:nsamp), 1, function(x) {(z[,x]/sum(z[,x])) * sum(t2[,x])}))
    }else{
      ScaledTPMCounts <- as.data.frame(apply(as.matrix(1:length(samps)), 1, function(x) {(z[,x]/sum(z[,x])) * sum(t2[,x])}))
    }
    
    #print("head of colnames of TPM before modification within cntsToTPM is ")
    #print(head(colnames(TPM)))
    
    rownames(ScaledTPMCounts) <- rownames(t2)
    if(is.null(samps)==TRUE){
      colnames(ScaledTPMCounts) <- paste0("Sample", 1:nsamp, "Cnt")
    }else{
      colnames(ScaledTPMCounts) <- paste0(samps, "Cnt")
    }
    
    #print("head of colnames of TPM after modification within cntsToTPM is ")
    #print(head(colnames(TPM)))
    ScaledTPMCounts$tx_id <- rownames(ScaledTPMCounts)
    return(ScaledTPMCounts)
    
  }else if(returnScaledTPMCounts==FALSE){
    if(is.null(samps)==TRUE){
      TPM <- as.data.frame(1e6 * apply(as.matrix(1:nsamp), 1, function(x) {z[,x]/sum(z[,x])}))
    }else{
      TPM <- as.data.frame(1e6 * apply(as.matrix(1:length(samps)), 1, function(x) {z[,x]/sum(z[,x])}))
    }
    
    #print("head of colnames of TPM before modification within cntsToTPM is ")
    #print(head(colnames(TPM)))
    
    rownames(TPM) <- rownames(t2)
    if(is.null(samps)==TRUE){
      colnames(TPM) <- paste0("Sample", 1:nsamp, "TPM")
    }else{
      colnames(TPM) <- paste0(samps, "TPM")
    }
    
    #print("head of colnames of TPM after modification within cntsToTPM is ")
    #print(head(colnames(TPM)))
    TPM$tx_id <- rownames(TPM)
    return(TPM)
  }
  
}


cntsToScaledTPMCounts <- function(cnts, nsamp, len = NULL, samps = NULL){
  scaledTPMCounts <- cntsToTPM(cnts = cnts, nsamp = nsamp, len = len, samps = samps, returnScaledTPMCounts = TRUE)
  return(scaledTPMCounts)
}

#x is the change value here
changeCountsAndTPM <- function(x, cntGene, abGene, key, Group, seedtoset, tx2gene, half = TRUE, genestochange = NULL, samps = NULL){
  
  nsamp <- nrow(key)
  
  #Change the counts for level 1 of the condition, "A" for the SQCC Data

  #Setting the seed everytime will result in the same sample chosen each time, which for now is what
  # is desired- can change this later if needed
  if(!is.null(seedtoset)){
    set.seed(seedtoset)
  }

  CompCntGene <- subset(cntGene, (cntGene$NTrans!=1 & cntGene$SumTGE!=0))
  fullgenenames <- sort(unique(CompCntGene$gene_id))

  #Only change for half of the genes to be able to generate a ROC curve
  if((half==TRUE & is.null(genestochange))){
    genestochange <- sort(sample(fullgenenames, ceiling(length(fullgenenames)/2)))
  }

  if(half==FALSE & is.null(genestochange)){
    genestochange <- fullgenenames
  }


  if(x==1){
    #Add information on which genes were "changed" to the dataframe
    #With a value of 1, nothing is actually changed, this is only added to not break code downstream
    cntGenecntsScaledTPM <- cntsToScaledTPMCounts(cnts = cntGene, nsamp = nsamp, samps = samps)
    attr(abGene, "geneschanged") <- genestochange
    attr(cntGene, "geneschanged") <- genestochange
    attr(cntGenecntsScaledTPM, "geneschanged") <- genestochange
    return(list(abGene = abGene, cntGene = cntGene, geneschanged = genestochange, cntGenecntsScaledTPM = cntGenecntsScaledTPM))
  }

  SampsCond1 <- subset(key$Identifier, Group==levels(key$Condition)[1])
  cols <- paste0(SampsCond1, "Cnt")
  ColsCond1 <- which(colnames(cntGene) %in% cols)
  
  if(is.null(samps)==TRUE){
    samps <- key$Identifier
  }

  cntGeneUpdated <- cntGene

  rws <- (cntGeneUpdated$MajorTrans==1 & cntGeneUpdated$gene_id %in% genestochange)

  #x is the change value
  cntGeneUpdated[rws,ColsCond1] <- cntGeneUpdated[rws,ColsCond1] * x

  #Calculate new TPMs from the modified counts as well as ScaledTPM Counts.  These need to be based on counts for all transcripts
    #and genes (ie from cntGeneUpdated, not CompCntGene)
  print(paste0("Number of rows for cntGeneUpdated that is used to calculate new TPMs is ", nrow(cntGeneUpdated)))
  TPM <- cntsToTPM(cnts = cntGeneUpdated, nsamp = nsamp, samps = samps)
  ScaledTPMCounts <- cntsToScaledTPMCounts(cnts = cntGeneUpdated, nsamp = nsamp, samps = samps)

  #Confirm that the transcripts are in the same order in the new TPM and abGene
  if(sum(rownames(TPM)!=rownames(abGene)) != 0){
    stop("Transcript row names are not aligned, dataframes are not in the same order")
  }
  
  if(sum(rownames(TPM)!=rownames(ScaledTPMCounts)) != 0){
    stop("Transcript row names are not aligned, dataframes are not in the same order")
  }


  #Now that we have the counts are updated, we need to update and recalculate the
    #total gene expression (TGE), RTA, etc
  newCntCols <- c("tx_id", paste0(samps, "Cnt"))
  newAbCols <- c("tx_id", paste0(samps, "TPM"))
  newLengthCols <- c("tx_id", paste0(samps, "Len"))

  #print("Head of AbCol Names is ")
  #print(head(newAbCols))

  #print("Head of CntCol Names is ")
  #print(head(newCntCols))

  #print("Head of Length Names is ")
  #print(head(newLengthCols))

  ab <- TPM[,newAbCols]
  cnt <- cntGeneUpdated[,newCntCols]
  leng <- cntGeneUpdated[,newLengthCols]
  
  cntScaledTPM <- ScaledTPMCounts[,newCntCols]

  cnam <- c("tx_id", samps)

  colnames(ab) <- cnam
  colnames(cnt) <- cnam
  colnames(leng) <- cnam
  colnames(cntScaledTPM) <- cnam

  #print(head(cnam))

  NewRes <- prepareData(abundance = ab, counts = cnt, lengths = leng, tx2gene = tx2gene, nsamp = nsamp, key = key, samps = samps)
  abGeneNewTemp <- NewRes$abGeneTempF
  cntGeneNewTemp <- NewRes$cntGeneTempF
  
  NewRes2 <- prepareData(abundance = ab, counts = cntScaledTPM, lengths = leng, tx2gene = tx2gene, nsamp = nsamp, key = key, samps = samps)
  cntGenecntsScaledTPMNewTemp <- NewRes2$cntGeneTempF
  
  #Now, add in the major transcript info, where the major transcript was found based on the original data
  majtcols <- c("tx_id", "MajorTransName", "MajorTrans")
  MajTransInfo <- abGene[,majtcols]

  abGeneNew <- merge(abGeneNewTemp, MajTransInfo, by = "tx_id")
  cntGeneNew <- merge(cntGeneNewTemp, MajTransInfo, by = "tx_id")
  cntGenecntsScaledTPMNew <- merge(cntGenecntsScaledTPMNewTemp, MajTransInfo, by = "tx_id")

  abGeneNew2 <- as.data.frame(abGeneNew[order(abGeneNew$gene_id, abGeneNew$tx_id),])
  rownames(abGeneNew2) <- abGeneNew2$tx_id

  cntGeneNew2 <- as.data.frame(cntGeneNew[order(cntGeneNew$gene_id, cntGeneNew$tx_id),])
  rownames(cntGeneNew2) <- cntGeneNew2$tx_id
  
  cntGenecntsScaledTPMNew2 <- as.data.frame(cntGenecntsScaledTPMNew[order(cntGenecntsScaledTPMNew$gene_id, cntGenecntsScaledTPMNew$tx_id),])
  rownames(cntGenecntsScaledTPMNew2) <- cntGenecntsScaledTPMNew2$tx_id

  if(!is.null(genestochange)){
    abGeneNew2$genechanged <- as.numeric(abGeneNew2$gene_id %in% genestochange)
    cntGeneNew2$genechanged <- as.numeric(cntGeneNew2$gene_id %in% genestochange)
    cntGenecntsScaledTPMNew2$genechanged <- as.numeric(cntGenecntsScaledTPMNew2$gene_id %in% genestochange)
  }

  #Add information on which genes were changed to the dataframe
  attr(abGeneNew2, "geneschanged") <- genestochange
  attr(cntGeneNew2, "geneschanged") <- genestochange
  attr(cntGenecntsScaledTPMNew2, "geneschanged") <- genestochange
  #save(abGeneNew, file = "Temp.RData")
  return(list(abGene = abGeneNew2, cntGene = cntGeneNew2, geneschanged = genestochange, cntGenecntsScaledTPM = cntGenecntsScaledTPMNew2))
}

#for(i in 1:length(changes)){
#x is a change quantity here
#abDatasetsOrig: Use original abDatasets from unmodified data to keep the majortrans and "Other" groups the same
#If the outputted abDatasets only needs to be for a subset of genes, specify these genes via the genestouse option
updateData <- function(x, cntGene, abGene, abDatasetsOrig = NULL, key, Group, seed, tx2gene, genestochange = NULL, useOtherGroups = TRUE, genestouse = NULL, samps = NULL, CompMI = FALSE,
                       DRIMSeqPowerAnalysis = FALSE){
  #To change back to a loop change x to changes[i]
  #Data names append the multiplicative factor that the counts are changed by to end
  #so abDatasets0.25 mults major trans by 0.25, abDatasets1 is just the observed data
  #since it involves mult the major trans counts by 1, etc
  NewData <- changeCountsAndTPM(x, cntGene = cntGene,
                                abGene = abGene, key = key, Group = Group, seedtoset = seed,
                                tx2gene = tx2gene, genestochange = genestochange, samps = samps)
  abGene <- NewData$abGene
  cntGene <- NewData$cntGene
  cntGenecntsScaledTPM <- NewData$cntGenecntsScaledTPM
  geneschanged <- NewData$geneschanged



  #Restrict to only genes that can be used in the compositional analyses via abDatasets and cntDatasets
  #ie only genes that have >1 trans and have TGE >0 across all samples
  #Use CompAbGene and CompCntGene to get abDatasets and cntDatasets and use CompCntGene in DRIMSeq Analysis to make results comparable

  CompAbGene <- subset(abGene, (abGene$NTrans!=1 & abGene$SumTGE!=0))
  CompCntGene <- subset(cntGene, (cntGene$NTrans!=1 & cntGene$SumTGE!=0))
  CompCntGenecntsScaledTPM <- subset(cntGenecntsScaledTPM, (cntGenecntsScaledTPM$NTrans!=1 & cntGenecntsScaledTPM$SumTGE!=0))

  fullgenenames <- sort(unique(CompAbGene$gene_id))

  if(is.null(genestouse)==TRUE){
    genestouse <- fullgenenames
  }
  
  
  if(DRIMSeqPowerAnalysis==TRUE){
    return(list(cntGenecntsScaledTPM = cntGenecntsScaledTPM, geneschanged = geneschanged))
  }else{
    abDatasetsNew <- lapply(genestouse, generateData, dat = CompAbGene,
                            nsamp = length(Group), abundance = TRUE, abData = CompAbGene,
                            abCompDatasets = abDatasetsOrig, useOtherGroups = useOtherGroups,
                            useExistingOtherGroups = TRUE, samps = samps, CompMI = CompMI)
    
    names(abDatasetsNew) <- genestouse
    
    return(list(CompAbGene = CompAbGene, CompCntGene = CompCntGene, abDatasets = abDatasetsNew, geneschanged = geneschanged))
  }
  #
  # abDatasetsNew <- laply(genestouse, generateData, dat = CompAbGene,
  #                        nsamp = length(Group), abundance = TRUE, abData = CompAbGene,
  #                        abCompDatasets = abDatasetsOrig, useOtherGroups = useOtherGroups,
  #                        useExistingOtherGroups = TRUE, samps = samps, .inform = TRUE)

  #Use the Existing Other Groups if you are using them - if not, useExistingOtherGroups will be ignored

  

  #Cnt Datasets is not currently used, don't calculate it for now
  # cntDatasets <- parLapply(clust, genestouse, generateData, dat = get(paste0("CompCntGene", x)),
  #                  nsamp = length(Group), abundance = FALSE, abData = get(paste0("CompAbGene", x)),
  #                  abCompDatasets = get(paste0("abDatasets", x)))
  # names(cntDatasets) <- genestouse
  
}




#q is df of condition variable
calcPillaiPval <- function(SigmaTildeNull, SigmaTildeAlt, res1 = NULL, q, nsamp, GibbsCovs = NULL, df_residual = NA, fullReturn = FALSE, returnTestStatOnly = FALSE){


  if(is.null(res1) & is.na(df_residual)){
    stop("df. residual must be specified to calcPillaiPval or else an lm object must specified in res1 to extract df.residual from")
  }

  Etilde <- nsamp * SigmaTildeAlt
  Htilde <- nsamp * (SigmaTildeNull - SigmaTildeAlt)

  vals <- tryCatch(diag((Htilde %*%solve(Etilde + Htilde))), error = function(x){})
  if(is.null(vals)){
    return(NA)
  }
  pill_stat <- sum(vals)


  #See the Multivariate ANOVA Testing pdf document (from the SAS help file) for the necessary formulas

  #v <- nsamp - ncol(model.matrix(res1))

  if(!is.na(df_residual)){
    v <- df_residual
  }else{
    v <- res1$df.residual
  }
  #v is the error/residual df- also extract from the r anova fit


  #p is the number of eigenvales (ie rank of Etilde + Htilde)
  p <- length(vals)

  s <- min(p,q)

  m <- 0.5 * (abs(p - q) - 1)
  n <- 0.5 * (v - p - 1)

  #Formulas come from the SAS help file Multivariate ANOVA Testing
  piece1 <- 2*n + s + 1
  piece2 <- 2*m + s + 1
  fstat_pillai <- (piece1/piece2) * (pill_stat/(s - pill_stat))

  if(fstat_pillai < 0){
    return(NA)
  }
  
  
  if(returnTestStatOnly==TRUE){
    return(fstat_pillai)
  }
  numdf_pillai <- s * piece2
  denomdf_pillai <- s * piece1
  pval_pillai <- 1-pf(fstat_pillai, df1 = numdf_pillai, df2 = denomdf_pillai)
  
  if(fullReturn==TRUE){
    return(list(pval_pillai = pval_pillai, pill_stat = pill_stat, df_residual = v, Etilde = Etilde, Htilde = Htilde, numdf_pillai = numdf_pillai, denomdf_pillai = denomdf_pillai))
  }else{
    return(pval_pillai)
  }
  
}




#This function runs for one biological sample at a Time, and thus needs to have the analysis repeated for each biological sample
  #Note that countsFromAbundance is not an argument because the Gibbs counts from Salmon don't change based on it
  #Set GibbsSamps to TRUE if the infReps were Gibbs Samples draws and FALSE if they're boostrap sample draws
SaveGibbsDataAsRData <- function(curr_samp, curr_file_loc, GibbsSamps = TRUE, countsFromAbundance = "no", direc_to_save_res = NULL){

    startTime <- proc.time()
    QuantSalmon <- tryCatch(tximport(curr_file_loc, type = "salmon", txOut = TRUE, countsFromAbundance = countsFromAbundance,
                                     ignoreTxVersion = FALSE, dropInfReps = F), error=function(e){})
    if(is.null(QuantSalmon)){
      stop("The Gibbs sampler seems to have failed for this sample, so no Gibbs output will be able to be saved")
    }


    t1 <- data.frame(QuantSalmon$infReps)

    #This could be nboot or ngibbs depending on whether gibbs samps or boot samps are used- just call it ninfreps here regardless
    ninfreps <- ncol(t1)
    colnames(t1) <- 1:ninfreps
    rownames(t1) <- rownames(QuantSalmon$abundance)

    t1$tx_id <- rownames(t1)

    #Add in gene and transcript information
    t2 <- merge(t1, tx2gene, by = "tx_id")
    t3 <- t2[order(t2$gene_id, t2$tx_id),]
    rownames(t3) <- t3$tx_id

    attr(t3, "Identifier") <- curr_samp

    #Want to move the tx_id column to be the last column- this is the easiest way to do this since the tx info is in the row names
    txcol <- which(colnames(t3) %in% "tx_id")
    t4 <- t3[,-txcol]
    t4$tx_id <- rownames(t3)

    #Extract length information and order it by gene/transcript
    l1 <- data.frame(QuantSalmon$length)
    colnames(l1) <- paste0(curr_samp, "Len")
    l1$tx_id <- rownames(l1)
    l2 <- merge(l1, tx2gene, by = "tx_id")

    l3 <- l2[order(l2$gene_id, l2$tx_id),]
    rownames(l3) <- l3$tx_id


    #Stack all gibbs/boot samples for this biological sample on top of each other
    outp <- rbindlist(apply(as.matrix(1:ninfreps), 1, SaveGibbsDataAsRDataHelper, t4 = t4, l3 = l3, curr_samp = curr_samp))
    #outp <- rbindlist(apply(as.matrix(1:1), 1, SaveGibbsDataAsRDataHelper, t4 = t4, l3 = l3, curr_samp = curr_samp))
    outp2 <- merge(outp, tx2gene, by = "tx_id")
    SaveInfRepsAsRCompTime <- proc.time() - startTime

  if(GibbsSamps==TRUE){
    if(is.null(direc_to_save_res)){
      save_dir <- paste0(getwd(), "/GibbsSamps/")
    }else{
      save_dir <- paste0(direc_to_save_res, "GibbsSamps/")
    }

    attr(outp2, "save_dir") <- save_dir
    attr(outp2, "ngibbs") <- ninfreps
    attr(outp2, "ninfreps") <- ninfreps
    attr(outp2, "Sample") <- curr_samp
    nam <- paste0("GibbsSamps", curr_samp)
    assign(nam, outp2)

    if(!dir.exists(save_dir)){
      dir.create(save_dir)
    }
    save(list = c(nam, "SaveInfRepsAsRCompTime", "countsFromAbundance"), file = paste0(save_dir, "GibbsSamps", curr_samp, ".RData" ))
  }else{
    if(is.null(direc_to_save_res)){
      save_dir <- paste0(getwd(), "/BootSamps/")
    }else{
      save_dir <- paste0(direc_to_save_res, "BootSamps/")
    }
    attr(outp2, "save_dir") <- save_dir
    attr(outp2, "nboot") <- ninfreps
    attr(outp2, "ninfreps") <- ninfreps
    attr(outp2, "Sample") <- curr_samp
    nam <- paste0("BootSamps", curr_samp)
    assign(nam, outp2)

    if(!dir.exists(save_dir)){
      dir.create(save_dir)
    }
    save(list = c(nam, "SaveInfRepsAsRCompTime", "countsFromAbundance"), file = paste0(save_dir, "BootSamps", curr_samp, ".RData" ))
  }
}

SaveGibbsDataAsRDataHelper <- function(x, t4, l3, curr_samp){
  cntcol <- as.character(x)
  colst5 <- c(cntcol, "gene_id", "NTrans", "tx_id")
  t5 <- t4[,colst5]

  colnames(t5)[colnames(t5)==cntcol] <- paste0(curr_samp, "Cnt")

  temp1 <- cntsToTPM(cnts = t5, nsamp = 1, len = l3, samps = curr_samp)
  temp1$infRepNum <- x

  t6 <- t5[,c("tx_id", paste0(curr_samp, "Cnt"))]
  temp2 <- merge(temp1, t6, by = "tx_id")
  return(temp2)

}



#investigateGibbsSamps is set to true only when looking at the covariance or the Gibbs samples to see how it work
  #Usually won't need this, just for SVB testing purposes
  #If you do have this, need to specify 1 biological sample at a time as invesSamp arg

calcIlrMeansCovs <- function(x, dat, nsamp, investigateGibbsSamps = FALSE, invesSamp = NA, CLE = FALSE){
  test1 <- dat[[x]]
  if(is.null(test1)==TRUE){
    return(NULL)
  }
  if(ncol(test1)==1){
    return(NULL)
  }
  test2 <- test1[mixedorder(rownames(test1)),]

  SampNames <- laply(rownames(test2), function(x){strsplit(x, "TPM")[[1]][1]})
  UniqueSampNames <- mixedsort(unique(SampNames))

  UniqueSampNumbers <- as.numeric(laply(UniqueSampNames, function(x){strsplit(x, "Sample")[[1]][2]}))
  #codetest <- apply(test2, 1, ilr)
  #Calculate the ilr values and their means and covs separately for each biological sample
  Means <- list()
  Covs <- list()

  if(investigateGibbsSamps==TRUE){
    vals <- invesSamp
  }else{
    vals <- UniqueSampNumbers
  }
  for (i in 1:length(vals)){
    #minrange <- (i-1)*100 + 1
    #maxrange <- i * 100

    curr_samp <- UniqueSampNames[i]
    #curr_samp <- SampNames[i]

    #test3 <- test2[minrange:maxrange,]
    #test3 <- test2[UniqueSampNames==curr_samp,]

    #Dimension of test3 is nsamp*ncol of the current abDataset file

    #Only select observations corresponding to the current sample
    if(length(SampNames)!=nrow(test2)){
      stop("Subsetting of test2 within calcilrMeansCovs has the incorrect number of dimensions, check it and try again")
    }
    test3 <- test2[SampNames==curr_samp,]
    
    test3SampNames <- laply(rownames(test3), function(x){strsplit(x, "TPM")[[1]][1]})
    if(length(unique(test3SampNames))!=1){
      stop("Something is wrong in calcIlrMeansCovs, observations from more than 1 sample are being used at the same time")
    }
    #If no rows, there are no Gibbs samples for this sample and we can move to the next sample
    if(nrow(test3)==0){
      Means[[i]] <- "No Gibbs/Boot Samples Available"
      Covs[[i]] <- "No Gibbs/Boot Samples Available"
      next
    }

    #Note that we can just apply the ilr transform on the TPM (and not on the relative proportion abundances) because the ilr values will be the same
      # in either case since the ilr of a normalized vector (or closed, which is normalilzed summing to 1) is the same as the ilr of the same vector
      # non-normalized.  Can confirm that test4T below is the same as test4 (except for small computation related error)

    # test3T <- ccomp(test3, total=1)
    # test4T <- apply(test3T, 1, ilr)

    #Dim of test 4 will be (ncol of current abDataset file - 1) * nsamp
      #This is not the final form, will be transposed by test5

    if(CLE==TRUE){
      test3_2 <- CorrectLowExpression(test3)
      test4 <- apply(test3_2, 1, function(x){ilr(x)})
    }else{
      test4 <- apply(test3, 1, function(x){ilr(x)})
    }
    


    #Dim of test 5 will be nsamp * (ncol of current abDataset file - 1), as intended
    if(is.null(ncol(test4))){
      test5 <- data.frame(as.matrix(test4, ncol = 1))
    }else{
      test5 <- data.frame(t(test4))
    }
    #Columns don't have a practical interpretation here, so change the names to null
    #See Analyzing Compositional Data with R book, p 45 for an explanation of this
    colnames(test5) <- NULL
    if(investigateGibbsSamps==TRUE){
      return(test5)
    }

    ilrMeans <- as.data.frame(t(as.matrix(colMeans(test5))))
    colnames(ilrMeans) <- NULL
    rownames(ilrMeans) <- NULL

    ilrCov <- cov(test5)

    Means[[i]] <- ilrMeans
    #names(Means)[i] <- paste0("Sample", i)
    Covs[[i]] <- ilrCov
    #names(Covs)[i] <- paste0("Sample", i)
  }
  #names(Means) <- paste0("Sample", 1:nsamp)
  names(Means) <- UniqueSampNames
  #names(Covs) <- paste0("Sample", 1:nsamp)
  names(Covs) <- UniqueSampNames
  #Means2 <- rbindlist(Means)
  #rownames(Means2) <- paste0("Sample", 1:nsamp)
  return(list(ilrMeans = Means, ilrCovs = Covs))
}

#Only specify key_full if you are using a random subset of the samples for the analysis
  #This is needed to extract the total number of samples, whcih is needed to be able to get all of the
  #Possible count columns of interest which will later be matched to the samples used for that specific analysis
PreDRIMSeq <- function(cntGene, key, failedgibbssamps = NULL, PreFiltered = FALSE){
  if(is.null(failedgibbssamps)==TRUE){
    failedgibbssamps <- integer(0)
  }

  # if(is.null(key_full)==FALSE){
  #   nsamp <- nrow(key_full)
  # }
  #Only include genes that have >1 trans, since only these can be used for differential transcript analysis
  if(PreFiltered == TRUE){
    cntGene2 <- cntGene
  }else{
    cntGene2 <- cntGene[(cntGene$NTrans > 1 & cntGene$SumTGE!=0),]
  }
  

  #Only using columns corresponding to counts that had gibbs replicates that worked for comparison with Compositional method
  #sampscols <- paste0("Sample", 1:nsamp, "Cnt")
  if(length(failedgibbssamps)==0){
    sampstouse <- key$Identifier
    sampstousecols <- paste0(sampstouse, "Cnt")
    print("First 10 samples to be used are given below")
    print(head(sampstouse, 10))
    cnts <- cntGene2[,c("tx_id", "gene_id", sampstousecols)]

    Group <- key$Condition[key$Identifier %in% sampstouse]

  }else{
    failedsamps <- paste0("Sample", failedgibbssamps)
    failedsampscols <- paste0(failedsamps, "Cnt")

    sampstouse <- key$Identifier[!(key$Identifier %in% failedsamps)]
    sampstousecols <- paste0(sampstouse, "Cnt")
    
    print("First 10 samples to be used are given below")
    print(head(sampstouse, 10))
    #sampstousecols <- sampscols[!(sampscols %in% failedsampscols)]
    cntGene3 <- cntGene2[,!(colnames(cntGene2) %in% failedsampscols)]

    #cntcolnums <- which(colnames(cntGene3) %in% sampstousecols)
    cnts <- cntGene3[,c("tx_id", "gene_id", sampstousecols)]

    Group <- key$Condition[key$Identifier %in% sampstouse]

  }

  colnames(cnts)[colnames(cnts) == "tx_id"] <- "feature_id"
  colnames(cnts)[colnames(cnts) %in% sampstousecols] <- sampstouse
  rownames(cnts) <- cnts$feature_id

  grp <- relevel(factor(Group), ref = 1)
  samp <- data.frame(sample_id = sampstouse, group = grp)

  return(list(cnts = cnts, samp = samp))
}


#x here is a group combination here
#Note that Filtering method here is often "none" because I prefilter the data to match my default filtering for
  # the compositional method and don't want to filter it again here
DRIMSeqObsRes <- function(cnts, samp, FilteringMethod = "none", fullReturn = FALSE, dataAlreadyFiltered = FALSE, Add1ToEveryCount = FALSE){
  #grp <- relevel(as.factor(AllGroupCombinations[i,]), ref = 1)
  #library(DRIMSeq)
  #start <- proc.time()
  library(DRIMSeq)
  if(colnames(cnts)[1]!= "feature_id" | colnames(cnts)[2]!= "gene_id"){
    stop("Ensure col of cnts argument input into DRIMSeqObsRes are in the right order with gene_id and feature_id as the first 2 columns")
  }

  StartTimeData <- proc.time()
  print(paste0("Add 1 to every count is ", Add1ToEveryCount))
  #browser()
  if(Add1ToEveryCount==TRUE){
    cnts2 <- cnts
    cnts2[,3:ncol(cnts2)] <-  cnts2[,3:ncol(cnts2)] + 1
    use_add_uniform <- FALSE
  }else{
    cnts2 <- cnts
    use_add_uniform <- TRUE
  }
  
  print(paste0("Use add_uniform is ", use_add_uniform))
  
  DRIMData <- dmDSdata(counts = cnts2, samples = samp)
  print("DRIMData object has been created")
  design_full <- model.matrix( ~ samp$group)
  print("The head of the samp object use din DRIMSeq's analysis is given below")
  print(head(samp))
  print("The head of the first 6 columns of the counts object used in the DRIMSeq analysis is given below")
  print(cnts2[1:5,1:6])
  print(paste0("The condition used in DRIMSeq's design matrix is below"))
  print(samp$group)
  print(paste0("Dimension of cnts object used on the DRIMSeq analysis is below"))
  print(dim(cnts))

  #set.seed(2323)
  #DRIMPres <- dmPrecision(DRIMData, design = design_full, common_precision = FALSE, BPPARAM = custBiocParallel)

  if(FilteringMethod == "DRIMSeq" & dataAlreadyFiltered==FALSE){
    n <- length(samp$group)
    n.small <- min(table(samp$group))
    DRIMData2 <- dmFilter(DRIMData,
                          min_samps_feature_expr=n.small, min_feature_expr=10,
                          min_samps_feature_prop=n.small, min_feature_prop=0.1,
                          min_samps_gene_expr=n, min_gene_expr=10)
  }else if(FilteringMethod == "none" | dataAlreadyFiltered==TRUE){
    DRIMData2 <- DRIMData
  }
  DRIMSeqDataTime <- proc.time() - StartTimeData
  
  StartRunTime <- proc.time()
  DRIMPres <- dmPrecision(DRIMData2, design = design_full, common_precision = TRUE, genewise_precision = TRUE, verbose = 1, add_uniform = use_add_uniform)
  print("dmPrecision step has completed")
  #DRIMFit <- dmFit(DRIMPres, design = design_full, verbose = 1, BPPARAM = custBiocParallel)
  DRIMFit <- dmFit(DRIMPres, design = design_full, bb_model = FALSE, add_uniform = use_add_uniform, verbose = 1)

  print("dmFit step has completed")
  design_null <- model.matrix( ~ 1, data = samp)

  #DRIMTest <- dmTest(DRIMFit, bb_model = TRUE, design = design_null, BPPARAM = custBiocParallel)
  DRIMTest <- dmTest(DRIMFit, bb_model = FALSE, design = design_null, verbose = 1)

  DRIMSeqRunTime <- proc.time() - StartRunTime
  DRIMres_gene <- results(DRIMTest, level = "gene")
  #DRIMres_feature <- results(DRIMTest, level = "feature")

  PvalsTemp <- DRIMres_gene[,c("pvalue", "adj_pvalue", "lr", "df")]
  #rownames(PvalsTemp) <- paste0(DRIMres_gene$gene_id, "Comb", i)
  rownames(PvalsTemp) <- paste0(DRIMres_gene$gene_id)


  # if(i==1){
  #   AllDRIMpvals <- PvalsTemp
  # }else{
  #   AllDRIMpvals <- rbind(AllDRIMpvals, PvalsTemp)
  # }
  #TimeThisGroupCombo <- proc.time() - start
  #return(list(Pvals = PvalsTemp, CompTime = TimeThisGroupCombo))

  if(fullReturn==TRUE){
    return(list(DRIMFit = DRIMFit, DRIMTest = DRIMTest, DRIMPvals = PvalsTemp, DRIMSeqDataTime = DRIMSeqDataTime, DRIMSeqRunTime = DRIMSeqRunTime))
  }else{
    return(PvalsTemp)
  }


}


DRIMSeqObservedAnalysis <- function(countsFromAbundance, FilteringMethod, genestouse, sampstouse = NULL, save_dir = NULL, fullReturn = TRUE){
  if(countsFromAbundance == "no" & FilteringMethod=="DRIMSeq"){
    #Results for Observed Data Only
    #load("cntGeneFiltered.RData")
    #cntGeneToUse <- cntGeneFiltered


    load("cntGeneFiltered.RData")
    cntGeneToUse <- cntGene
    # Load the list of Gibbs samples that failed because these won't be used in the compositional analysis and
    #this will make sure DRIMSeq is only used on the same samples
    #load("failedgibbssamps.RData")
    #failedgibbssamps <- failedgibbs
    failedgibbssamps <- NULL

  }else if(countsFromAbundance == "no" & FilteringMethod=="OtherGroups"){

    if(!file.exists("CntsOtherGroups.RData")){
      load("cntGeneFiltered.RData")
      load("cntDatasetsNoOtherGroupsFiltered.RData")

      CntsOtherGroups <- convertcntGeneToOtherGroupsForDRIMSeq(cntGene = cntGene, cntDatasets = cntDatasets)

      save(CntsOtherGroups, file = "CntsOtherGroups.RData")
      cntGeneToUse <- CntsOtherGroups
    }else{
      load("CntsOtherGroups.RData")
      cntGeneToUse <- CntsOtherGroups
    }

    # Load the list of Gibbs samples that failed because these won't be used in the compositional analysis and
    #this will make sure DRIMSeq is only used on the same samples
    #load("failedgibbssamps.RData")
    #failedgibbssamps <- failedgibbs
    failedgibbssamps <- NULL

  }else if(countsFromAbundance == "scaledTPM" & FilteringMethod=="DRIMSeq"){
    #load("cntGenecntsScaledTPMFiltered.RData")
    #cntGeneToUse <- cntGeneFiltered

    load("cntGenecntsScaledTPMFiltered.RData")
    cntGeneToUse <- cntGeneFiltered


    #load("failedgibbssampsCountsScaledTPM.RData")
    #failedgibbssamps <- failedgibbs

    failedgibbssamps <- NULL

  }else if(countsFromAbundance == "scaledTPM" & FilteringMethod=="OtherGroups"){
    if(!file.exists("CntsOtherGroupscntsScaledTPM.RData")){
      load("cntGenecntsScaledTPM.RData")
      load("cntDatasetscntsScaledTPM.RData")

      CntsOtherGroupscntsScaledTPM <- convertcntGeneToOtherGroupsForDRIMSeq(cntGene = cntGene, cntDatasets = cntDatasets)

      save(CntsOtherGroupscntsScaledTPM, file = "CntsOtherGroupscntsScaledTPM.RData")
      cntGeneToUse <- CntsOtherGroupscntsScaledTPM
    }else{
      load("CntsOtherGroupscntsScaledTPM.RData")
      cntGeneToUse <- CntsOtherGroupscntsScaledTPM
    }

    # Load the list of Gibbs samples that failed because these won't be used in the compositional analysis and
    #this will make sure DRIMSeq is only used on the same samples
      #With a newer version of Salmon this shouldn't be a problem as I don't think any of them failed anymore
    #load("failedgibbssampsCountsScaledTPM.RData")
    #failedgibbssamps <- failedgibbs
    failedgibbssamps <- NULL
  }
  cntGeneToUse2 <- subset(cntGeneToUse, cntGeneToUse$gene_id %in% genestouse)
  startTime <- proc.time()
  if(is.null(sampstouse)){
    key_to_use <- key
  }else{
    key_to_use <- subset(key, key$Identifier %in% sampstouse)
  }
  
  temp1 <- PreDRIMSeq(cntGene = cntGeneToUse2, failedgibbssamps = failedgibbssamps, key = key_to_use)
  cnts <- temp1$cnts
  samp <- temp1$samp
  #browser()
  if(is.null(sampstouse)){
    samp2 <- samp
    cnts_to_use <- cnts
  }else{
    samp2 <- subset(samp, samp$sample_id %in% sampstouse)
    samp2$group <- relevel(factor(samp2$group), ref = 1)
    cnts_to_use <- cnts[,c("feature_id", "gene_id", paste0(sampstouse))]
  }

  DRIMRes <- DRIMSeqObsRes(cnts = cnts_to_use, samp = samp2, FilteringMethod = FilteringMethod, fullReturn = fullReturn, dataAlreadyFiltered = TRUE, Add1ToEveryCount = TRUE)
  DRIMCompTime <- proc.time() - startTime
  DRIMFit <- DRIMRes$DRIMFit
  DRIMTest <- DRIMRes$DRIMTest
  DRIMPvals <- DRIMRes$DRIMPvals
  DRIMSeqDataTime <- DRIMRes$DRIMSeqDataTime
  DRIMSeqRunTime <- DRIMRes$DRIMSeqRunTime
  mem_info <- gc()
  if(!dir.exists("DRIMSeqObsRes")){
    dir.create("DRIMSeqObsRes")
  }

  if(!dir.exists("DRIMSeqObsResSubset")){
    dir.create("DRIMSeqObsResSubset")
  }

  if(is.null(save_dir)){
    if(countsFromAbundance == "no" & FilteringMethod=="OtherGroups" & is.null(sampstouse)){
      save_fil <- "DRIMSeqObsRes/DRIMSeqObsResOtherGroups.RData"
    }else if(countsFromAbundance == "scaledTPM" & FilteringMethod=="OtherGroups" & is.null(sampstouse)){
      save_fil <- "DRIMSeqObsRes/DRIMSeqObsRescntsScaledTPMOtherGroups.RData"
    }else if(countsFromAbundance == "no" & FilteringMethod=="DRIMSeq" & is.null(sampstouse)){
      save_fil <- "DRIMSeqObsRes/DRIMSeqObsResDRIMSeqFiltering.RData"
    }else if(countsFromAbundance == "scaledTPM" & FilteringMethod=="DRIMSeq" & is.null(sampstouse)){
      save_fil <- "DRIMSeqObsRes/DRIMSeqObsRescntsScaledTPMDRIMSeqFiltering.RData"
    }else if(countsFromAbundance == "no" & FilteringMethod=="OtherGroups" & !(is.null(sampstouse))){
      save_fil <- "DRIMSeqObsResSubset/DRIMSeqObsResOtherGroups.RData"
    }else if(countsFromAbundance == "scaledTPM" & FilteringMethod=="OtherGroups" & !(is.null(sampstouse))){
      save_fil <- "DRIMSeqObsResSubset/DRIMSeqObsRescntsScaledTPMOtherGroups.RData"
    }else if(countsFromAbundance == "no" & FilteringMethod=="DRIMSeq" & !(is.null(sampstouse))){
      save_fil <- "DRIMSeqObsResSubset/DRIMSeqObsResDRIMSeqFiltering.RData"
    }else if(countsFromAbundance == "scaledTPM" & FilteringMethod=="DRIMSeq" & !(is.null(sampstouse))){
      save_fil <- "DRIMSeqObsResSubset/DRIMSeqObsRescntsScaledTPMDRIMSeqFiltering.RData"
    }

    save(DRIMFit, DRIMTest, DRIMPvals, DRIMCompTime, DRIMSeqDataTime, DRIMSeqRunTime, mem_info, countsFromAbundance,
         file = save_fil)
  }else{
    if(!dir.exists(save_dir)){dir.create(save_dir)}
    if(FilteringMethod=="DRIMSeq"){
      # save(DRIMFit, DRIMTest, DRIMPvals, DRIMCompTime, DRIMSeqDataTime, DRIMSeqRunTime, mem_info, countsFromAbundance,
      #      file = paste0(save_dir, "DRIMSeqResDRIMSeqFiltering.RData"))
      save(DRIMPvals, DRIMCompTime, DRIMSeqDataTime, DRIMSeqRunTime, mem_info, countsFromAbundance,
           file = paste0(save_dir, "DRIMSeqResDRIMSeqFiltering.RData"))
    }else if(FilteringMethod=="OtherGroups"){
     # save(DRIMFit, DRIMTest, DRIMPvals, DRIMCompTime, DRIMSeqDataTime, DRIMSeqRunTime, mem_info, countsFromAbundance,
      #     file = paste0(save_dir, "DRIMSeqRes.RData"))
      save(DRIMPvals, DRIMCompTime, DRIMSeqDataTime, DRIMSeqRunTime, mem_info, countsFromAbundance,
           file = paste0(save_dir, "DRIMSeqRes.RData"))
    }

  }

}

#y is a Group combination here
DRIMSeqPower <- function(y, change, cntGene, abGene, cntDatasets, cntGeneFiltered, key, seed, tx2gene, failedgibbssamps = NULL,
                         genestochange = NULL, useOtherGroups = TRUE, Add1ToEveryCount = FALSE){
  
  samps <- key$Identifier
  startTime <- proc.time()
  print("The First 10 samples used in the analysis are printed below")
  print(head(samps, 10))
  print(paste0("Length of samps to use is ", length(samps)))
  if(sum(key$Condition!=y)!=0){stop("Check the specification of key and Group (y), these should be the same but do not match here")}
  NewData <- updateData(x = change, cntGene = cntGene, abGene = abGene, abDatasetsOrig = NULL,
                        key = key, Group = y, seed = seed, tx2gene = tx2gene, genestochange = genestochange,
                        DRIMSeqPowerAnalysis = TRUE, samps = samps)
  print("The Data has Been Updated.  The ScaledTPM counts will be used in the analysis.")
  cntGenecntsScaledTPM <- NewData$cntGenecntsScaledTPM
  if(!is.null(genestochange)){
    if(sum(NewData$geneschanged!=genestochange)!=0){
      stop("The genes modified for the power analysis should match genestochange but do not.")
    }
  }

  nsamp <- nrow(key)

  print(paste0("Nrow of cntGenecntsScaledTPM after power update is ", nrow(cntGenecntsScaledTPM)))

  #This function applies the other group conversion to the raw data to make it inline with the Compositional approach
    #ie it combines the appropriate transcripts into an other group to make sure the "filtering" is the same between the different approaches
    #Thus for these comparisons DRIMSeq pre built filters are not needed or used

  #CompCntGene <- subset(cntGenecntsScaledTPM, cntGenecntsScaledTPM$NTrans!=1 & cntGenecntsScaledTPM$SumTGE!=0)
  print(paste0("useOtherGroups is ", useOtherGroups))
  if(useOtherGroups==TRUE){
    CompCntGene2 <- convertcntGeneToOtherGroupsForDRIMSeq(cntGene = cntGenecntsScaledTPM, cntDatasets = cntDatasets, samps = samps)
  }else if(useOtherGroups==FALSE){
    #This gives a list of transcripts that have previously passed filtering
    CompCntGene2 <- subset(cntGenecntsScaledTPM, cntGenecntsScaledTPM$tx_id %in% cntGeneFiltered$tx_id)
  }
  if("gene_id" %in% colnames(CompCntGene2)){
    CompCntGene3 <- CompCntGene2
  }else{
    CompCntGene3 <- merge(CompCntGene2, tx2gene, by = "tx_id")
  }
  
  print(paste0("The number of transcripts used in this DRIMSeq Analysis is ", length(CompCntGene3$tx_id)))
  print(paste0("The number of genes used in this DRIMSeq Analysis is ", length(unique(CompCntGene3$gene_id))))


  #Generate a new key file corresponding to the new group permutation
    #This is only commented out in case it is needed later, but this is done in the individual power file as new_key so this isnt needed as is here
  # new_key <- data.frame(key$Identifier, y)
  # colnames(new_key) <- c("Identifier", "Condition")
  tempdrim <- PreDRIMSeq(cntGene = CompCntGene3, key = key, failedgibbssamps = failedgibbssamps, PreFiltered = TRUE)
  print("PreDRIMSeq is complete")
  #Only using columns corresponding to counts
  cnts <- tempdrim$cnts
  samp <- tempdrim$samp
  
  print(paste0("Ncol of cnts is ", ncol(cnts)))
  print(paste0("Ncol of samp is ", ncol(samp)))
  
  print(paste0("Nrow of cnts is ", nrow(cnts)))
  print(paste0("Nrow of samp is ", nrow(samp)))

  #Note that the filtering here is always set to none, as dataset will be filteredpreviously, 
    #whether through the use of Other Groups or by loading the prefiltered cntGene file based on DRIMSeq Filtering
  res <- DRIMSeqObsRes(cnts = cnts, samp = samp, FilteringMethod = "none", Add1ToEveryCount = Add1ToEveryCount)
  attr(res, "geneschanged") <- NewData$geneschanged
  return(list(res = res, Group = y, geneschanged = NewData$geneschanged, change = change))
}


calcOverlap2 <- function(x, tx2gene, exon_ranges){
  #Restrict analysis to a specific gene and its transcripts
  sub <- subset(tx2gene, tx2gene$gene_id==x)

  #Get sequences only from the transcripts associated with the current gene
  curr_exon_ranges <- exon_ranges[names(exon_ranges) %in% sub$tx_id]

  #res1 <- lapply(1:length(curr_exon_ranges), function(n) setdiff(curr_exon_ranges[n], unlist(curr_exon_ranges[-n])))
  res2 <- lapply(1:length(curr_exon_ranges), function(n) sum(width(setdiff(curr_exon_ranges[n], unlist(curr_exon_ranges[-n])))))

  tot_length_unique <- sum(as.numeric(res2))
  tot_length <- sum(width(unlist(curr_exon_ranges)))

  score <- tot_length_unique/tot_length

  #A score of 1 here would be no overlap, so return (1 - this value) to have 0 be no overlap and 1 be perfect overlap
  overlap <- 1-score
  return(overlap)
}

#Neat little function that will load an .RData object and save its contents to whatever as the results of the output
  #For example can do d <- loadRData("file.RData") and whatever is loaded is saved as d
loadRData <- function(fileName, objNameToGet = NULL){
  #loads an RData file, and returns it
  load(fileName)
  #print(ls()[ls() != "fileName"])
  if(is.null(objNameToGet)){
    rm(objNameToGet)
    #print(ls()[ls() != "fileName"])
    return(get(ls()[ls() != "fileName"]))
  }else{
    return(get(objNameToGet))
  }
  
}


#Extracts all pvalues for one group combination
#x is the group combination number
extractPowerPvals <- function(x, change, kind, direc, nparts = NULL, infReps = NULL, GEUV1Data = FALSE){
  if(x==1){
    print(paste0("Change value to be used within extractPowerPvals is ", change))
  }

  func1 <- function(x, CompMI = FALSE){
    d <- loadRData(x)
    if(CompMI==TRUE){
      #return(d)
      return(d$res)
    }else{
      return(d$res)
    }

  }

  if(kind == "DRIMSeq"){
    #setwd(paste0(def_wd, "SQCCDRIMSeqROCRes/", "Change", change))
    load(paste0(direc, "Change", change, "/DRIMSeqROCResGroupCombo", x, ".RData"))
    res <- get(paste0("DRIMSeqROCResGroupCombo", x))[[1]]
    res$gene_id <- rownames(res)
  }

  if(kind == "Comp"){
    #setwd(paste0(def_wd, "SQCCCompositionalMIROCRes/", "Change", change))
    if(!is.numeric(nparts)){
      re <- paste0(direc, "Change", change, "/CompROCResGroupCombo", x, ".RData")
      #print(paste0("file to load is ", re))
      if(file.exists(re)==TRUE){
        load(re)
        geneschanged <- attr(get(paste0("CompROCResGroupCombo", x))[[1]], "geneschanged")
        res <- get(paste0("CompROCResGroupCombo", x))[[1]]
      }
    }else{
      fils <- list()
      for (i in (1:nparts)){
        fils[[i]] <- paste0(direc, "Change", change, "/GroupCombo", x, "/CompROCResGroupCombo", x, "Part", i, ".RData")
      }
      reslist <- lapply(fils, func1)
      res <- rbindlist(reslist, fill = TRUE)
    }

  }

  if(kind == "CompMI"){
    fil <- paste0(direc, "Change", change, "/CompMIROCResGroupNum", x, ".RData")

    if(file.exists(fil)==FALSE){
      fil <- paste0(direc, "Change", change, "/CompMIROCResGroupCombo", x, ".RData")
    }


    if(file.exists(fil)==FALSE){
      fil <- paste0(direc, infReps, "/", "Change", change, "/CompMIROCResGroupNum", x, ".RData")
    }

    if(file.exists(fil)==FALSE){
      fil <- paste0(direc, infReps, "/", "Change", change, "/CompMIROCResGroupCombo", x, ".RData")
    }

    #If the directory Gibbs does not exist, it is because only Thin100 is used- so use this directory in this case
    if(file.exists(fil)==FALSE){
      fil <- paste0(direc, "GibbsThin100/", "Change", change, "/CompMIROCResGroupCombo", x, ".RData")
    }

    res <- func1(x = fil, CompMI = TRUE)
    #setwd(paste0(def_wd, "SQCCCompositionalMIROCRes/", "Change", change))
    # fils <- list()
    # for(i in 1:nparts){
    #   fils[[i]] <- paste0(direc, "Change", change, "/CompMIROCResGroupNum", i, ".RData")
    # }
    #
    # reslist <- lapply(fils, func1, CompMI = TRUE)
    # restemp <- rbindlist(reslist, fill = TRUE)
    # res <- subset(restemp, restemp$GroupNum==x)
  }

  if(kind == "CompMIThin100"){
      fil <- paste0(direc, "GibbsThin100/", "Change", change, "/CompMIROCResGroupNum", x, ".RData")

    if(file.exists(fil)==FALSE){
      fil <- paste0(direc, "GibbsThin100/", "Change", change, "/CompMIROCResGroupCombo", x, ".RData")
    }
    res <- func1(x = fil, CompMI = TRUE)
    #setwd(paste0(def_wd, "SQCCCompositionalMIROCRes/", "Change", change))
    # fils <- list()
    # for(i in 1:nparts){
    #   fils[[i]] <- paste0(direc, "Change", change, "/CompMIROCResGroupNum", i, ".RData")
    # }
    #
    # reslist <- lapply(fils, func1, CompMI = TRUE)
    # restemp <- rbindlist(reslist, fill = TRUE)
    # res <- subset(restemp, restemp$GroupNum==x)
  }

  if(kind == "CompMIThin16"){
    fil <- paste0(direc, "GibbsThin16/", "Change", change, "/CompMIROCResGroupNum", x, ".RData")

    if(file.exists(fil)==FALSE){
      fil <- paste0(direc, "GibbsThin16/", "Change", change, "/CompMIROCResGroupCombo", x, ".RData")
    }
    res <- func1(x = fil, CompMI = TRUE)
    #setwd(paste0(def_wd, "SQCCCompositionalMIROCRes/", "Change", change))
    # fils <- list()
    # for(i in 1:nparts){
    #   fils[[i]] <- paste0(direc, "Change", change, "/CompMIROCResGroupNum", i, ".RData")
    # }
    #
    # reslist <- lapply(fils, func1, CompMI = TRUE)
    # restemp <- rbindlist(reslist, fill = TRUE)
    # res <- subset(restemp, restemp$GroupNum==x)
  }

  if(kind == "RATs"){
    fil <- paste0(direc, "Change", change, "/RATsPowerResGroupCombo", x, ".RData")
    load(fil)

    restemp <- get(paste0("RATsPowerResGroupCombo", x))
    res <- restemp$res
    #setwd(paste0(def_wd, "SQCCCompositionalMIROCRes/", "Change", change))
    # fils <- list()
    # for(i in 1:nparts){
    #   fils[[i]] <- paste0(direc, "Change", change, "/CompMIROCResGroupNum", i, ".RData")
    # }
    #
    # reslist <- lapply(fils, func1, CompMI = TRUE)
    # restemp <- rbindlist(reslist, fill = TRUE)
    # res <- subset(restemp, restemp$GroupNum==x)
  }
  
  
  if(kind=="CompDTUMethod" | kind=="CompMIMethod"){
    if(GEUV1Data==TRUE){
      np <- 100
      ninf <- 100
    }else{
      np <- 10
      ninf <- 500
    }
    
    if(x==1){
      print(paste0("This code is assuming there are ", np, " parts for each group combination for the RanksRes and is always run on ", ninf, " inferential replicates.  Change this warning if this changes and ensure the code is handling this correctly"))
      print("In particular nparts and ninfreps are hard coded here")
    }
    # if(change!=4){
    #   stop("Only change 4 results have been run for the CompImputePermute.  If this has changed delete this line and ensure the code is properly handling the other change values.")
    # }
    
    resT <- vector(mode = "list", length = np)
    for(l in 1:np){
      #print(paste0("GroupCombo is ", x, " and Part is ", l))
      if(!file.exists(paste0(direc, "Change", change, "/GroupCombo", x, "/", "GroupCombo", x, "PermuteImputeRes", l, ".RData"))){
        stop(paste0("File does not exist for GroupCombo", x, " and Part ", l))
      }
      dat <- loadRData(paste0(direc, "Change", change, "/GroupCombo", x, "/", "GroupCombo", x, "PermuteImputeRes", l, ".RData"))
      if(!is.null(dat)){
        dat$gene_id <- rownames(dat)
      }
      resT[[l]] <- dat
    }
    resT2 <- data.table::rbindlist(resT, fill = TRUE)
    res <- data.frame(resT2)
    
    if(length(unique(res$gene_id))!=length(res$gene_id)){
      stop("The results have genes repeated across different parts and this shouldn't be the case.  Figure out why.")
    }
    
  }

  res$ComboID <- x
  return(res)

}


#c(seq(0.00000001, 0.0001, 0.00000025), seq(0.0001, 0.001, 0.0001), seq(0.001, 0.30, 0.001), 0.01, seq(0.30, 1, 0.01))
SVBROCCurves <- function(val, AllPvals, cutoffs = c(seq(0.00000001, 0.0001, 0.00000025), seq(0.0001, 0.001, 0.0001), seq(0.001, 0.30, 0.0005), 0.01, seq(0.30, 1, 0.01)), curr_methd){
  #load(paste0("CompositionalResAllGroupCombos", val, ".RData"))
  cutoffs <- as.matrix(cutoffs, ncol = 1)

  #assign(paste0("AllPvalsComp", val), ExtractAllPvals(get(paste0("CompAllGroupCombos", val)), Comp = TRUE, tx2gene = tx2gene, geneschanged = geneschanged))
  t2 <- data.frame(t(apply(cutoffs[1:nrow(cutoffs), , drop = FALSE], 1, calcSensFpr, AllPvals = AllPvals, curr_methd = curr_methd)))
  colnames(t2) <- c("sens", "fpr", "fdr", "acc", "sens_padjfdr", "fpr_padjfdr", "fdr_padjfdr", "acc_padjfdr", "pvalcutoff", "TotalNumCombosUsed", "MaxPossibleCombos")
  return(t2)

}

calcSensFpr <- function(x, AllPvals, curr_methd){

    AllPvals$reject <- as.numeric(AllPvals$p < x)
    AllPvals$reject_padjfdr <- as.numeric(AllPvals$padjfdr < x)


  # if(is.null(AllPvals$qvalues)){
  #   stop("Qvalues are Null, they are not being calculated and inputted to calcSensFpr properly")
  # }


  #rm(qvalsobj)
  #gc()


    Vals <- AllPvals


  ValsTruePos <- subset(Vals, (Vals$genechanged==1 & Vals$reject==1))
  ValsFalseNeg <- subset(Vals, (Vals$genechanged==1 & Vals$reject==0))

  ValsTrueNeg <- subset(Vals, (Vals$genechanged==0 & Vals$reject==0))
  ValsFalsePos <- subset(Vals, (Vals$genechanged==0 & Vals$reject==1))

  nTruePos <- nrow(ValsTruePos)
  nFalseNeg <- nrow(ValsFalseNeg)

  nTrueNeg <- nrow(ValsTrueNeg)
  nFalsePos <- nrow(ValsFalsePos)


  #Proportion of genes that have been changed that actually reject
  sens <- nTruePos/(nTruePos + nFalseNeg)

  #Proportion of non-changed genes that reject
  fpr <- nFalsePos/(nFalsePos + nTrueNeg)

  #False discovery rate- proportion of rejections that are mistakes
  fdr <- nFalsePos/(nFalsePos + nTruePos)

  #Accuracy- proportion of all tests that are correct
  acc <- (nTruePos + nTrueNeg)/(nTruePos + nFalsePos + nTrueNeg + nFalseNeg)


 #Now, results based on qvalues instead of pvalues

  ValsTruePos_padjfdr <- subset(Vals, (Vals$genechanged==1 & Vals$reject_padjfdr==1))
  ValsFalseNeg_padjfdr <- subset(Vals, (Vals$genechanged==1 & Vals$reject_padjfdr==0))

  ValsTrueNeg_padjfdr <- subset(Vals, (Vals$genechanged==0 & Vals$reject_padjfdr==0))
  ValsFalsePos_padjfdr <- subset(Vals, (Vals$genechanged==0 & Vals$reject_padjfdr==1))


  nTruePos_padjfdr <- nrow(ValsTruePos_padjfdr)
  nFalseNeg_padjfdr <- nrow(ValsFalseNeg_padjfdr)

  nTrueNeg_padjfdr <- nrow(ValsTrueNeg_padjfdr)
  nFalsePos_padjfdr <- nrow(ValsFalsePos_padjfdr)


  #Proportion of genes that have been changed that actually reject
  sens_padjfdr <- nTruePos_padjfdr/(nTruePos_padjfdr + nFalseNeg_padjfdr)

  #Proportion of non-changed genes that reject
  fpr_padjfdr <- nFalsePos_padjfdr/(nFalsePos_padjfdr + nTrueNeg_padjfdr)

  #False discovery rate- proportion of rejections that are mistakes
  fdr_padjfdr <- nFalsePos_padjfdr/(nFalsePos_padjfdr + nTruePos_padjfdr)


  #Accuracy- proportion of all tests that are correct
  acc_padjfdr <- (nTruePos_padjfdr + nTrueNeg_padjfdr)/(nTruePos_padjfdr + nFalsePos_padjfdr + nTrueNeg_padjfdr + nFalseNeg_padjfdr)

  #Total Gene/ComboID combinations that the results are based on
  totcombos <- nrow(Vals)

  maxcombos <- nrow(AllPvals)

  return(c(sens, fpr, fdr, acc, sens_padjfdr, fpr_padjfdr, fdr_padjfdr, acc_padjfdr, x, totcombos, maxcombos))
  #return(list(sens = sens, fpr = fpr))
}

#x is the method to use
SensFprRes <- function(x){
  curr_methd <- x
    print(paste0("Current method is ", curr_methd))
    val <- curr_change
    if(curr_methd=="DRIMSeq"){
      if(DRIMSeqFiltering==TRUE & TwentySamplesOnlyCase==TRUE){
        direc <- paste0(save_dir, "DRIMSeqROCResDRIMSeqFilteringTwentySamples/")
      }else if(DRIMSeqFiltering==TRUE & TwentySamplesOnlyCase==FALSE){
        direc <- paste0(save_dir, "DRIMSeqROCResDRIMSeqFiltering/")
      }else if(DRIMSeqFiltering==FALSE){
        direc <- paste0(save_dir, "DRIMSeqROCRes/")
      }
      kind <- "DRIMSeq"
      
      #The below method is run on the scaledTPM counts (since it was run after Aug2019)
    }else if(curr_methd=="DRIMSeqAdd1ToEveryCount"){
      if(DRIMSeqFiltering==TRUE & TwentySamplesOnlyCase==TRUE){
        direc <- paste0(save_dir, "DRIMSeqROCResDRIMSeqFilteringAdd1ToEveryCountTwentySamples/")
      }else if(DRIMSeqFiltering==TRUE & TwentySamplesOnlyCase==FALSE){
        direc <- paste0(save_dir, "DRIMSeqROCResDRIMSeqFilteringAdd1ToEveryCount/")
      }
      
      if(DRIMSeqFiltering==FALSE){
        stop("The directory is not set for DRIMSeqFiltering=FALSE and Add1ToEveryCount being TRUE.  Change this and ensure the code is working if this changes.")
      }
      kind <- "DRIMSeq"
      
    }else if(curr_methd %in% c("CompMI", "CompMICombineCoefs")){
      #CompMI results for the GEUV1 data were run in the ImputePermute format, for the SEQC data the other older format
      if(GEUV1Data==TRUE){
        if(DRIMSeqFiltering==TRUE & infReps=="Boot" & TwentySamplesOnlyCase==TRUE){
          direc <- paste0(pvals_save_dir, "CompDTUMethodsPowerAnalysisResBootDRIMSeqFilteringTwentySamples/ActualDataCompMI/")
        }else if(DRIMSeqFiltering==TRUE & infReps=="Boot" & TwentySamplesOnlyCase==FALSE){
          direc <- paste0(pvals_save_dir, "CompDTUMethodsPowerAnalysisResBootDRIMSeqFiltering/ActualDataCompMI/")
        }else if(DRIMSeqFiltering==TRUE & infReps=="Gibbs" & TwentySamplesOnlyCase==TRUE){
          direc <- paste0(pvals_save_dir, "CompDTUMethodsPowerAnalysisResGibbsDRIMSeqFilteringTwentySamples/ActualDataCompMI/")
        }else if(DRIMSeqFiltering==TRUE & infReps=="Gibbs" & TwentySamplesOnlyCase==FALSE){
          direc <- paste0(pvals_save_dir, "CompDTUMethodsPowerAnalysisResGibbsDRIMSeqFiltering/ActualDataCompMI/")
        }else{
          stop("Other options for the MI methods from the ImputeandPermute results have not been implemented properly yet.")
        }
        kind <- "CompMIMethod"
      }else{
        if(DRIMSeqFiltering==TRUE){
          direc <- paste0(save_dir, "CompositionalMIROCResDRIMSeqFiltering/")
        }else if(DRIMSeqFiltering==FALSE){
          direc <- paste0(save_dir, "CompositionalMIROCRes/")
        }
        kind <- "CompMI"
      }

    }else if(curr_methd=="RATs" | curr_methd=="RATsInfReps"){
      if(infReps=="Gibbs" | infReps=="GibbsThin100"){
        if(DRIMSeqFiltering==TRUE){
          direc <- paste0(save_dir, "RATsPowerResDRIMSeqFiltering/")
        }else if(DRIMSeqFiltering==FALSE){
          direc <- paste0(save_dir, "RATsPowerRes/")
        }
      }else if(infReps=="Boot"){
        if(GEUV1Data==TRUE & TwentySamplesOnlyCase==TRUE & DRIMSeqFiltering==TRUE){
          direc <- paste0(save_dir, "RATsPowerResBootDRIMSeqFilteringTwentySamples/")
        }else if(GEUV1Data==TRUE & TwentySamplesOnlyCase==FALSE & DRIMSeqFiltering==TRUE){
          direc <- paste0(save_dir, "RATsPowerResBootDRIMSeqFiltering/")
        }else if (GEUV1Data==TRUE & DRIMSeqFiltering==FALSE){
          stop("RATs results were never run on the GEUV1 data for OtherGroups, check argument specification")
        }else{
          if(DRIMSeqFiltering==TRUE){
            direc <- paste0(save_dir, "RATsPowerResBootDRIMSeqFiltering/")
          }else if(DRIMSeqFiltering==FALSE){
            direc <- paste0(save_dir, "RATsPowerResBoot/")
          }
        }
        
      }
      kind <- "RATs"
    }else if(curr_methd %in% c("CompDTU", "CompDTUme", "CompDTUAbRowMeanInfRepBoot")){
      
      AbFromInfRepsMethods <- c("CompDTUAbRowMeanInfRepBoot")
      
      if(curr_methd %in% AbFromInfRepsMethods){
        print("For this method, DRIMSeqFiltering is set in the file SummarizePowerResNewAug2019 to TRUE because the results have only been run for DRIMSeq filtering and infReps is set to Boot because
              these results do not depend on infReps such that the results will be located in the directory corresponding to BootDRIMSeqFiltering")
      }
      if(GEUV1Data==TRUE){
        print(paste0("Current method is one of the CompImputePermute ones"))
      }else{
        print(paste0("Current method is one of the CompImputePermute ones, which is only run using DRIMSeqFiltering with Bootstrap samples for now for the SEQC data.  If this changes remove/modify this print statement."))
      }
      print(paste0("The directory to load these pvalues from is hardcoded and does not change when changing save_dir or pvals_save_dir options.  Ensure this is what you want/expect"))
      print(paste0("These results only correspond to results run on Acutal Data (not simulated data)  Modify this if needed"))
      if(GEUV1Data==TRUE){
        #First, intercept the AbFromInfRepsMethods since those methods don't differ depending on infReps being Gibbs or Boot
        if(DRIMSeqFiltering==TRUE & (curr_methd %in% AbFromInfRepsMethods) & TwentySamplesOnlyCase==TRUE){
          direc <- paste0(pvals_save_dir, "CompDTUMethodsPowerAnalysisResCompDTUResForAbFromInfRepsTwentySamples/ActualDataCompDTUResForAbFromInfReps/")
        }else if(DRIMSeqFiltering==TRUE & (curr_methd %in% AbFromInfRepsMethods) & TwentySamplesOnlyCase==FALSE){
          direc <- paste0(pvals_save_dir, "CompDTUMethodsPowerAnalysisResCompDTUResForAbFromInfReps/ActualDataCompDTUResForAbFromInfReps/")
        }else if(DRIMSeqFiltering==TRUE & infReps=="Boot" & TwentySamplesOnlyCase==TRUE){
          direc <- paste0(pvals_save_dir, "CompDTUMethodsPowerAnalysisResBootDRIMSeqFilteringTwentySamples/ActualData/")
        }else if(DRIMSeqFiltering==TRUE & infReps=="Boot" & TwentySamplesOnlyCase==FALSE){
          direc <- paste0(pvals_save_dir, "CompDTUMethodsPowerAnalysisResBootDRIMSeqFiltering/ActualData/")
        }else if(DRIMSeqFiltering==TRUE & infReps=="Gibbs" & TwentySamplesOnlyCase==TRUE){
          direc <- paste0(pvals_save_dir, "CompDTUMethodsPowerAnalysisResGibbsDRIMSeqFilteringTwentySamples/ActualData/")
        }else if(DRIMSeqFiltering==TRUE & infReps=="Gibbs" & TwentySamplesOnlyCase==FALSE){
          direc <- paste0(pvals_save_dir, "CompDTUMethodsPowerAnalysisResGibbsDRIMSeqFiltering/ActualData/")
        }else if(DRIMSeqFiltering==FALSE & infReps=="Gibbs"){
          direc <- paste0(pvals_save_dir, "CompDTUMethodsPowerAnalysisRes/ActualData/")
        }else{
          stop("You have specified a CompImputePermute Result for a combination of inputs that has not yet been run.  If this changes remove this line and ensure the particular case is being handled correctly.")
        }
        
      }
      
      kind <- "CompDTUMethod"
    }

    print(paste0("Directory the pvalues will be loaded from is ", direc))
    if(dir.exists(direc)){
      print('The above directory exists')
    }else{
      stop("The above directory does not exist")
    }
    print(paste0("Current kind is ", kind))
    print(paste0("Current value of DRIMSeqFiltering is ", DRIMSeqFiltering))
    print(paste0("Current value of TwentySamplesOnlyCase is ", TwentySamplesOnlyCase))
    print(paste0("Current value of infReps is ", infReps))

    if(dir.exists(direc)==FALSE){
      stop("direc of power files does not exist, check things are specified correctly or else will need to modify SensFprRes function")
    }

    assign(paste0("AllPvals", curr_methd, val),
           rbindlist(lapply(1:ncombos, extractPowerPvals, change = val,
                            kind = kind, direc = direc, nparts = nparts,
                            infReps = infReps, GEUV1Data = GEUV1Data), fill = TRUE))
    #library(plyr)
    # assign(paste0("AllPvals", curr_methd, val),
    #        rbindlist(laply(1:ncombos, extractPowerPvals, change = val,
    #                         kind = kind, direc = direc, nparts = nparts, infReps = infReps, .inform = T), fill = TRUE))

    ValsTemp <- get(paste0("AllPvals", curr_methd, val))

    rm(list = c(paste0("AllPvals", curr_methd, val)))
    gc()

      ValsTemp$genechanged <- as.numeric(ValsTemp$gene_id %in% geneschanged)
    
    ValsTemp2 <- merge(ValsTemp, GeneNTrans, by = "gene_id")

    rm(ValsTemp)
    gc()
    print(gc())

    
   if(curr_methd=="CompMI" | curr_methd=="CompMIThin16" | curr_methd=="CompMIThin100"){
      ValsTemp2$p <- ValsTemp2$pvalt
    }else if(curr_methd=="CompMICombineCoefs" | curr_methd=="CompMICombineCoefsThin16" | curr_methd=="CompMICombineCoefsThin100"){
      if("pval_combineCoefs" %in% colnames(ValsTemp2)){
        ValsTemp2$p <- ValsTemp2$pval_combineCoefs
      }else{
        ValsTemp2$p <- ValsTemp2$pval_combineStats
      }
    }else if(curr_methd=="CompMI"){
      ValsTemp2$p <- ValsTemp2$pvalt
    }else if(curr_methd=="RATs"){
      #Modify RATs Pvalues to incorporate the "elig_fx" and "rep_reprod" post-hoc filters the paper recommends to make comparison as fair as possible to the full
      #recommendations from RATs 
      #This code/procedure was taken from Love (2018) "Swimdown" paper
      #In particular see the file "dtu_analysis" from swimdown code available at https://github.com/mikelove/swimdown/blob/master/dtu/dtu_analysis.R
      ValsTemp2$p <- ValsTemp2$pval
      ValsTemp2$p[!(ValsTemp2$elig_fx & ValsTemp2$rep_reprod)] <- 1
      ValsTemp2$padjfdr <- ValsTemp2$adj_pvalue
      ValsTemp2$padjfdr[!(ValsTemp2$elig_fx & ValsTemp2$rep_reprod)] <- 1
    }else if(curr_methd=="RATsInfReps"){
      #Modify RATs Pvalues to incorporate the "elig_fx" and "rep_reprod" post-hoc filters the paper recommends to make comparison as fair as possible to the full
      #recommendations from RATs 
      #This code/procedure was taken from Love (2018) "Swimdown" paper
      #In particular see the file "dtu_analysis" from swimdown code available at https://github.com/mikelove/swimdown/blob/master/dtu/dtu_analysis.R
      ValsTemp2$p <- ValsTemp2$pval_infReps
      ValsTemp2$p[!(ValsTemp2$elig_fx_infReps & ValsTemp2$rep_reprod_infReps & ValsTemp2$quant_reprod_infReps)] <- 1
      ValsTemp2$padjfdr <- ValsTemp2$adj_pval_infReps
      ValsTemp2$padjfdr[!(ValsTemp2$elig_fx_infReps & ValsTemp2$rep_reprod_infReps & ValsTemp2$quant_reprod_infReps)] <- 1
    }else if(curr_methd=="CompDTU"){
      ValsTemp2$p <- ValsTemp2$pval_CompDTU
    }else if(curr_methd=="CompDTUme"){
      ValsTemp2$p <- ValsTemp2$pval_CompDTUme
    }else if(curr_methd=="CompDTUAbRowMeanInfRepBoot"){
      ValsTemp2$p <- ValsTemp2$pval_CompDTUFromMeanBoot
    }else{
      stop("No pvalue found for the current method specification.  Check that the method specification is correct or add the proper pvalue to look for. ")
    }
    


    #Only use genes that are in genestouse
    Vals3T <- ValsTemp2[ValsTemp2$gene_id %in% genestouse,]

    rm(ValsTemp2)
    gc()

    Vals3T2 <- Vals3T
    rm(Vals3T)
    gc()
    print(gc())

    #Make sure NTrans exists properly in the Vals3 file
    if(is.null(Vals3T2$NTrans)){
      if(!is.null(Vals3T2$NTrans.x) & !is.null(Vals3T2$NTrans.y)){
        stopifnot(sum(Vals3T2$NTrans.x!=Vals3T2$NTrans.y)==0)
        Vals3T2$NTrans <- Vals3T2$NTrans.x
      }else if(!is.null(Vals3T2$NTrans.x) & is.null(Vals3T2$NTrans.y)){
        Vals3T2$NTrans <- Vals3T2$NTrans.x
      }else if (is.null(Vals3T2$NTrans.x) & !is.null(Vals3T2$NTrans.y)){
        Vals3T2$NTrans <- Vals3T2$NTrans.y
      }
    }

    Vals3T3 <- Vals3T2[order(Vals3T2$ComboID, Vals3T2$gene_id),]

    rm(Vals3T2)
    gc()

    #Note that the adjp value (based on fpr correction) has to be done this way because you really want to
      #correct within each group combination, not combining all pvalues across all group combinations
    if(kind=="RATs"){
      Vals3 <- Vals3T3
    }else{
      NewVals <- list()
      for(k in 1:ncombos){
        #print(paste0("k is ", k))
        ValsTTTT <- subset(Vals3T3, Vals3T3$ComboID==k)
        if(kind == "Comp" | kind == "DRIMSeq" | kind == "CompMI"){
          ValsTTTT$p[ValsTTTT$p==0] <- .Machine$double.eps
          ValsTTTT$padjfdr <- p.adjust(ValsTTTT$p, method = "fdr")
        }
        
        
        NewVals[[k]] <- ValsTTTT
      }

      Vals3 <- data.table::rbindlist(NewVals)
      rm(NewVals)
      rm(ValsTTTT)
      gc()
    }

    rm(Vals3T3)
    gc()

    #Save All Pvalues for a given change/method to be able to use them later
      AllPvals <- Vals3

      print(paste0("Current value of infReps is ", infReps))
      
      if(TwentySamplesOnlyCase==TRUE){
        dirrr_modif <- "TwentySamples"
      }else{
        dirrr_modif <- ""
      }
      
      if(infReps=="Gibbs"){
        if(DRIMSeqFiltering==TRUE){
          direc <- paste0(pvals_save_dir, "PowerResDRIMSeqFiltering", dirrr_modif, "/AllPvals/")
        }else if(DRIMSeqFiltering==FALSE){
          direc <- paste0(pvals_save_dir, "PowerRes", dirrr_modif, "/AllPvals/")
        }

      }else if(infReps=="Boot"){
        if(DRIMSeqFiltering==TRUE){
          direc <- paste0(pvals_save_dir, "PowerResBootDRIMSeqFiltering", dirrr_modif, "/AllPvals/")
        }else if(DRIMSeqFiltering==FALSE){
          direc <- paste0(pvals_save_dir, "PowerResBoot", dirrr_modif, "/AllPvals/")
        }
      }

      if(!dir.exists(direc)) {dir.create(direc, recursive = TRUE)}
      print(paste0("Directory in which AllPvals are saved is ", direc))
      if(TwentySamplesOnlyCase==TRUE){
        fil_modifier <- "TwentySamples"
      }else{
        fil_modifier <- ""
      }
      if(FilterGenes==TRUE){
          fil_name_to_save <- paste0(direc, "AllPvalsFilter", FilterCutoff, "Change", curr_change, curr_methd, fil_modifier, ".RData")
          save(AllPvals, file = fil_name_to_save, version = 2)
      }else{
          fil_name_to_save <-paste0(direc, "AllPvalsChange", curr_change, curr_methd, fil_modifier, ".RData")
          save(AllPvals, file = fil_name_to_save, version = 2)
      }
      print(paste0("File saved is ", fil_name_to_save))
      print(paste0("Saving of AllPvals is complete"))
      rm(AllPvals)
      gc()

    #browser()
    Vals31 <- subset(Vals3, Vals3$gene_id %in% genesSmallNTrans)
    Vals32 <- subset(Vals3, Vals3$gene_id %in% genesMediumNTrans)
    Vals33 <- subset(Vals3, Vals3$gene_id %in% genesLargeNTrans)

    assign(paste0("MissingPropSmallNTrans", curr_methd), mean(is.na(Vals31$p)))
    assign(paste0("MissingPropMedNTrans", curr_methd), mean(is.na(Vals32$p)))
    assign(paste0("MissingPropLargeNTrans", curr_methd), mean(is.na(Vals33$p)))

    ValsLowOverlap <- subset(Vals3, Vals3$gene_id %in% genesLowOverlap)
    ValsMedOverlap <- subset(Vals3, Vals3$gene_id %in% genesMedOverlap)
    ValsHighOverlap <- subset(Vals3, Vals3$gene_id %in% genesHighOverlap)

    ValsLowerThirdOverlap <- subset(Vals3, Vals3$gene_id %in% genesLowerThirdOverlap)
    ValsMidThirdOverlap <- subset(Vals3, Vals3$gene_id %in% genesMidThirdOverlap)
    ValsHighThirdOverlap <- subset(Vals3, Vals3$gene_id %in% genesHighThirdOverlap)

    ValsSmallNTransLowerThirdOverlap <- subset(Vals31, Vals31$gene_id %in% genesLowerThirdOverlap)
    ValsSmallNTransMidThirdOverlap <- subset(Vals31, Vals31$gene_id %in% genesMidThirdOverlap)
    ValsSmallNTransHighThirdOverlap <- subset(Vals31, Vals31$gene_id %in% genesHighThirdOverlap)

    ValsSmallNTransLowTGE <- subset(Vals31, Vals31$gene_id %in% genesLowTGE)
    ValsSmallNTransMediumTGE <- subset(Vals31, Vals31$gene_id %in% genesMediumTGE)
    ValsSmallNTransHighTGE <- subset(Vals31, Vals31$gene_id %in% genesHighTGE)

    ValsMediumNTransLowerThirdOverlap <- subset(Vals32, Vals32$gene_id %in% genesLowerThirdOverlap)
    ValsMediumNTransMidThirdOverlap <- subset(Vals32, Vals32$gene_id %in% genesMidThirdOverlap)
    ValsMediumNTransHighThirdOverlap <- subset(Vals32, Vals32$gene_id %in% genesHighThirdOverlap)

    ValsLargeNTransLowerThirdOverlap <- subset(Vals33, Vals33$gene_id %in% genesLowerThirdOverlap)
    ValsLargeNTransMidThirdOverlap <- subset(Vals33, Vals33$gene_id %in% genesMidThirdOverlap)
    ValsLargeNTransHighThirdOverlap <- subset(Vals33, Vals33$gene_id %in% genesHighThirdOverlap)

    ValsLowTGE <- subset(Vals3, Vals3$gene_id %in% genesLowTGE)
    ValsMediumTGE <- subset(Vals3, Vals3$gene_id %in% genesMediumTGE)
    ValsHighTGE <- subset(Vals3, Vals3$gene_id %in% genesHighTGE)

    ValsGenes2Columns <- subset(Vals3, Vals3$gene_id %in% genes2Columns)
    ValsGenesWithNoOtherGroup <- subset(Vals3, Vals3$gene_id %in% genesWithNoOtherGroup)
    ValsGenesWithNoOtherGroupAnd2Columns <- subset(Vals3, Vals3$gene_id %in% genesWithNoOtherGroupAnd2Columns)

    ValsGenes2ColumnsHighTGE <- subset(Vals3, Vals3$gene_id %in% genes2ColumnsHighTGE)
    ValsGenesWithNoOtherGroupHighTGE <- subset(Vals3, Vals3$gene_id %in% genesWithNoOtherGroupHighTGE)
    ValsGenesWithNoOtherGroupAnd2ColumnsHighTGE <- subset(Vals3, Vals3$gene_id %in% genesWithNoOtherGroupAnd2ColumnsHighTGE)

    ValsGenesWithNoOtherGroupSmallNTrans <- subset(Vals3, Vals3$gene_id %in% genesWithNoOtherGroupSmallNTrans)
    ValsGenesWithNoOtherGroupMediumNTrans <- subset(Vals3, Vals3$gene_id %in% genesWithNoOtherGroupMediumNTrans)
    ValsGenesWithNoOtherGroupLargeNTrans <- subset(Vals3, Vals3$gene_id %in% genesWithNoOtherGroupLargeNTrans)

    ValsGenes2ColumnsHighOverlap <- subset(Vals3, Vals3$gene_id %in% genes2ColumnsHighOverlap)
    ValsGenesWithNoOtherGroupsHighOverlap <- subset(Vals3, Vals3$gene_id %in% genesWithNoOtherGroupHighOverlap)
    ValsGenesWithNoOtherGroupsAnd2ColumnsHighOverlap <- subset(Vals3, Vals3$gene_id %in% genesWithNoOtherGroupAnd2ColumnsHighOverlap)

    ValsGenes2ColumnsHighThirdOverlap <- subset(Vals3, Vals3$gene_id %in% genes2ColumnsHighThirdOverlap)
    ValsGenesWithNoOtherGroupsHighThirdOverlap <- subset(Vals3, Vals3$gene_id %in% genesWithNoOtherGroupHighThirdOverlap)
    ValsGenesWithNoOtherGroupsAnd2ColumnsHighThirdOverlap <- subset(Vals3, Vals3$gene_id %in% genesWithNoOtherGroupAnd2ColumnsHighThirdOverlap)

    ValsGenes3Columns <- subset(Vals3, Vals3$gene_id %in% genes3Columns)
    ValsGenesGr3Columns <- subset(Vals3, Vals3$gene_id %in% genesGr3Columns)
    
    
    if(useInfRV==TRUE){
      ValsgenesMeanGeneMaxInfRVTPMLow <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMaxInfRVTPMLow)
      ValsgenesMeanGeneMaxInfRVTPMMedium <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMaxInfRVTPMMedium)
      ValsgenesMeanGeneMaxInfRVTPMHigh <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMaxInfRVTPMHigh)
      
      ValsgenesMeanGeneMedianInfRVTPMLow <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMedianInfRVTPMLow)
      ValsgenesMeanGeneMedianInfRVTPMMedium <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMedianInfRVTPMMedium)
      ValsgenesMeanGeneMedianInfRVTPMHigh <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMedianInfRVTPMHigh)
      
      ValsgenesMeanGeneMeanInfRVTPMLow <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMeanInfRVTPMLow)
      ValsgenesMeanGeneMeanInfRVTPMMedium <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMeanInfRVTPMMedium)
      ValsgenesMeanGeneMeanInfRVTPMHigh <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMeanInfRVTPMHigh)
      
      ValsgenesMeanGeneMaxInfRVTPMEighty <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMaxInfRVTPMEighty)
      ValsgenesMeanGeneMaxInfRVTPMGrEighty <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMaxInfRVTPMGrEighty)
      ValsgenesMeanGeneMaxInfRVTPMGrNinety <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMaxInfRVTPMGrNinety)
      
      ValsgenesMeanGeneMeanInfRVTPMEighty <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMeanInfRVTPMEighty)
      ValsgenesMeanGeneMeanInfRVTPMGrEighty <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMeanInfRVTPMGrEighty)
      ValsgenesMeanGeneMeanInfRVTPMGrNinety <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMeanInfRVTPMGrNinety)
      
      ValsgenesMeanGeneMedianInfRVTPMEighty <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMedianInfRVTPMEighty)
      ValsgenesMeanGeneMedianInfRVTPMGrEighty <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMedianInfRVTPMGrEighty)
      ValsgenesMeanGeneMedianInfRVTPMGrNinety <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMedianInfRVTPMGrNinety)
      
      
      
      ValsgenesMeanGeneMaxVarilrScaleLow <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMaxVarilrScaleLow)
      ValsgenesMeanGeneMaxVarilrScaleMedium <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMaxVarilrScaleMedium)
      ValsgenesMeanGeneMaxVarilrScaleHigh <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMaxVarilrScaleHigh)
      
      
      ValsgenesMeanGeneMaxVarilrScaleGr80 <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMaxVarilrScaleGr80)
      ValsgenesMeanGeneMaxVarilrScaleGr90 <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMaxVarilrScaleGr90)
      
      
      
      
      ValsgenesMeanGeneMeanVarilrScaleLow <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMeanVarilrScaleLow)
      ValsgenesMeanGeneMeanVarilrScaleMedium <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMeanVarilrScaleMedium)
      ValsgenesMeanGeneMeanVarilrScaleHigh <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMeanVarilrScaleHigh)
      
      
      ValsgenesMeanGeneMeanVarilrScaleGr80 <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMeanVarilrScaleGr80)
      ValsgenesMeanGeneMeanVarilrScaleGr90 <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMeanVarilrScaleGr90)
      
      
      
      ValsgenesMeanGeneMedianVarilrScaleLow <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMedianVarilrScaleLow)
      ValsgenesMeanGeneMedianVarilrScaleMedium <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMedianVarilrScaleMedium)
      ValsgenesMeanGeneMedianVarilrScaleHigh <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMedianVarilrScaleHigh)
      
      
      ValsgenesMeanGeneMedianVarilrScaleGr80 <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMedianVarilrScaleGr80)
      ValsgenesMeanGeneMedianVarilrScaleGr90 <- subset(Vals3, Vals3$gene_id %in% genesMeanGeneMedianVarilrScaleGr90)
    }
    
    
    if(useInfRV==TRUE){
  
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMaxVarilrScaleLow"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMaxVarilrScaleLow, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMaxVarilrScaleMedium"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMaxVarilrScaleMedium, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMaxVarilrScaleHigh"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMaxVarilrScaleHigh, curr_methd = curr_methd))
      
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMaxVarilrScaleGr80"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMaxVarilrScaleGr80, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMaxVarilrScaleGr90"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMaxVarilrScaleGr90, curr_methd = curr_methd))
      
      
      
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMeanVarilrScaleLow"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMeanVarilrScaleLow, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMeanVarilrScaleMedium"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMeanVarilrScaleMedium, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMeanVarilrScaleHigh"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMeanVarilrScaleHigh, curr_methd = curr_methd))
      
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMeanVarilrScaleGr80"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMeanVarilrScaleGr80, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMeanVarilrScaleGr90"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMeanVarilrScaleGr90, curr_methd = curr_methd))
      
      
      
      
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMedianVarilrScaleLow"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMedianVarilrScaleLow, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMedianVarilrScaleMedium"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMedianVarilrScaleMedium, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMedianVarilrScaleHigh"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMedianVarilrScaleHigh, curr_methd = curr_methd))
      
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMedianVarilrScaleGr80"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMedianVarilrScaleGr80, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMedianVarilrScaleGr90"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMedianVarilrScaleGr90, curr_methd = curr_methd))
      
      
      
      
      
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMaxInfRVTPMLow"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMaxInfRVTPMLow, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMaxInfRVTPMMedium"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMaxInfRVTPMMedium, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMaxInfRVTPMHigh"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMaxInfRVTPMHigh, curr_methd = curr_methd))
      
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMedianInfRVTPMLow"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMedianInfRVTPMLow, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMedianInfRVTPMMedium"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMedianInfRVTPMMedium, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMedianInfRVTPMHigh"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMedianInfRVTPMHigh, curr_methd = curr_methd))
      
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMeanInfRVTPMLow"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMeanInfRVTPMLow, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMeanInfRVTPMMedium"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMeanInfRVTPMMedium, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMeanInfRVTPMHigh"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMeanInfRVTPMHigh, curr_methd = curr_methd))
      
      
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMaxInfRVTPMEighty"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMaxInfRVTPMEighty, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMaxInfRVTPMGrEighty"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMaxInfRVTPMGrEighty, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMaxInfRVTPMGrNinety"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMaxInfRVTPMGrNinety, curr_methd = curr_methd))
      
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMedianInfRVTPMEighty"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMedianInfRVTPMEighty, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMedianInfRVTPMGrEighty"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMedianInfRVTPMGrEighty, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMedianInfRVTPMGrNinety"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMedianInfRVTPMGrNinety, curr_methd = curr_methd))
      
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMeanInfRVTPMEighty"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMeanInfRVTPMEighty, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMeanInfRVTPMGrEighty"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMeanInfRVTPMGrEighty, curr_methd = curr_methd))
      assign(paste0("SensFpr", curr_methd, val, "genesMeanGeneMeanInfRVTPMGrNinety"), SVBROCCurves(val, AllPvals = ValsgenesMeanGeneMeanInfRVTPMGrNinety, curr_methd = curr_methd))
      
    }
    

    assign(paste0("SensFpr", curr_methd, val, "Genes3Columns"), SVBROCCurves(val, AllPvals = ValsGenes3Columns, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "GenesGr3Columns"), SVBROCCurves(val, AllPvals = ValsGenesGr3Columns, curr_methd = curr_methd))
    
    assign(paste0("SensFpr", curr_methd, val), SVBROCCurves(val, AllPvals = Vals3, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "SmallNTrans"), SVBROCCurves(val, AllPvals = Vals31, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "MediumNTrans"), SVBROCCurves(val, AllPvals = Vals32, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "LargeNTrans"), SVBROCCurves(val, AllPvals = Vals33, curr_methd = curr_methd))


    assign(paste0("SensFpr", curr_methd, val, "SmallNTransLowerThirdOverlap"), SVBROCCurves(val, AllPvals = ValsSmallNTransLowerThirdOverlap, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "SmallNTransMidThirdOverlap"), SVBROCCurves(val, AllPvals = ValsSmallNTransMidThirdOverlap, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "SmallNTransHighThirdOverlap"), SVBROCCurves(val, AllPvals = ValsSmallNTransHighThirdOverlap, curr_methd = curr_methd))


    assign(paste0("SensFpr", curr_methd, val, "MediumNTransLowerThirdOverlap"), SVBROCCurves(val, AllPvals = ValsMediumNTransLowerThirdOverlap, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "MediumNTransMidThirdOverlap"), SVBROCCurves(val, AllPvals = ValsMediumNTransMidThirdOverlap, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "MediumNTransHighThirdOverlap"), SVBROCCurves(val, AllPvals = ValsMediumNTransHighThirdOverlap, curr_methd = curr_methd))


    assign(paste0("SensFpr", curr_methd, val, "LargeNTransLowerThirdOverlap"), SVBROCCurves(val, AllPvals = ValsLargeNTransLowerThirdOverlap, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "LargeNTransMidThirdOverlap"), SVBROCCurves(val, AllPvals = ValsLargeNTransMidThirdOverlap, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "LargeNTransHighThirdOverlap"), SVBROCCurves(val, AllPvals = ValsLargeNTransHighThirdOverlap, curr_methd = curr_methd))

    assign(paste0("SensFpr", curr_methd, val, "SmallNTransLowTGE"), SVBROCCurves(val, AllPvals = ValsSmallNTransLowTGE, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "SmallNTransMediumTGE"), SVBROCCurves(val, AllPvals = ValsSmallNTransMediumTGE, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "SmallNTransHighTGE"), SVBROCCurves(val, AllPvals = ValsSmallNTransHighTGE, curr_methd = curr_methd))

    assign(paste0("SensFpr", curr_methd, val, "LowOverlap"), SVBROCCurves(val, AllPvals = ValsLowOverlap, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "MedOverlap"), SVBROCCurves(val, AllPvals = ValsMedOverlap, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "HighOverlap"), SVBROCCurves(val, AllPvals = ValsHighOverlap, curr_methd = curr_methd))

    assign(paste0("SensFpr", curr_methd, val, "LowerThirdOverlap"), SVBROCCurves(val, AllPvals = ValsLowerThirdOverlap, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "MidThirdOverlap"), SVBROCCurves(val, AllPvals = ValsMidThirdOverlap, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "HighThirdOverlap"), SVBROCCurves(val, AllPvals = ValsHighThirdOverlap, curr_methd = curr_methd))

    assign(paste0("SensFpr", curr_methd, val, "LowTGE"), SVBROCCurves(val, AllPvals = ValsLowTGE, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "MediumTGE"), SVBROCCurves(val, AllPvals = ValsMediumTGE, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "HighTGE"), SVBROCCurves(val, AllPvals = ValsHighTGE, curr_methd = curr_methd))

    assign(paste0("SensFpr", curr_methd, val, "Genes2Columns"), SVBROCCurves(val, AllPvals = ValsGenes2Columns, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "GenesWithNoOtherGroup"), SVBROCCurves(val, AllPvals = ValsGenesWithNoOtherGroup, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "GenesWithNoOtherGroupAnd2Columns"), SVBROCCurves(val, AllPvals = ValsGenesWithNoOtherGroupAnd2Columns, curr_methd = curr_methd))

    assign(paste0("SensFpr", curr_methd, val, "Genes2ColumnsHighTGE"), SVBROCCurves(val, AllPvals = ValsGenes2ColumnsHighTGE, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "GenesWithNoOtherGroupHighTGE"), SVBROCCurves(val, AllPvals = ValsGenesWithNoOtherGroupHighTGE, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "GenesWithNoOtherGroupAnd2ColumnsHighTGE"), SVBROCCurves(val, AllPvals = ValsGenesWithNoOtherGroupAnd2ColumnsHighTGE, curr_methd = curr_methd))

    assign(paste0("SensFpr", curr_methd, val, "GenesWithNoOtherGroupSmallNTrans"), SVBROCCurves(val, AllPvals = ValsGenesWithNoOtherGroupSmallNTrans, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "GenesWithNoOtherGroupMediumNTrans"), SVBROCCurves(val, AllPvals = ValsGenesWithNoOtherGroupMediumNTrans, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "GenesWithNoOtherGroupLargeNTrans"), SVBROCCurves(val, AllPvals = ValsGenesWithNoOtherGroupLargeNTrans, curr_methd = curr_methd))

    assign(paste0("SensFpr", curr_methd, val, "Genes2ColumnsHighOverlap"), SVBROCCurves(val, AllPvals = ValsGenes2ColumnsHighOverlap, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "GenesWithNoOtherGroupsHighOverlap"), SVBROCCurves(val, AllPvals = ValsGenesWithNoOtherGroupsHighOverlap, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "GenesWithNoOtherGroupsAnd2ColumnsHighOverlap"), SVBROCCurves(val, AllPvals = ValsGenesWithNoOtherGroupsAnd2ColumnsHighOverlap, curr_methd = curr_methd))

    assign(paste0("SensFpr", curr_methd, val, "Genes2ColumnsHighThirdOverlap"), SVBROCCurves(val, AllPvals = ValsGenes2ColumnsHighThirdOverlap, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "GenesWithNoOtherGroupsHighThirdOverlap"), SVBROCCurves(val, AllPvals = ValsGenesWithNoOtherGroupsHighThirdOverlap, curr_methd = curr_methd))
    assign(paste0("SensFpr", curr_methd, val, "GenesWithNoOtherGroupsAnd2ColumnsHighThirdOverlap"), SVBROCCurves(val, AllPvals = ValsGenesWithNoOtherGroupsAnd2ColumnsHighThirdOverlap, curr_methd = curr_methd))
    



    #Remove files that are not needed and that will just end up taking up memory
    rm(list = c(ls(pattern = "AllPvals"), ls(pattern = "ValsTemp")))

  #assign(paste0("AllPvals", curr_methd, val), Vals3)
  #I wanted to return all pvalues this way, but it takes up way too much memory so not practical
  #objs <- c(ls(pattern = "SensFpr"), paste0("AllPvals", curr_methd, val))

  objs <- c(ls(pattern = "SensFpr"))
  a1 <- mget(objs, inherits = TRUE)
  names(a1) <- objs
  return(a1)
}





SummarizePowerRes <- function(def_wd, save_dir, data_name, curr_change, curr_methd, geneschanged, nparts, ncombos, tx2gene = NULL, FilterGenes = TRUE, FilterCutoff = 0.50, useOtherGroups = TRUE, SQCCData = FALSE, pvals_save_dir = NULL, cntsScaledTPM = FALSE, infReps = NULL, DRIMSeqFiltering = FALSE, NumYsDat, useInfRV = FALSE, GEUV1Data = FALSE, infRVDRIMSeqFilteringBoot = FALSE, minimum_possible_pval_IP = (1/10000), TwentySamplesOnlyCase = FALSE){
  setwd(def_wd)
  if(is.null(tx2gene)){
    load("tx2gene.RData")
  }

  if(is.null(pvals_save_dir)){
    saveAllPvals <- FALSE
  }else{
    saveAllPvals <- TRUE
  }

  if(is.null(infReps)){
    infReps <- "Gibbs"
    print(paste0("infReps argument was not specified so it was set to a default value of Gibbs"))
  }

  if(DRIMSeqFiltering==TRUE){
    if(cntsScaledTPM==TRUE){
      load("cntGenecntsScaledTPMFiltered.RData")
      cntGene <- cntGeneFiltered
    }else{
      load("cntGeneFiltered.RData")
    }

    load("abGeneFiltered.RData")
    load("abDatasetsNoOtherGroupsFiltered.RData")

    abGene <- abGeneFiltered
    abDatasets <- abDatasetsFiltered
  }else if(DRIMSeqFiltering==FALSE){
    if(cntsScaledTPM==TRUE){
      load("cntGenecntsScaledTPM.RData")
    }else{
      load("cntGene.RData")
    }

    load("abGene.RData")
    load("abDatasets.RData")
  }

  load("OverlapRes.RData")




  #Only factor in genes from abDatasets into power analysis since only these were run using the compositional method
  genenames <- names(abDatasets)

  #Load the necessary files for the Power analysis
  #load("FilesForPowerAnalysis1.RData")

  print(paste0("FilterGenes is ", FilterGenes))
  if(FilterGenes==TRUE){
    #This will give a vector
    print(paste0("FilterGenes value is ", FilterCutoff))
    genestouset <- unique(abGene$gene_id[abGene$MeanTGE > FilterCutoff & abGene$NTrans > 1])
    genestouse <- genestouset[genestouset %in% genenames]
  }else{
    genestouse <- genenames
  }
  
  print(paste0("Length of genestouse is ", length(genestouse)))

  #Generate GeneNTrans, a data frame that has the number of transcripts for a given gene from tx2gene
  GeneNTransTemp <- data.frame(table(tx2gene$gene_id), stringsAsFactors = FALSE)
  colnames(GeneNTransTemp) <- c("gene_id", "NTrans")

  #Restrict GeneNTrans to only be from genes that will actually be used in the analysis
  GeneNTrans <- subset(GeneNTransTemp, GeneNTransTemp$gene_id %in% genestouse)


  overlapres$gene_id <- rownames(overlapres)

  transQuan <- quantile(GeneNTrans$NTrans, probs = c((1/3), (2/3), 1))
  NTranscutoff1 <- transQuan[1]
  NTranscutoff2 <- transQuan[2]

  genesSmallNTrans <- as.character(subset(GeneNTrans$gene_id,
                                          GeneNTrans$NTrans >=2 & GeneNTrans$NTrans <= NTranscutoff1))
  genesMediumNTrans <- as.character(subset(GeneNTrans$gene_id, GeneNTrans$NTrans > NTranscutoff1 & GeneNTrans$NTrans <= NTranscutoff2))
  genesLargeNTrans <- as.character(subset(GeneNTrans$gene_id, GeneNTrans$NTrans > NTranscutoff2))


  #Extract List of genes that have exactly 2 columns after collapsing to other groups/filtering
    #ie removing transcripts in the case of DRIMSeq filtering
  Ncolumns <- laply(abDatasets, function(x){dim(x)}[2])
  NcolumnsDatF <- data.frame(names(abDatasets), Ncolumns, stringsAsFactors = F)

  colnames(NcolumnsDatF) <- c("gene_id", "Ncolumns")

  genes2Columns <- as.character(subset(NcolumnsDatF$gene_id, NcolumnsDatF$Ncolumns==2))
  genes3Columns <- as.character(subset(NcolumnsDatF$gene_id, NcolumnsDatF$Ncolumns==3))
  genesGr3Columns <- as.character(subset(NcolumnsDatF$gene_id, NcolumnsDatF$Ncolumns > 3))

  #Genes that do not have any other groups despite them being used
  GenesWithNoOtherGroup <- laply(abDatasets, function(x){(length(attr(x, "OtherTrans"))==0) & !is.null(dim(x))})
  GenesWithNoOtherGroupDatF <- data.frame(names(abDatasets), GenesWithNoOtherGroup, stringsAsFactors = F)
  colnames(GenesWithNoOtherGroupDatF) <- c("gene_id", "GeneHasNoOtherGroup")

  genesWithNoOtherGroup <- as.character(subset(GenesWithNoOtherGroupDatF$gene_id, GenesWithNoOtherGroupDatF$GeneHasNoOtherGroup))

  genesWithNoOtherGroupAnd2Columns <- intersect(genes2Columns, genesWithNoOtherGroup)

  #Names of genes done by overlap at logical overlap cutoffs of 1/3, 2/3, and 1
  genesLowOverlapT <- overlapres$gene_id[overlapres$overlap <= (1/3)]
  genesMedOverlapT <- overlapres$gene_id[overlapres$overlap > (1/3) & overlapres$overlap <= (2/3)]
  genesHighOverlapT <- overlapres$gene_id[overlapres$overlap > (2/3)]

  genesLowOverlap <- genesLowOverlapT[genesLowOverlapT %in% genestouse]
  genesMedOverlap <- genesMedOverlapT[genesMedOverlapT %in% genestouse]
  genesHighOverlap <- genesHighOverlapT[genesHighOverlapT %in% genestouse]


  #Names of genes done by splitting up into three equal groups by overlap quantile
  quan <- quantile(overlapres$overlap, probs = c((1/3), (2/3), 1))

  genesLowerThirdOverlapT <- overlapres$gene_id[overlapres$overlap <= quan[1]]
  genesMidThirdOverlapT <- overlapres$gene_id[overlapres$overlap > quan[1] & overlapres$overlap <= quan[2]]
  genesHighThirdOverlapT <- overlapres$gene_id[overlapres$overlap > quan[2]]

  genesLowerThirdOverlap <- genesLowerThirdOverlapT[genesLowerThirdOverlapT %in% genestouse]
  genesMidThirdOverlap <- genesMidThirdOverlapT[genesMidThirdOverlapT %in% genestouse]
  genesHighThirdOverlap <- genesHighThirdOverlapT[genesHighThirdOverlapT %in% genestouse]

  #Names of genes splitting the sum of the total gene expression into three groups by quantile
  TGEbyGenetemp <- aggregate(abGene$SumTGE, by = list(abGene$gene_id), FUN = "mean", na.rm = TRUE)
  colnames(TGEbyGenetemp) <- c("gene_id", "SumTGE")
  TGEbyGene <- subset(TGEbyGenetemp, TGEbyGenetemp$gene_id %in% genestouse)

  quanTGE <- quantile(TGEbyGene$SumTGE, probs = c((1/3), (2/3), 1))
  genesLowTGE <- TGEbyGene$gene_id[TGEbyGene$SumTGE < quanTGE[1]]
  genesMediumTGE <- TGEbyGene$gene_id[TGEbyGene$SumTGE >= quanTGE[1] & TGEbyGene$SumTGE < quanTGE[2]]
  genesHighTGE <- TGEbyGene$gene_id[TGEbyGene$SumTGE >= quanTGE[2]]

  genesSmallNTransLowerThirdOverlap <- intersect(genesSmallNTrans, genesLowerThirdOverlap)
  genesSmallNTransMidThirdOverlap <- intersect(genesSmallNTrans, genesMidThirdOverlap)
  genesSmallNTransHighThirdOverlap <- intersect(genesSmallNTrans, genesHighThirdOverlap)

  genesSmallNTransLowTGE <- intersect(genesSmallNTrans, genesLowTGE)
  genesSmallNTransMediumTGE <- intersect(genesSmallNTrans, genesMediumTGE)
  genesSmallNTransHighTGE <- intersect(genesSmallNTrans, genesHighTGE)


  genesMediumNTransLowerThirdOverlap <- intersect(genesMediumNTrans, genesLowerThirdOverlap)
  genesMediumNTransMidThirdOverlap <- intersect(genesMediumNTrans, genesMidThirdOverlap)
  genesMediumNTransHighThirdOverlap <- intersect(genesMediumNTrans, genesHighThirdOverlap)

  genesLargeNTransLowerThirdOverlap <- intersect(genesLargeNTrans, genesLowerThirdOverlap)
  genesLargeNTransMidThirdOverlap <- intersect(genesLargeNTrans, genesMidThirdOverlap)
  genesLargeNTransHighThirdOverlap <- intersect(genesLargeNTrans, genesHighThirdOverlap)


  genes2ColumnsHighTGE <- intersect(genes2Columns, genesHighTGE)
  genesWithNoOtherGroupHighTGE <-  intersect(genesWithNoOtherGroup, genesHighTGE)
  genesWithNoOtherGroupAnd2ColumnsHighTGE <- intersect(genesWithNoOtherGroupAnd2Columns, genesHighTGE)

  genesWithNoOtherGroupSmallNTrans <- intersect(genesWithNoOtherGroup, genesSmallNTrans)
  genesWithNoOtherGroupMediumNTrans <- intersect(genesWithNoOtherGroup, genesMediumNTrans)
  genesWithNoOtherGroupLargeNTrans <- intersect(genesWithNoOtherGroup, genesLargeNTrans)

  genes2ColumnsHighOverlap <- intersect(genes2Columns, genesHighOverlap)
  genesWithNoOtherGroupHighOverlap <-  intersect(genesWithNoOtherGroup, genesHighOverlap)
  genesWithNoOtherGroupAnd2ColumnsHighOverlap <- intersect(genesWithNoOtherGroupAnd2Columns, genesHighOverlap)

  genes2ColumnsHighThirdOverlap <- intersect(genes2Columns, genesHighThirdOverlap)
  genesWithNoOtherGroupHighThirdOverlap <-  intersect(genesWithNoOtherGroup, genesHighThirdOverlap)
  genesWithNoOtherGroupAnd2ColumnsHighThirdOverlap <- intersect(genesWithNoOtherGroupAnd2Columns, genesHighThirdOverlap)


  
  if(useInfRV==TRUE){
    if(DRIMSeqFiltering==TRUE & infReps=="Boot"){
      infRVDRIMSeqFilteringBoot <- TRUE
    }
    if(infRVDRIMSeqFilteringBoot==TRUE){
      InfRVRes <- loadRData("InfRVResAllForBootSamps.RData")
    }else{
      InfRVRes <- loadRData("InfRVResAllForGibbsSamps.RData")
    }
    
    quanMeanGeneMaxInfRVTPM <- quantile(InfRVRes$MeanGeneMaxInfRVTPM, probs = c((1/3), (2/3), 1))
    quanMeanGeneMedianInfRVTPM <- quantile(InfRVRes$MeanGeneMedianInfRVTPM, probs = c((1/3), (2/3), 1))
    quanMeanGeneMeanInfRVTPM <- quantile(InfRVRes$MeanGeneMeanInfRVTPM, probs = c((1/3), (2/3), 1))
    
    genesMeanGeneMaxInfRVTPMLow <- InfRVRes$gene_id[InfRVRes$MeanGeneMaxInfRVTPM <= quanMeanGeneMaxInfRVTPM[1]]
    genesMeanGeneMaxInfRVTPMMedium <- InfRVRes$gene_id[InfRVRes$MeanGeneMaxInfRV > quanMeanGeneMaxInfRVTPM[1] & InfRVRes$MeanGeneMaxInfRVTPM <= quanMeanGeneMaxInfRVTPM[2]]
    genesMeanGeneMaxInfRVTPMHigh <- InfRVRes$gene_id[InfRVRes$MeanGeneMaxInfRVTPM > quanMeanGeneMaxInfRVTPM[2]]
    
    genesMeanGeneMedianInfRVTPMLow <- InfRVRes$gene_id[InfRVRes$MeanGeneMedianInfRVTPM <= quanMeanGeneMedianInfRVTPM[1]]
    genesMeanGeneMedianInfRVTPMMedium <- InfRVRes$gene_id[InfRVRes$MeanGeneMedianInfRVTPM > quanMeanGeneMedianInfRVTPM[1] & InfRVRes$MeanGeneMedianInfRVTPM <= quanMeanGeneMedianInfRVTPM[2]]
    genesMeanGeneMedianInfRVTPMHigh <- InfRVRes$gene_id[InfRVRes$MeanGeneMedianInfRVTPM> quanMeanGeneMedianInfRVTPM[2]]
    
    genesMeanGeneMeanInfRVTPMLow <- InfRVRes$gene_id[InfRVRes$MeanGeneMeanInfRVTPM <= quanMeanGeneMeanInfRVTPM[1]]
    genesMeanGeneMeanInfRVTPMMedium <- InfRVRes$gene_id[InfRVRes$MeanGeneMeanInfRVTPM > quanMeanGeneMeanInfRVTPM[1] & InfRVRes$MeanGeneMeanInfRVTPM <= quanMeanGeneMeanInfRVTPM[2]]
    genesMeanGeneMeanInfRVTPMHigh <- InfRVRes$gene_id[InfRVRes$MeanGeneMeanInfRVTPM > quanMeanGeneMeanInfRVTPM[2]]
    
    
    quanMeanGeneMaxInfRVTPM2 <- quantile(InfRVRes$MeanGeneMaxInfRVTPM, probs = c(0.80, 0.90))
    quanMeanGeneMedianInfRVTPM2 <- quantile(InfRVRes$MeanGeneMedianInfRVTPM, probs = c(0.80, 0.90))
    quanMeanGeneMeanInfRVTPM2 <- quantile(InfRVRes$MeanGeneMeanInfRVTPM, probs = c(0.80, 0.90))
    
    
    genesMeanGeneMaxInfRVTPMEighty <- InfRVRes$gene_id[InfRVRes$MeanGeneMaxInfRVTPM > quanMeanGeneMaxInfRVTPM2[1] & InfRVRes$MeanGeneMaxInfRVTPM <= quanMeanGeneMaxInfRVTPM2[2]]
    genesMeanGeneMedianInfRVTPMEighty <- InfRVRes$gene_id[InfRVRes$MeanGeneMedianInfRVTPM > quanMeanGeneMedianInfRVTPM2[1] & InfRVRes$MeanGeneMedianInfRVTPM <= quanMeanGeneMedianInfRVTPM2[2]]
    genesMeanGeneMeanInfRVTPMEighty <- InfRVRes$gene_id[InfRVRes$MeanGeneMeanInfRVTPM > quanMeanGeneMeanInfRVTPM2[1] & InfRVRes$MeanGeneMeanInfRVTPM <= quanMeanGeneMeanInfRVTPM2[2]]
    
    genesMeanGeneMaxInfRVTPMGrEighty <- InfRVRes$gene_id[InfRVRes$MeanGeneMaxInfRVTPM > quanMeanGeneMaxInfRVTPM2[1]]
    genesMeanGeneMedianInfRVTPMGrEighty <- InfRVRes$gene_id[InfRVRes$MeanGeneMedianInfRVTPM > quanMeanGeneMedianInfRVTPM2[1]]
    genesMeanGeneMeanInfRVTPMGrEighty <- InfRVRes$gene_id[InfRVRes$MeanGeneMeanInfRVTPM > quanMeanGeneMeanInfRVTPM2[1]]
    
    genesMeanGeneMaxInfRVTPMGrNinety <- InfRVRes$gene_id[InfRVRes$MeanGeneMaxInfRVTPM > quanMeanGeneMaxInfRVTPM2[2]]
    genesMeanGeneMedianInfRVTPMGrNinety <- InfRVRes$gene_id[InfRVRes$MeanGeneMedianInfRVTPM > quanMeanGeneMedianInfRVTPM2[2]]
    genesMeanGeneMeanInfRVTPMGrNinety <- InfRVRes$gene_id[InfRVRes$MeanGeneMeanInfRVTPM > quanMeanGeneMeanInfRVTPM2[2]]
    
    quanMeanGeneMaxilrScale <- quantile(InfRVRes$MeanGeneMaxVarInfRepilrScale, probs = c((1/3), (2/3), 0.80, 0.90))
    quanMeanGeneMedianilrScale <- quantile(InfRVRes$MeanGeneMedianVarInfRepilrScale, probs = c((1/3), (2/3), 0.80, 0.90))
    quanMeanGeneMeanilrScale <- quantile(InfRVRes$MeanGeneMeanVarInfRepilrScale, probs = c((1/3), (2/3), 0.80, 0.90))
    
    
    genesMeanGeneMaxVarilrScaleLow <- InfRVRes$gene_id[InfRVRes$MeanGeneMaxVarInfRepilrScale <= quanMeanGeneMaxilrScale[1]]
    genesMeanGeneMaxVarilrScaleMedium <- InfRVRes$gene_id[InfRVRes$MeanGeneMaxVarInfRepilrScale > quanMeanGeneMaxilrScale[1] & InfRVRes$MeanGeneMaxVarInfRepilrScale <= quanMeanGeneMaxilrScale[2]]
    genesMeanGeneMaxVarilrScaleHigh <- InfRVRes$gene_id[InfRVRes$MeanGeneMaxVarInfRepilrScale > quanMeanGeneMaxilrScale[2]]
    
    genesMeanGeneMaxVarilrScaleGr80 <- InfRVRes$gene_id[InfRVRes$MeanGeneMaxVarInfRepilrScale > quanMeanGeneMaxilrScale[3]]
    genesMeanGeneMaxVarilrScaleGr90 <- InfRVRes$gene_id[InfRVRes$MeanGeneMaxVarInfRepilrScale > quanMeanGeneMaxilrScale[4]]
    
    
    
    
    genesMeanGeneMeanVarilrScaleLow <- InfRVRes$gene_id[InfRVRes$MeanGeneMeanVarInfRepilrScale <= quanMeanGeneMeanilrScale[1]]
    genesMeanGeneMeanVarilrScaleMedium <- InfRVRes$gene_id[InfRVRes$MeanGeneMeanVarInfRepilrScale > quanMeanGeneMeanilrScale[1] & InfRVRes$MeanGeneMeanVarInfRepilrScale <= quanMeanGeneMeanilrScale[2]]
    genesMeanGeneMeanVarilrScaleHigh <- InfRVRes$gene_id[InfRVRes$MeanGeneMeanVarInfRepilrScale > quanMeanGeneMeanilrScale[2]]
    
    genesMeanGeneMeanVarilrScaleGr80 <- InfRVRes$gene_id[InfRVRes$MeanGeneMeanVarInfRepilrScale > quanMeanGeneMeanilrScale[3]]
    genesMeanGeneMeanVarilrScaleGr90 <- InfRVRes$gene_id[InfRVRes$MeanGeneMeanVarInfRepilrScale > quanMeanGeneMeanilrScale[4]]
    
    
    
    
    genesMeanGeneMedianVarilrScaleLow <- InfRVRes$gene_id[InfRVRes$MeanGeneMedianVarInfRepilrScale <= quanMeanGeneMedianilrScale[1]]
    genesMeanGeneMedianVarilrScaleMedium <- InfRVRes$gene_id[InfRVRes$MeanGeneMedianVarInfRepilrScale > quanMeanGeneMedianilrScale[1] & InfRVRes$MeanGeneMedianVarInfRepilrScale <= quanMeanGeneMedianilrScale[2]]
    genesMeanGeneMedianVarilrScaleHigh <- InfRVRes$gene_id[InfRVRes$MeanGeneMedianVarInfRepilrScale > quanMeanGeneMedianilrScale[2]]
    
    genesMeanGeneMedianVarilrScaleGr80 <- InfRVRes$gene_id[InfRVRes$MeanGeneMedianVarInfRepilrScale > quanMeanGeneMedianilrScale[3]]
    genesMeanGeneMedianVarilrScaleGr90 <- InfRVRes$gene_id[InfRVRes$MeanGeneMedianVarInfRepilrScale > quanMeanGeneMedianilrScale[4]]
    
  }

  curr_env <- environment()
  print(curr_env)

  environment(SensFprRes) <- curr_env
  #restemp <- lapply(methds, SensFprRes)
  #names(restemp) <- methds
  #browser()
  restemp <- SensFprRes(curr_methd)
  #names(restemp) <- curr_methd


  # for(j in 1:length(restemp)){
  #   curr_methd <- names(restemp)[j]
  #
  # }
  #AllPvals <- restemp$AllPvals


  assign(paste0("ResListChange", curr_change, curr_methd), restemp)
  #list2env(restemp[[1]], envir = .GlobalEnv)
  #list2env(restemp, envir = .GlobalEnv)


  #rm calcSensFpr function just because it matches the pattern expression below and there's no need to save it again
  rm(calcSensFpr)
  rm(SensFprRes)
  
  fil <- c(ls(pattern="ResListChange"), ls(pattern="genes"))

  if(TwentySamplesOnlyCase==TRUE){
    dir_modif <- "TwentySamples"
  }else{
    dir_modif <- ""
  }
  
  if(infReps=="Gibbs" | infReps=="GibbsThin100"){
    if(DRIMSeqFiltering==TRUE){
      dir_piece_3 <- "PowerResDRIMSeqFiltering/"
    }else if(DRIMSeqFiltering==FALSE){
      dir_piece_3 <- "PowerRes/"
    }
  }else if(infReps=="GibbsThin16"){
    if(DRIMSeqFiltering==TRUE){
      dir_piece_3 <- "PowerResThin16DRIMSeqFiltering/"
    }else if(DRIMSeqFiltering==FALSE){
      dir_piece_3 <- "PowerResThin16/"
    }
  }else if(infReps=="Boot"){
    if(DRIMSeqFiltering==TRUE){
      dir_piece_3 <- "PowerResBootDRIMSeqFiltering/"
    }else if(DRIMSeqFiltering==FALSE){
      dir_piece_3 <- "PowerResBoot/"
    }

  }
    if(FilterGenes==TRUE){
      direc <- paste0(def_wd, data_name, dir_piece_3, "AllResultsFilter", FilterCutoff, dir_modif, "/")
      if(!dir.exists(direc)) {dir.create(direc, recursive = TRUE)}
      save(list = fil,file = paste0(direc, data_name, "PowerResAllResultsFilter", FilterCutoff, "Change", curr_change, curr_methd, ".RData"), version = 2)
    }else{
      direc <- paste0(def_wd, data_name, dir_piece_3, "AllResults", dir_modif, "/")
      if(!dir.exists(direc)) {dir.create(direc, recursive = TRUE)}
      save(list = fil,file = paste0(direc, data_name, "PowerResAllResultsChange", curr_change, curr_methd, ".RData"), version = 2)
    }
  
}







SaveFullinfRepDat <- function(CreateabAndcntDatasets = TRUE){

  if(GibbsSamps==TRUE){
    load_dir <- paste0(def_wd1, "GibbsSamps")
    infReps <- "Gibbs"
  }else{
    load_dir <- paste0(def_wd1, "BootSamps")
    infReps <- "Boot"
  }


  #Now, need to generate dataframe that has the following columns: tx_id, Sample1TPM, Sample2TPM, ..., GibbsRep#
  #To do this, need to load in the GibbsSamps dataframe for each sample, select only a current subset of genes
  #(since using all of them at once would take up too much memory) and merge them all together
  #and then run generateData as done below

  genestouse <- filteredgenenames

  #Function to split character vector into n (approximately) equally sized chunks
  #Taken from https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
  #Originally by mathheadinclouds and edited by Dis Shishkov
  chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
  genestouse_split <- chunk2(genestouse, nparts)
  all_genes_annotation_split <- chunk2(unique(sort(cntGene$gene_id)), nparts)
  
  print(paste0("Number of rows in cntGene to use to save FullInfRepDat is ", nrow(cntGene)))
  print(paste0("Total number of genes for the entire annotation is ", length(unique(sort(cntGene$gene_id)))))

  #Gene list only for those genes that will be used in the compositional analysis
  curr_genes <- genestouse_split[[curr_part_num]]

  #Gene list for all genes, even ones that only have 1 trans, etc
  #Needed for the compositionalMI Power analysis to calculate updated TPMs, etc
  curr_genes_all_genes_annotation <- all_genes_annotation_split[[curr_part_num]]


  infRepFiles <- mixedsort(list.files(load_dir, pattern = ".RData", full.names = TRUE))

  #Load the first sample to get the ngibbs/nboot value and assign that value fo ninfreps
  d <- loadRData(infRepFiles[1], objNameToGet = paste0(infReps, "SampsSample1"))
  if(GibbsSamps==TRUE){
    ninfreps <- attr(d, "ngibbs")
  }else{
    ninfreps <- attr(d, "nboot")
  }


  lengths <- data.frame(QuantSalmon$length)
  lengths$tx_id <- rownames(lengths)

  InfRepsKey <- data.frame(rep(key$Sample, ninfreps), rep(key$Condition, ninfreps), rep(key$Identifier, ninfreps))
  colnames(InfRepsKey) <- c("Sample", "Condition", "Identifier")



  #x is which of the samples from infRepFiles is being used
  f1 <- function(x, genes, infRepFiles){
    print(paste0("Currently loading file ", x))
    curr_sampt1 <- strsplit(x, "Sample")[[1]][2]
    curr_sampt2 <- strsplit(curr_sampt1, ".RData")[[1]][1]
    print(paste0("Current Sample then is ", curr_sampt2))
    d <- loadRData(x, objNameToGet = paste0(infReps, "SampsSample", curr_sampt2))
    #Apply over each of the genestouse_split subsets
    d2 <- subset(d, d$gene_id %in% genes)
    #d2 <- lapply(genestouse_split, function(x, d){return(subset(d, d$gene_id %in% x))}, d=d)
    return(d2)
  }

  #Combine the data for the current subset of genes to all biological samples
  #Note that this subset comes from all genes from cntGene, not just the genes that will be used for the Comp analysis
  #This could be useful if you ever wanted a dataset that combined inferential replicates for genes that even have
  #only one transcript and for the power analysis, where updated TPMs need to be calculated based on all genes
    #not just the ones that pass filtering
  AllDat <- lapply(infRepFiles, f1, genes = curr_genes_all_genes_annotation, infRepFiles = infRepFiles)

  #Create data file, and call it datt
  datt <- AllDat[[1]]

  #i here indexes biological samples, not the different parts the array_val indexes
  for (i in 2:length(AllDat)){
    print(paste0("Currently processing biological sample ", i, " out of ", length(AllDat)))
    curr <- AllDat[[i]]

    #Check and make sure the files are in the same order such that the trans/Gibbs rep combos match up properly
    if(sum(paste0(datt$tx_id, "GibbsRep", datt$GibbsRep)!=paste0(curr$tx_id, "GibbsRep", curr$GibbsRep))!=0){
      stop("Dataframes are not in the right order and the code will not work properly")
    }
    col <- colnames(curr)[grep("TPM", colnames(curr))]
    datt$newcol <-  curr[,col, with = F]
    colnames(datt)[colnames(datt)=="newcol"] <- col

    col2 <- colnames(curr)[grep("Cnt", colnames(curr))]
    datt$newcol2 <-  curr[,col2, with = F]
    colnames(datt)[colnames(datt)=="newcol2"] <- col2

    rm(curr)
    gc()
    print(gc())

  }

  #datt2 <- datt[order(datt$gene_id, datt$tx_id, datt$GibbsRep),]
  if(!dir.exists(paste0(save_dir, sd3))){
    dir.create(paste0(save_dir, sd3))
  }

  assign(paste0("CompTimedattPart", curr_part_num), proc.time() - startTime)
  nam <- c("datt", paste0("CompTimedattPart", curr_part_num))
  print(paste0("Directory where files will be saved is ", paste0(save_dir, sd3)))
  save(list = nam, file = paste0(save_dir, sd3, "FullinfRepDatPart", curr_part_num, ".RData"))


  rm(datt)
  rm(AllDat)
  gc()
  print(gc())
  
  if(CreateabAndcntDatasets==FALSE){
    print("Saving of the FullinfRep Datasets is complete and the file stopped because CreateabAndcntDatasets=FALSE.  IF you want these files to be created change this argument to TRUE")
    return(NULL)
  }


  #stop("Code below generated abDatasets/cntDatasets for all Gibbs replicates, only run above this line for now")

  #Combine the data for the current subset of genes to all biological samples
  #Note that this subset is only genes that are used in abDatasets after filtering, unlike the same lines of code above which are 
    #for all genes in the annotation
  #ie the results below are on the genes that pass filtering, while the ones above are on all genes
  AllDat2 <- lapply(infRepFiles, f1, genes = curr_genes, infRepFiles = infRepFiles)

  #Create data file, datt
  datt2 <- AllDat2[[1]]

  #i here indexes biological samples, not the different parts the array_val indexes
  for (i in 2:length(AllDat2)){
    print(paste0("Currently processing biological sample ", i, " out of ", length(AllDat2)))
    curr <- AllDat2[[i]]

    #Check and make sure the files are in the same order such that the trans/Gibbs rep combos match up properly
    if(sum(paste0(datt2$tx_id, "GibbsRep", datt2$GibbsRep)!=paste0(curr$tx_id, "GibbsRep", curr$GibbsRep))!=0){
      stop("Dataframes are not in the right order and the code will not work properly")
    }
    col <- colnames(curr)[grep("TPM", colnames(curr))]
    datt2$newcol <-  curr[,col, with = F]
    colnames(datt2)[colnames(datt2)=="newcol"] <- col

    col2 <- colnames(curr)[grep("Cnt", colnames(curr))]
    datt2$newcol2 <-  curr[,col2, with = F]
    colnames(datt2)[colnames(datt2)=="newcol2"] <- col2

    rm(curr)
    gc()
    print(gc())

  }

  if(DRIMSeqFiltering==TRUE & CreateabAndcntDatasets==TRUE){

    abGNO <- lapply(curr_genes, generateData, dat = datt2, nsamp = nrow(key), abundance = TRUE, abCompDatasets = abDatasetsFiltered, useOtherGroups = FALSE, useExistingOtherGroups = FALSE, infReps = infReps, ninfreps = ninfreps)
    names(abGNO) <- curr_genes

    assign(paste0("abDatasets", type, "NoOtherGroupsFilteredPart", curr_part_num), abGNO)
    obj2 <- c(paste0("abDatasets", type, "NoOtherGroupsFilteredPart", curr_part_num))

    rm(abGNO)
    gc()

    if(!dir.exists(paste0(save_dir, sd1))){
      dir.create(paste0(save_dir, sd1))
    }

    save(InfRepsKey, ninfreps, nsamp, list = obj2, file = paste0(save_dir, sd1, "abDatasets", type, "NoOtherGroupsFilteredPart", curr_part_num, ".RData"))


    rm(list = paste0("abDatasets", type, "NoOtherGroupsFilteredPart", curr_part_num))
    gc()


    cntGNO <- lapply(curr_genes, generateData, dat = datt2, nsamp = nrow(key), abundance = FALSE, abCompDatasets = abDatasetsFiltered, useOtherGroups = FALSE, useExistingOtherGroups = FALSE, infReps = infReps, ninfreps = ninfreps)
    names(cntGNO) <- curr_genes

    assign(paste0("cntDatasets", type, "NoOtherGroupsFilteredPart", curr_part_num), cntGNO)
    obj2 <- c(paste0("cntDatasets", type, "NoOtherGroupsFilteredPart", curr_part_num))

    rm(cntGNO)
    gc()

    if(!dir.exists(paste0(save_dir, sd2))){
      dir.create(paste0(save_dir, sd2))
    }

    #save(abDatasetsGibbsRepsNoOtherGroups, InfRepsKey, ninfreps, nsamp, file = "abDatasetsGibbsRepsNoOtherGroups.RData")
    save(InfRepsKey, ninfreps, nsamp, list = obj2, file = paste0(save_dir, sd2, "cntDatasets", type,  "NoOtherGroupsFilteredPart", curr_part_num, ".RData"))
  }

  if(DRIMSeqFiltering==FALSE & CreateabAndcntDatasets==TRUE){
    abG <- lapply(curr_genes, generateData, dat = datt2, nsamp = nrow(key), abundance = TRUE, abCompDatasets = abDatasets, useExistingOtherGroups = TRUE, infReps = infReps, ninfreps = ninfreps)
    # abG <- laply(curr_genes[100:200], generateData, dat = datt, nsamp = nrow(key), abundance = TRUE, abCompDatasets = abDatasets, useExistingOtherGroups = TRUE, infReps = infReps, ninfreps = ninfreps, .inform = TRUE, .progress = "text")
    names(abG) <- curr_genes

    assign(paste0("abDatasets", type, "Part", curr_part_num), abG)

    rm(abG)
    gc()
    obj1 <- c(paste0("abDatasets", type, "Part", curr_part_num))
    #save(abDatasetsGibbsReps, InfRepsKey, ninfreps, nsamp, file = "abDatasetsGibbsReps.RData")

    if(!dir.exists(paste0(save_dir, sd1))){
      dir.create(paste0(save_dir, sd1))
    }

    save(InfRepsKey, ninfreps, nsamp, list = obj1, file = paste0(save_dir, sd1, "abDatasets", type, "Part", curr_part_num, ".RData"))


    rm(list = paste0("abDatasets", type, "Part", curr_part_num))
    gc()

    cntG <- lapply(curr_genes, generateData, dat = datt2, nsamp = nrow(key), abundance = FALSE, abCompDatasets = abDatasets, useExistingOtherGroups = TRUE, infReps = infReps, ninfreps = ninfreps)
    names(cntG) <- curr_genes

    assign(paste0("cntDatasets", type, "Part", curr_part_num), cntG)
    obj1 <- c(paste0("cntDatasets", type, "Part", curr_part_num))

    rm(cntG)
    gc()

    if(!dir.exists(paste0(save_dir, sd2))){
      dir.create(paste0(save_dir, sd2))
    }

    save(InfRepsKey, ninfreps, nsamp, list = obj1, file = paste0(save_dir, sd2, "cntDatasets", type, "Part", curr_part_num, ".RData"))

    rm(list = paste0("cntDatasets", type, "Part", curr_part_num))
    gc()


    abGNO <- lapply(curr_genes, generateData, dat = datt2, nsamp = nrow(key), abundance = TRUE, abCompDatasets = abDatasetsNoOtherGroups, useOtherGroups = FALSE, useExistingOtherGroups = FALSE, infReps = infReps, ninfreps = ninfreps)
    names(abGNO) <- curr_genes

    assign(paste0("abDatasets", type, "NoOtherGroupsPart", curr_part_num), abGNO)
    obj2 <- c(paste0("abDatasets", type, "NoOtherGroupsPart", curr_part_num))

    rm(abGNO)
    gc()

    #save(abDatasetsGibbsRepsNoOtherGroups, InfRepsKey, ninfreps, nsamp, file = "abDatasetsGibbsRepsNoOtherGroups.RData")
    save(InfRepsKey, ninfreps, nsamp, list = obj2, file = paste0(save_dir, sd1, "abDatasets", type, "NoOtherGroupsPart", curr_part_num, ".RData"))


    rm(list = paste0("abDatasets", type, "NoOtherGroupsPart", curr_part_num))
    gc()


    cntGNO <- lapply(curr_genes, generateData, dat = datt2, nsamp = nrow(key), abundance = FALSE, abCompDatasets = abDatasetsNoOtherGroups, useOtherGroups = FALSE, useExistingOtherGroups = FALSE, infReps = infReps, ninfreps = ninfreps)
    names(cntGNO) <- curr_genes

    assign(paste0("cntDatasets", type, "NoOtherGroupsPart", curr_part_num), cntGNO)
    obj2 <- c(paste0("cntDatasets", type, "NoOtherGroupsPart", curr_part_num))

    rm(cntGNO)
    gc()

    #save(abDatasetsGibbsRepsNoOtherGroups, InfRepsKey, ninfreps, nsamp, file = "abDatasetsGibbsRepsNoOtherGroups.RData")
    save(InfRepsKey, ninfreps, nsamp, list = obj2, file = paste0(save_dir, sd2, "cntDatasets", type,  "NoOtherGroupsPart", curr_part_num, ".RData"))
  }


}




SaveInfRepSpecificPartFullinfRepDatFiles <- function(curr_part_num, load_dir, nparts, ninfreps){
    #print(paste0("Current Part Number is ", curr_part_num, " out of ", nparts))
    curr <- loadRData(paste0(load_dir, "FullinfRepDatPart", curr_part_num, ".RData"), objNameToGet = "datt")
    for(j in 1:ninfreps){
      infRepPartData <- subset(curr, curr$infRepNum==j)
      if(!dir.exists(paste0(load_dir, "FullinfRepDatPart", curr_part_num, "/"))){
        dir.create(paste0(load_dir, "FullinfRepDatPart", curr_part_num, "/"))
      }
      save(infRepPartData, file = paste0(load_dir, "FullinfRepDatPart", curr_part_num, "/", "FullinfRepDatPart", curr_part_num, "InfRep", j, ".RData"))
    }
    rm(curr)
    gc()
  
}


SaveWithinSubjCovMatrices <- function(){
  if(DRIMSeqFiltering==FALSE){
    load(paste0(save_dir, abDatasetsSubDir, "abDatasets", dirpiece, "Part", curr_part, ".RData"))
    load(paste0(save_dir, abDatasetsSubDir, "abDatasets", dirpiece, "NoOtherGroupsPart", curr_part, ".RData"))

    abDatasetsToUse1 <- get(paste0("abDatasets", dirpiece, "Part", curr_part))
    genestouse1 <- names(abDatasetsToUse1)

    abDatasetsToUse2 <- get(paste0("abDatasets", dirpiece, "NoOtherGroupsPart", curr_part))
    genestouse2 <- names(abDatasetsToUse2)
  }else{
    load(paste0(save_dir, abDatasetsSubDir, "abDatasets", dirpiece, "NoOtherGroupsFilteredPart", curr_part, ".RData"))

    abDatasetsToUse3 <- get(paste0("abDatasets", dirpiece, "NoOtherGroupsFilteredPart", curr_part))
    genestouse3 <- names(abDatasetsToUse3)
  }





  if(DRIMSeqFiltering==FALSE){
    #ilrMeansCovsOtherGroups <- laply(genestouse, calcIlrMeansCovs, dat = abDatasetsToUse1, nsamp = nsamp, .progress = "text", .drop = F)
    ilrMeansCovsOtherGroups <- lapply(genestouse1, calcIlrMeansCovs, dat = abDatasetsToUse1, nsamp = nsamp)
    names(ilrMeansCovsOtherGroups) <- genestouse1

    #These lines are only used to test the covariance matricies of the Gibbs/Bootstrap samples for issues, usually not needed
    # ilrMeansCovsOtherGroups <- lapply(genestouse1[1:10], calcIlrMeansCovs, dat = abDatasetsToUse1, nsamp = nsamp, investigateGibbsSamps = TRUE, invesSamp = 3)
    #names(ilrMeansCovsOtherGroups) <- genestouse1[1:10]

    assign(paste0("ilrMeansCovsOtherGroupsPart", curr_part), ilrMeansCovsOtherGroups)

    if(!dir.exists(paste0(save_dir, ilrMeansCovsSubDir))){
      dir.create(paste0(save_dir, ilrMeansCovsSubDir))
    }

    save(list = paste0("ilrMeansCovsOtherGroupsPart", curr_part), file = paste0(save_dir, ilrMeansCovsSubDir,"ilrMeansCovsOtherGroupsPart", curr_part, ".RData"))


    # ilrMeansCovsNoOtherGroups <- laply(genestouse2, calcIlrMeansCovs, dat = abDatasetsToUse2, nsamp = nsamp, .progress = "text", .drop = F)
    ilrMeansCovsNoOtherGroups <- lapply(genestouse2, calcIlrMeansCovs, dat = abDatasetsToUse2, nsamp = nsamp)
    names(ilrMeansCovsNoOtherGroups) <- genestouse2

    assign(paste0("ilrMeansCovsNoOtherGroupsPart", curr_part), ilrMeansCovsNoOtherGroups)
    save(list = paste0("ilrMeansCovsNoOtherGroupsPart", curr_part), file = paste0(save_dir, ilrMeansCovsSubDir,"ilrMeansCovsNoOtherGroupsPart", curr_part, ".RData"))
  }

  if(DRIMSeqFiltering==TRUE){
    load(paste0(dir1, "cntGenecntsScaledTPMFiltered.RData"))
    #genen <- names(abDatasetsToUse2)
    filterabDatasets <- function(x, cntGeneFiltered, abDatasetsToUse2){
      curr_gene <- x
      txs <- cntGeneFiltered$tx_id[cntGeneFiltered$gene_id==curr_gene]
      curr <- abDatasetsToUse2[[curr_gene]]
      curr2 <- curr[,colnames(curr) %in% txs, drop = FALSE]
      if(is.null(ncol(curr2))){
        return(NULL)
      }

      #if(ncol(curr2)!=length(txs)){stop("ree")}

      if(ncol(curr2)==0){
        curr2 <- NULL
      }
      return(curr2)
    }

    #abDatasetsToUse3 <- lapply(genen, filterabDatasets, cntGeneFiltered = cntGeneFiltered, abDatasetsToUse2 = abDatasetsToUse2)
    #abDatasetsToUse3 <- laply(genen, filterabDatasets, cntGeneFiltered = cntGeneFiltered, abDatasetsToUse2 = abDatasetsToUse2, .inform=T, .progress = "text")

    #names(abDatasetsToUse3) <- genen

    genestouse3 <- names(abDatasetsToUse3)

    d1 <- proc.time()
    ilrMeansCovsNoOtherGroupsFiltered <- lapply(genestouse3, calcIlrMeansCovs, dat = abDatasetsToUse3, nsamp = nsamp)
    proc.time() - d1
    #ilrMeansCovsNoOtherGroupsFiltered <- laply(genestouse3, calcIlrMeansCovs, dat = abDatasetsToUse3, nsamp = nsamp, .inform=T, .progress = T)
    names(ilrMeansCovsNoOtherGroupsFiltered) <- genestouse3

    assign(paste0("ilrMeansCovsNoOtherGroupsFilteredPart", curr_part), ilrMeansCovsNoOtherGroupsFiltered)

    if(!dir.exists(paste0(save_dir, ilrMeansCovsSubDir))){
      dir.create(paste0(save_dir, ilrMeansCovsSubDir))
    }
    save(list = paste0("ilrMeansCovsNoOtherGroupsFilteredPart", curr_part),
         file = paste0(save_dir, ilrMeansCovsSubDir,"ilrMeansCovsNoOtherGroupsFilteredPart", curr_part, ".RData"))
  }
}




updateilrMeansCovsAndabDatasets <- function(subset_data = FALSE, generateNewIlrMeansCovs = TRUE, InfRepSpecificPartFullinfRepDatFilesExist = FALSE,
                                            TwentySamplesTotalAnalysis = FALSE, ModifyAbundancesFromInfRepsForPowerAnalysis = FALSE, AbundancesFromInfRepsDir = NULL){
  print(paste0("Full infRep data is being loaded from ", load_dir))
  startTime <- proc.time()

  if(subset_data==TRUE){
    samps <- sub_key$Identifier
    names(y) <- samps
    key_to_use <- sub_key
  }else{
    samps <- key$Identifier
    names(y) <- samps
    key_to_use <- key
  }
  
  if(TwentySamplesTotalAnalysis==TRUE & (nrow(key_to_use)!=20 | length(samps)!=20)){
    stop("Specification of the input arguments is incorrect with TwentySamplesTotalAnalysis=TRUE")
  }


  #If the files already exist, skip it and move on

  if(TwentySamplesTotalAnalysis==TRUE){
    direc_modifier <- "TwentySamples"
  }else{
    direc_modifier <- ""
  }

  if(DRIMSeqFiltering==TRUE){
    if(infReps=="GibbsThin100"){
      curr_direc <- paste0(wd2, "UpdatedPowerDataDRIMSeqFiltering", direc_modifier, "/Change", change, "/", "GroupCombo", GroupNum, "/")
    }else if(infReps=="Boot"){
      curr_direc <- paste0(wd2, "UpdatedPowerDataBootDRIMSeqFiltering", direc_modifier, "/Change", change, "/", "GroupCombo", GroupNum, "/")
    }else if(infReps=="GibbsThin16"){
      curr_direc <- paste0(wd2, "UpdatedPowerDataGibbsThin16DRIMSeqFiltering", direc_modifier, "/Change", change, "/", "GroupCombo", GroupNum, "/")
    }else if (infReps=="Gibbs"){
      stop("For the purpose of this file specify Gibbs with the thinning, ie GibbsThin100, GibbsThin16, etc")
    }

  }else if(DRIMSeqFiltering==FALSE){
    if(infReps=="GibbsThin100"){
      curr_direc <- paste0(wd2, "UpdatedPowerData", direc_modifier, "/Change", change, "/", "GroupCombo", GroupNum, "/")
    }else if(infReps=="Boot"){
      curr_direc <- paste0(wd2, "UpdatedPowerDataBoot", direc_modifier, "/Change", change, "/", "GroupCombo", GroupNum, "/")
    }else if(infReps=="GibbsThin16"){
      curr_direc <- paste0(wd2, "UpdatedPowerDataGibbsThin16", direc_modifier, "/Change", change, "/", "GroupCombo", GroupNum, "/")
    }else if (infReps=="Gibbs"){
      stop("For the purpose of this file specify Gibbs with the thinning, ie GibbsThin100, GibbsThin16, etc")
    }
  }

  print(paste0("curr_direc is ", curr_direc))
  filteredgenenames <- names(abDatasetsFiltered)
  
  #Function to split character vector into n (approximately) equally sized chunks
  #Taken from https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
  #Originally by mathheadinclouds and edited by Dis Shishkov
  chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))


  
  filteredgenenames_split <- chunk2(filteredgenenames, nparts)
  
  print(paste0("Number of rows in cntGene to use to save FullInfRepDat is ", nrow(cntGene)))
  print(paste0("Total number of genes for the entire annotation is ", length(unique(sort(cntGene$gene_id)))))
  
  if(ModifyAbundancesFromInfRepsForPowerAnalysis==TRUE){
    if(DRIMSeqFiltering==FALSE){
      stop("You are attempting to run results for abundance values calculated from Inf Reps for the case not using DRIMSeq Filtering.  As of now those results have not been run.  If that changes modify this statement and make sure to 
           verify directory structure is working as expected.")
    }
    direc_to_save_resAbFromInfReps <- paste0(wd2, "UpdatedAbDatasetsAbFromInfRepsDRIMSeqFiltering", direc_modifier, "/Change", change, "/", "GroupCombo", GroupNum, "/")
    if(!dir.exists(direc_to_save_resAbFromInfReps)){dir.create(direc_to_save_resAbFromInfReps, recursive = T)}
    
    print(paste0("Directory to save updated abDatasets for abundances calculated from InfReps is ", direc_to_save_resAbFromInfReps))
    print("Calculating New AbDatasets for Abundances Calculated from InfReps")
    ModAbFromInfRepsForPowerAnalysis <- function(abFromInfRepFuncT, GibbsSampsT, countsFromAbundance = "scaledTPM", AbundancesFromInfRepsDir,  direc_to_save_resAbFromInfReps){
      abFromInfRepFunc <- abFromInfRepFuncT
      GibbsSamps <- GibbsSampsT
      fil_mod <- returnFilMod(abFromInfRepFunc = abFromInfRepFunc, GibbsSamps = GibbsSamps)
      
      if(file.exists(paste0(direc_to_save_resAbFromInfReps, "UpdatedabDatasetsChange", change, "GroupCombo", GroupNum, fil_mod, ".RData"))){
        return(NULL)
      }
      
      abGeneCurr <- loadRData(paste0(AbundancesFromInfRepsDir, "abGene", fil_mod, ".RData"), objNameToGet = "abGene")
      
      #The cntGene object that is loaded here must be regular counts always and not countsScaledTPM because the counts 
        #have to be updated and the TPM values calculated using the total library size (ie total count) of the 
        #sample.  Using the scaled TPM counts would result in calculating the TPM values using the wrong library size
      cntGeneCurr <- loadRData(paste0(AbundancesFromInfRepsDir, "cntGene", fil_mod, ".RData"), objNameToGet = "cntGene")
      
      print(paste0("Number of rows in cntGene and abGene are ", nrow(cntGene), " and ", nrow(abGene)))
      NewDataAb <- updateData(x = change, cntGene = cntGeneCurr, abGene = abGeneCurr,
                              abDatasetsOrig = abDatasetsFiltered, key = key_to_use, Group = y,
                              seed = NULL, tx2gene = tx2gene, useOtherGroups = useOtherGroups,
                              genestouse = names(abDatasetsFiltered), samps = samps,
                              genestochange = genestochange, CompMI = TRUE)
      
      #Make sure to only use genes contained in the abDatasets inputted to the function (abDatasets) because
      #only there have ilrMeansCovs data for now
      
      assign(paste0("UpdatedabDatasetsChange", change, "GroupCombo", GroupNum, fil_mod), NewDataAb$abDatasets)
      print(paste0("Length of updated abDatasets is ", length(NewDataAb$abDatasets), " (Should match dimensions of filtered abDatasets)"))
      print(paste0("Length of filtered abDatasets is ", length(abDatasetsFiltered)))
      
      rm(NewDataAb)
      gc()
      
      if(!dir.exists(curr_direc)){dir.create(curr_direc, recursive = TRUE)}
      
      save(list = paste0("UpdatedabDatasetsChange", change, "GroupCombo", GroupNum, fil_mod), file = paste0(direc_to_save_resAbFromInfReps, "UpdatedabDatasetsChange", change, "GroupCombo", GroupNum, fil_mod, ".RData"))
      proc.time() - startTime
      print(gc())
    }
    
    print("Calculating New AbDatasets for Abundances Calculated from InfReps Number 1 of 4")
    ModAbFromInfRepsForPowerAnalysis(abFromInfRepFuncT = "median", GibbsSampsT = TRUE, AbundancesFromInfRepsDir = AbundancesFromInfRepsDir,
                                     direc_to_save_resAbFromInfReps = direc_to_save_resAbFromInfReps)
    print(gc())
    
    print("Calculating New AbDatasets for Abundances Calculated from InfReps Number 2 of 4")
    ModAbFromInfRepsForPowerAnalysis(abFromInfRepFuncT = "mean", GibbsSampsT = TRUE, AbundancesFromInfRepsDir = AbundancesFromInfRepsDir,
                                     direc_to_save_resAbFromInfReps = direc_to_save_resAbFromInfReps)
    print(gc())
    
    print("Calculating New AbDatasets for Abundances Calculated from InfReps Number 3 of 4")
    ModAbFromInfRepsForPowerAnalysis(abFromInfRepFuncT = "median", GibbsSampsT = FALSE, AbundancesFromInfRepsDir = AbundancesFromInfRepsDir,
                                     direc_to_save_resAbFromInfReps = direc_to_save_resAbFromInfReps)
    print(gc())
    
    print("Calculating New AbDatasets for Abundances Calculated from InfReps Number 4 of 4")
    ModAbFromInfRepsForPowerAnalysis(abFromInfRepFuncT = "mean", GibbsSampsT = FALSE, AbundancesFromInfRepsDir = AbundancesFromInfRepsDir,
                                     direc_to_save_resAbFromInfReps = direc_to_save_resAbFromInfReps)
    print(gc())
    
    print("With ModifyAbundancesFromInfRepsForPowerAnalysis=TRUE only generate the new files corresponding to abundances 
         calculated using the mean/median of infreps and not the files for the infReps themselves. To get these other files rerun with ModifyAbundancesFromInfRepsForPowerAnalysis=FALSE")
    return(NULL)
    
  }

  if(!file.exists(paste0(curr_direc, "UpdatedabDatasetsChange", change, "GroupCombo", GroupNum, ".RData"))){
    #cntGene <- cntGenePart1
    #Note, set CompMI to TRUE even though this is not necessarily CompMI because there can be a problem with certain transcripts existing in the
    #observed data that do not exist in the gibbs samples, which would cause the code to break if this option is not used
    
    print(paste0("Number of rows in cntGene and abGene are ", nrow(cntGene), " and ", nrow(abGene)))
    NewData <- updateData(x = change, cntGene = cntGene, abGene = abGene,
                          abDatasetsOrig = abDatasetsFiltered, key = key_to_use, Group = y,
                          seed = NULL, tx2gene = tx2gene, useOtherGroups = useOtherGroups,
                          genestouse = names(abDatasetsFiltered), samps = samps,
                          genestochange = genestochange, CompMI = TRUE)

    #Make sure to only use genes contained in the abDatasets inputted to the function (abDatasets) because
    #only there have ilrMeansCovs data for now
    
    assign(paste0("UpdatedabDatasetsChange", change, "GroupCombo", GroupNum), NewData$abDatasets)
    print(paste0("Length of updated abDatasets is ", length(NewData$abDatasets), " (Should match dimensions of filtered abDatasets)"))
    print(paste0("Length of filtered abDatasets is ", length(abDatasetsFiltered)))



    if(!dir.exists(curr_direc)){dir.create(curr_direc, recursive = TRUE)}

    save(list = paste0("UpdatedabDatasetsChange", change, "GroupCombo", GroupNum), file = paste0(curr_direc, "UpdatedabDatasetsChange", change, "GroupCombo", GroupNum, ".RData"))
    proc.time() - startTime
  }
  

  



  #Group <- y
  names(Group) <- samps
  


  ###################
  #Generate new ilrMeansCovs objects corresponding to the updated data as well to ensure the level of measurement
  #error in the ilrCovs of the Gibbs samples is the same as the level in the regular counts
  #ie also multiply the counts in the gibbs samples for the necessary transcripts by the change amount
  stT <- proc.time()

  MajorTrans <- unique(subset(abGene[,c("gene_id", "MajorTransName")]))
  rownames(MajorTrans) <- MajorTrans$gene_id

  lencols <- c(which(colnames(cntGene) %in% c("tx_id")), grep("Len", colnames(cntGene)))
  lenTemp <- cntGene[,lencols]
  print(paste0("nrow of lenTemp is ", nrow(lenTemp)))

  templist <- list()
  updated_data <- list()

  NewilrMeansCovs <- list()
  newAbDatasetsGibbs <- list()
  #nparts
  startT <- proc.time()
  print(paste0("ninfreps is ", ninfreps))
  max_num_run <- 1
  for(j in 1:ninfreps){
    curr_infRep <- j
    
    if(file.exists(paste0(curr_direc, "UpdatedDataChange", change, "GroupCombo", GroupNum, "InfRepNum", curr_infRep, ".RData"))){
      max_num_run <- j
    }
  }  
  
  #Load the max_num_run from above but re-run that largest one just in case it got cut off while it was being saved
  for(j in max_num_run:ninfreps){
    curr_infRep <- j

    print(paste0("infRep number is ", j))
    Sys.sleep(sample(1:10, 1))

    sdd <- proc.time()
    for(i in 1:nparts){

      if(i %% 25 ==0){
        print(paste0("Current part is ", i, " within infrep num ", j))
      }
      
      if(InfRepSpecificPartFullinfRepDatFilesExist==TRUE){
        templist[[i]] <- loadRData(paste0(load_dir, "FullinfRepDatPart", i, "/", "FullinfRepDatPart", i, "InfRep", j, ".RData"), objNameToGet = "infRepPartData")
        gc()
      }else{
        load(paste0(load_dir, "FullinfRepDatPart", i, ".RData"))
        curr_part_dat <- subset(datt, datt$infRepNum==curr_infRep)
        curr_part_dat
        #genes_part[[i]] <- sort(unique(curr_part_dat$gene_id))
        rm(datt)
        rm(curr_part_dat)
        gc()
      }


    }

    new_dat <- data.frame(rbindlist(templist))
    #commentmarker
    #browser()
    #print(paste0("nrow of new_dat is ", nrow(new_dat)))

    new_dat2 <- merge(new_dat, lenTemp, by = "tx_id")
    rm(new_dat)
    gc()
    new_dat2TPMcols <- grep("TPM", colnames(new_dat2))
    new_dat3 <- new_dat2[,-new_dat2TPMcols]
    
    #print(paste0("nrow of new_dat2 is ", nrow(new_dat2)))
    #print(paste0("nrow of new_dat3 is ", nrow(new_dat3)))

    rm(new_dat2)
    gc()

    rownames(new_dat3) <- new_dat3$tx_id



    nsamp <- nrow(key_to_use)
    if(is.null(samps)==TRUE){
      samps <- key_to_use$Identifier
    }

    new_dat4 <- merge(new_dat3, MajorTrans, by = "gene_id")
    
    #print(paste0("nrow of new_dat4 is ", nrow(new_dat4)))

    rm(new_dat3)
    gc()
    new_dat4$MajorTrans <- as.numeric(new_dat4$tx_id==new_dat4$MajorTransName)

    rws <- (new_dat4$MajorTrans==1 & new_dat4$gene_id %in% genestochange)

    #print(paste0("Current Group is ", Group))
    
    #Extract the sample names and only modify those that correspond to the first condition
    SampsCond1 <- subset(key_to_use$Identifier, Group==levels(key_to_use$Condition)[1])
    
    #print(paste0("Current Samples in Condition 1 include ", SampsCond1))
    cols <- paste0(SampsCond1, "Cnt")
    
    ColsCond1 <- which(colnames(new_dat4) %in% cols)

    new_dat5 <- new_dat4
    new_dat5[rws,ColsCond1] <- new_dat5[rws,ColsCond1] * change
    rownames(new_dat5) <- new_dat5$tx_id
    
    print(paste0("nrow of new_dat5 is ", nrow(new_dat5)))

    rm(new_dat4)
    gc()


    new_TPM <- cntsToTPM(cnts = new_dat5, nsamp = length(Group), samps = samps)

    new_dat6 <- merge(new_dat5, new_TPM, by = "tx_id")
    
    #print(paste0("nrow of new_dat6 is ", nrow(new_dat6)))

    rm(new_dat5)
    rm(new_TPM)
    gc()

    colstouse <- c("tx_id", "gene_id", "infRepNum", "NTrans", paste0(samps, "TPM"))
    new_dat7 <- new_dat6[,colstouse]
    rm(new_dat6)
    gc()
    
    print(paste0("nrow of new_dat7 is ", nrow(new_dat7)))

    assign(paste0("UpdatedDataChange", change, "GroupCombo", GroupNum, "InfRepNum", curr_infRep), new_dat7)

    save(list = paste0("UpdatedDataChange", change, "GroupCombo", GroupNum, "InfRepNum", curr_infRep), file = paste0(curr_direc, "UpdatedDataChange", change, "GroupCombo", GroupNum, "InfRepNum", curr_infRep, ".RData"))

    rm(list = paste0("UpdatedDataChange", change, "GroupCombo", GroupNum, "InfRepNum", curr_infRep))

    rm(new_dat7)
    gc()

    print(paste0("Time for infRep num ", j, " is ", proc.time()[3] - sdd[3]))
  } #Ending j loop (indexing over number of infreps to save abGene like files for each infRep)
  
  
  print(paste0("Total time to read in updated data is ", proc.time()[3] - startT[3]))

  
  genes_part <- filteredgenenames_split
  for(k in 1:nparts){
    #Now, split each infrep updated dataset into parts and
    part <- k
    print(paste0("Part is ", part))
    
    if(generateNewIlrMeansCovs==TRUE){
      if(file.exists(paste0(curr_direc, "NewilrMeansCovsChange", change, "GroupCombo", GroupNum, "Part", part, ".RData"))){
        next
      }
    }

    
    Sys.sleep(sample(1:10, 1))
    print(gc())

    
    
    print(paste0("Length of genes_part for the current part (Part ", part, ") is ", length(genes_part[[k]])))
    print(paste0("Head of current genes_part is ", head(genes_part[[k]])))
    #genes_part will only have length 0 if the code is being rerun after not finishing the first time
    #and the code is rerun on the cases that didn't finish
    #In this case need to regenerate genes_part for the parts that didn't finish
    # if(length(genes_part)==0){
    #   genes_part <- list()
    #   for(u in part:nparts){
    #     print(paste0("This must be a rerun of code that didn't finish. genes_part is being calculated for part ", u, " out of ", nparts))
    #     Sys.sleep(sample(1:10, 1))
    #     load(paste0(load_dir, "FullinfRepDatPart", u, ".RData"))
    #     #The subset of genes wont differ based on the current infrep being used, so just restrict it to 1 to
    #     #get the list needed
    #     curr_part_dat <- subset(datt, datt$infRepNum==1)
    #     gene_names_temp <- sort(unique(curr_part_dat$gene_id))
    #     genes_part[[u]] <- intersect(gene_names_temp, filteredgenenames_split[[u]])
    #     print(paste0("length of genenames for part ", u, " is ", length(genes_part[[u]])))
    #     rm(datt)
    #     rm(curr_part_dat)
    #     rm(gene_names_temp)
    #     gc()
    #   }
    # 
    # }
    
    if(!file.exists(paste0(curr_direc, "NewabDatasetsGibbsChange", change, "GroupCombo", GroupNum, "Part", part, ".RData"))){
      
      for(j in 1:ninfreps){
        #print(paste0("Current infrep number is ", j))
        curr_infRep <- j
        load(paste0(curr_direc, "UpdatedDataChange", change, "GroupCombo", GroupNum, "InfRepNum", curr_infRep, ".RData"))
        curr_data <- get(paste0("UpdatedDataChange", change, "GroupCombo", GroupNum, "InfRepNum", curr_infRep))
        rm(list = paste0("UpdatedDataChange", change, "GroupCombo", GroupNum, "InfRepNum", curr_infRep))
        new_dat7 <- subset(curr_data, curr_data$gene_id %in% genes_part[[k]])
        updated_data[[j]] <- new_dat7
        rm(new_dat7)
        rm(curr_data)
        gc()
      }
    
      print(gc())
      updated_data2 <- rbindlist(updated_data)

      abGNew <- lapply(genes_part[[k]], generateData, dat = updated_data2, nsamp = nrow(key_to_use), abundance = TRUE,
                       abCompDatasets = abDatasetsFiltered, useOtherGroups = useOtherGroups, useExistingOtherGroups = TRUE,
                       samps = samps, infReps = infReps, ninfreps = ninfreps, CompMI = TRUE)
      
      
      names(abGNew) <- genes_part[[k]]
      
      assign(paste0("NewabDatasetsGibbsChange", change, "GroupCombo", GroupNum, "Part", part), abGNew)
      
      save(list = paste0("NewabDatasetsGibbsChange", change, "GroupCombo", GroupNum, "Part", part),
           file = paste0(curr_direc, "NewabDatasetsGibbsChange", change, "GroupCombo", GroupNum, "Part", part, ".RData"))
      rm(list = paste0("NewabDatasetsGibbsChange", change, "GroupCombo", GroupNum, "Part", part))
      rm(updated_data2)
      gc()
      print(gc())
    }else{
      if(generateNewIlrMeansCovs==TRUE){
        abGNew <- loadRData(paste0(curr_direc, "NewabDatasetsGibbsChange", change, "GroupCombo", GroupNum, "Part", part, ".RData"))
      }
    }
    
    

    print(paste0("Saving of new abDatasets with Gibbs Reps is complete for part ", part))

    #newAbDatasetsGibbs[[k]] <- subset(abGNew, names(abGNew) %in% filteredgenenames)
    if(generateNewIlrMeansCovs==TRUE){
      NewilrMCovs <- lapply(genes_part[[k]], calcIlrMeansCovs, dat = abGNew, nsamp = length(Group), CLE = TRUE)
      #NewilrMCovs <- laply(genes_part[[k]], calcIlrMeansCovs, dat = abGNew, nsamp = length(Group), .progress = "text")
      names(NewilrMCovs) <- genes_part[[k]]
      
      
      assign(paste0("NewilrMeansCovsChange", change, "GroupCombo", GroupNum, "Part", part), NewilrMCovs)
      
      save(list = paste0("NewilrMeansCovsChange", change, "GroupCombo", GroupNum, "Part", part),
           file = paste0(curr_direc, "NewilrMeansCovsChange", change, "GroupCombo", GroupNum, "Part", part, ".RData"))
      #NewilrMeansCovs[[k]] <- subset(NewilrMCovs, names(NewilrMCovs) %in% filteredgenenames)
      
      
      rm(list = paste0("NewilrMeansCovsChange", change, "GroupCombo", GroupNum, "Part", part))
      rm(list = paste0("NewabDatasetsGibbsChange", change, "GroupCombo", GroupNum, "Part", part))
      rm(NewilrMCovs)
    }

    rm(abGNew)
    gc()
    print(gc())


  }


  # NewilrMeansCovsFinal <- NewilrMeansCovs[[1]]
  # newAbDatasetsGibbsFinal <- abGNew[[1]]
  # for(i in 2:nparts){
  #   NewilrMeansCovsFinal <- append(NewilrMeansCovsFinal, NewilrMeansCovs[[i]])
  #   newAbDatasetsGibbsFinal <- append(newAbDatasetsGibbsFinal, abGNew[[i]])
  # }


  #print(paste0("total computation time is ", proc.time() - stT))
}



#Need to get the "bootstrap" samples ready to input into this function
  #I think it will be easier to start from "GibbsSampsSample1" type files instead of any more processed files
  #But, maybe the size of the necessary objects with many samples would make this difficult if not impossible?
RATsObsAnalysis <- function(cntGene, cntDatasets, key, infRepDataA, infRepDataB, fullReturn = FALSE, DRIMSeqFiltering){

  #browser()
  if(DRIMSeqFiltering==TRUE){
    samps <- key$Identifier
    cntDataTemp <- cntGene[,c("gene_id", "tx_id", paste0(samps, "Cnt"))]
    rm(samps)
  }else if(DRIMSeqFiltering==FALSE){
    cntDataTemp <- convertcntGeneToOtherGroupsForDRIMSeq(cntGene = cntGene, cntDatasets = cntDatasets, samps = key$Identifier)
    print("Conversion of cntGene to OtherGroups is Complete")
  }
  gc()
  print(gc())
  


  samps_condA <- key$Identifier[key$Condition==levels(key$Condition)[1]]
  samps_condB <- key$Identifier[key$Condition==levels(key$Condition)[2]]

  cnt_data_A <- cntDataTemp[,c("tx_id", paste0(samps_condA, "Cnt"))]
  cnt_data_B <- cntDataTemp[,c("tx_id", paste0(samps_condB, "Cnt"))]
  
  print(paste0("The first 10 samps in condition A are below"))
  print(head(samps_condA, 10))
  
  print(paste0("The first 10 samps in condition B are below"))
  print(head(samps_condB, 10))
  
  print(paste0("The dimensions of the cnt_data_A object are given below (should be ntranscripts by (nsamp in current analysis + 1)"))
  print(dim(cnt_data_A))
  
  print(paste0("The dimensions of the cnt_data_B object are given below (should be ntranscripts by (nsamp in current analysis + 1)"))
  print(dim(cnt_data_B))

  colnames(cnt_data_A)[colnames(cnt_data_A)=="tx_id"] <- "target_id"
  colnames(cnt_data_B)[colnames(cnt_data_B)=="tx_id"] <- "target_id"

  #Need to take the transcript names from here because the presence of "other" groups means
  # the annotation in tx2gene will not be what is needed here
  annot <- cntDataTemp[,c("gene_id", "tx_id")]
  colnames(annot) <- c("parent_id", "target_id")

  for(i in 1:length(infRepDataA)){
    colnames(infRepDataA[[i]])[colnames(infRepDataA[[i]])=="feature_id"] <- "target_id"
  }

  for(i in 1:length(infRepDataB)){
    colnames(infRepDataB[[i]])[colnames(infRepDataB[[i]])=="feature_id"] <- "target_id"
  }
  #annot2 <- subset(annot, annot$target_id %in% infRepDataA[[1]]$target_id)
  count_data_A <- data.table(cnt_data_A)
  count_data_B <- data.table(cnt_data_B)

  print(paste0("The Number of transcripts used in the non-bootstrap analysis is ", nrow(count_data_A)))
  gc()
  print(gc())

  StartTime1 <- proc.time()
  rats_res <- call_DTU(annot=annot, count_data_A=count_data_A, count_data_B=count_data_B, qboot = FALSE)
  RatsCompTime <- proc.time() - StartTime1
  gc()
  print(gc())
  print("RATs run on the non-infRep data is complete")

  StartTime2 <- proc.time()
  
  #Remove any extra transcripts that are in the inferential replicate data but not in the regular count data
  infRepDataA2 <- list()
  for(i in 1:length(infRepDataA)){
    infRepDataA2[[i]] <- subset(infRepDataA[[i]], infRepDataA[[i]]$target_id %in% annot$target_id) 
  }
  
  infRepDataB2 <- list()
  for(i in 1:length(infRepDataB)){
    infRepDataB2[[i]] <- subset(infRepDataB[[i]], infRepDataB[[i]]$target_id %in% annot$target_id) 
  }
  
  rm(infRepDataA)
  rm(infRepDataB)
  gc()

  print(paste0("The Number of transcripts used in the bootstrap analysis is ", nrow(infRepDataA2[[1]])))
  rats_res_infReps <- call_DTU(annot=annot, boot_data_A = infRepDataA2, boot_data_B = infRepDataB2)

  RatsCompTimeInfReps <- proc.time() - StartTime2
  gc()
  print(gc())
  print("RATs run on the infRep data is complete")

  res_reg <- data.frame(rats_res$Genes$parent_id, rats_res$Genes$pval, rats_res$Genes$pval_corr, as.numeric(rats_res$Genes$DTU),
                        as.numeric(rats_res$Genes$elig_fx), as.numeric(rats_res$Genes$rep_reprod))
  colnames(res_reg) <- c("gene_id", "pval", "adj_pvalue", "gene_DTU", "elig_fx", "rep_reprod")


  res_infRep <- data.frame(rats_res_infReps$Genes$parent_id, rats_res_infReps$Genes$pval,
                           rats_res_infReps$Genes$pval_corr, as.numeric(rats_res_infReps$Genes$DTU),
                           as.numeric(rats_res_infReps$Genes$elig_fx), as.numeric(rats_res_infReps$Genes$quant_reprod),
                           as.numeric(rats_res_infReps$Genes$rep_reprod))
  colnames(res_infRep) <- c("gene_id", "pval_infReps", "adj_pval_infReps", "gene_DTU_infReps", "elig_fx_infReps", "quant_reprod_infReps", "rep_reprod_infReps")

  res <- merge(res_reg, res_infRep, by = "gene_id", all = TRUE)

  if(fullReturn==TRUE){
    return(list(rats_res = rats_res, rats_res_infReps = rats_res_infReps, res = res, RatsCompTime = RatsCompTime, RatsCompTimeInfReps = RatsCompTimeInfReps))
  }else{
    return(list(res = res, RatsCompTime = RatsCompTime, RatsCompTimeInfReps = RatsCompTimeInfReps))
  }

}



#x is the gene name
calcOtherGroupValsForRATs <- function(x, curr, cntDatasets, curr_samp){
  sub_curr <- subset(curr, curr$gene_id==x)
  sub_cntDatasets <- cntDatasets[[x]]$Counts
  if(is.null(sub_cntDatasets)){
    return(NULL)
  }

  OtherTrans <- attr(sub_cntDatasets, "OtherTrans")

  if(length(OtherTrans)==0){
    return(sub_curr)
  }
  cntCol <- paste0(curr_samp, "Cnt")

  sub_curr_OtherTrans <- subset(sub_curr, sub_curr$tx_id %in% OtherTrans)
  OtherTransVals <- aggregate(sub_curr_OtherTrans[,cntCol], FUN = sum, by = list(infRepNum = sub_curr_OtherTrans$infRepNum))
  OtherTransVals$gene_id <- x
  OtherTransVals$tx_id <- paste0("Other", OtherTransVals$gene_id)
  colnames(OtherTransVals) <- c("infRepNum", cntCol, "gene_id", "tx_id")

  OtherTransVals2 <- OtherTransVals[,c("tx_id", cntCol, "infRepNum", "gene_id")]
  sub_curr2 <- subset(sub_curr, !(sub_curr$tx_id %in% OtherTrans))

  sub_curr3 <- rbind(sub_curr2, OtherTransVals2)
  sub_curr4 <- data.table(sub_curr3)
  return(sub_curr4)
}

#y is the current sample to use
reshapeinfRepDataforRATs <- function(y, cntDatasets, cntGene, infReps, infReps_dir, useOtherGroups){
  curr_samp <- y
  print(paste0("Currently Processing ", curr_samp))

  if(infReps=="GibbsThin100"){
    curr_temp <- loadRData(paste0(infReps_dir, "GibbsSamps", curr_samp, ".RData"))
  }else if(infReps=="Boot"){
    curr_temp <- loadRData(paste0(infReps_dir, "BootSamps", curr_samp, ".RData"))
  }


  cntCol <- paste0(curr_samp, "Cnt")
  curr_data <- data.frame(curr_temp)


  curr <- curr_data[,c("tx_id", cntCol, "infRepNum", "gene_id")]

  
  # OtherGroupDatList <- laply(names(cntDatasets), calcOtherGroupValsForRATs, curr = curr, cntDatasets = cntDatasets, curr_samp = curr_samp,
  #                            .inform = TRUE, .progress = "text")

  if(useOtherGroups==TRUE){
    stop("If this is for the power analysis the counts left needs to include all counts, not just for trasncripts that pass filtering.  If you still want to do this you will have to check code and ensure it is doing what you want it to.")
    #OtherGroupDatList <- lapply(names(cntDatasets), calcOtherGroupValsForRATs, curr = curr, cntDatasets = cntDatasets, curr_samp = curr_samp)
    #infRepDataT <- rbindlist(OtherGroupDatList)
    infRepDataT <- subset(curr, curr$tx_id %in% cntGene$tx_id)
  }else if(useOtherGroups==FALSE){
    infRepDataT <- subset(curr, curr$tx_id %in% cntGene$tx_id)
  }


  t3 <- spread(infRepDataT, key = "infRepNum", value = cntCol)
  t3$gene_id <- NULL

  colnames(t3)[colnames(t3)=="tx_id"] <- "feature_id"

  if(is.data.table(t3)==FALSE){
    t4 <- data.table(t3)
  }else{
    t4 <- t3
  }

  return(t4)
}



calcScaledCountsRATsPower <- function(curr_col_curr_dat_changed, curr_len2, curr_samp){
  curr_col_len <- curr_len2[,paste0(curr_samp, "Len")]
  
  z <- curr_col_curr_dat_changed/curr_col_len
  
  ScaledTPMCounts <- (z/sum(z)) * sum(curr_col_curr_dat_changed)
  return(ScaledTPMCounts) 
}



CorrectLowExpression <- function(y, CLEParam = 0.05){
  if(is.null(y)){
    return(NULL)
  }
  
  if(ncol(y)==1){
    y[y==0] <- CLEParam
    return(y)
  }else{
    curr_dat <- y
    
    lowExp_corrected_dat <- t(apply(curr_dat, 1, CorrectLowExpressionHelper, CLEParam = CLEParam))
    return(lowExp_corrected_dat)
  }
  
  
}




CorrectLowExpressionHelper <- function(x, CLEParam = 0.05){
  #print(x)
  # if(nrow(x)!=1){
  #   stop("Number of rows is not 1")
  # }
  curr_rowSum <- sum(x)
  if(is.na(curr_rowSum)){
    return(x)
  }else if(curr_rowSum==0){
    return(x)
  }else{
    curr_dat2 <- x
    curr_dat2[curr_dat2 < CLEParam*curr_rowSum] <- CLEParam*curr_rowSum
    return(curr_dat2)
  }
}




calcInfRVRes <- function(type = "Cnt", nsamp, samps, abGeneFiltered, tx2gene, GibbsSamps, def_wd1){
  
  res <- vector(mode = "list", length = nsamp)
  for(i in 1:nsamp){
    sts <- proc.time()
    print(paste0("Current sample being processed is ", i, " out of ", nsamp))
    curr_sampT <- samps[i]
    curr_samp <- as.numeric(strsplit(curr_sampT, "Sample")[[1]][2])
    
    if(GibbsSamps==TRUE){
      curr_datT <- loadRData(paste0(def_wd1,"GibbsSampsSample", curr_samp, ".RData"), objNameToGet = paste0("GibbsSampsSample", curr_samp))
    }else{
      curr_datT <- loadRData(paste0(def_wd1,"BootSampsSample", curr_samp, ".RData"), objNameToGet = paste0("BootSampsSample", curr_samp))
    }
    
    #Restrict to only those transcripts that pass filtering
    curr_dat <- data.frame(subset(curr_datT, curr_datT$tx_id %in% abGeneFiltered$tx_id))
    
    rm(curr_datT)
    gc()
    #cooljerk
    gene_names <- unique(curr_dat$gene_id)
    #st1 <- proc.time()
    if(type=="ilr"){
      res_ilr <- vector(mode = "list", length = length(gene_names))
      for(g in 1:length(gene_names)){
        #if(g %% 100==0){
        #  print(paste0("Current gene number is ", g, " out of ", length(gene_names)))
        #}
        curr_g <-gene_names[g]
        
        curr_dat_sub <- subset(curr_dat, curr_dat$gene_id==curr_g)
        curr_dat_sub$gene_id <- NULL
        curr_dat_sub[,paste0("Sample", curr_samp, "Cnt")] <- NULL
        
        curr_dat_sub2 <- dcast(curr_dat_sub, infRepNum ~ tx_id, value.var = paste0("Sample", curr_samp, "TPM"))
        curr_dat_sub2$infRepNum <- NULL
        
        if(dim(curr_dat_sub2)[2]==1){
          ff1 <- data.frame(gene_id = curr_g, GeneMaxVarInfRepilrScale = NA, GeneMeanVarInfRepilrScale = NA, 
                            GeneMedianVarInfRepilrScale = NA, GeneMaxVarPropScale = NA,
                            GeneMeanVarPropScale = NA,
                            GeneMedianVarPropScale = NA, Sample = curr_samp)
          
          res_ilr[[g]] <- ff1
          next
        }else if(dim(curr_dat_sub2)[2]==2){
          props <- ccomp(CorrectLowExpression(curr_dat_sub2), total = 1)
          props_variances <- diag(cov(props))
          #props_MAD <- apply(props, 2, FUN = mad)
          
          ilr_vals <- apply(props, 1, FUN = ilr)
          ilr_vals_variances <- var(ilr_vals, na.rm = T)
          ilr_vals_means <- mean(ilr_vals, na.rm = T)
        }else{
          props <- ccomp(CorrectLowExpression(curr_dat_sub2), total = 1)
          props_variances <- diag(cov(props))
          #props_MAD <- apply(props, 2, FUN = mad)
          
          ilr_vals <- t(apply(props, 1, FUN = ilr))
          ilr_vals_variances <- diag(cov(ilr_vals))
          ilr_vals_means <- colMeans(ilr_vals, na.rm = T)
        }
        
        
        
        
        ff1 <- data.frame(gene_id = curr_g, GeneMaxVarInfRepilrScale = max(ilr_vals_variances, na.rm = T), GeneMeanVarInfRepilrScale = mean(ilr_vals_variances, na.rm = T), 
                          GeneMedianVarInfRepilrScale = median(ilr_vals_variances, na.rm = T), 
                          GeneMaxVarPropScale = max(props_variances, na.rm = T),
                          GeneMeanVarPropScale = mean(props_variances, na.rm = T),
                          GeneMedianVarPropScale = median(props_variances, na.rm = T),
                          Sample = curr_samp)
        
        res_ilr[[g]] <- ff1
      }
      
      res[[i]] <- data.frame(data.table::rbindlist(res_ilr, fill = TRUE))
      cttime <- (proc.time() - sts)[3]
      print(paste0("Total time for Sample ", i, " is ", cttime , " seconds"))
      next
    }
    #ct <- proc.time() - st1
    #print(ct)
    
    #The code below here up to the next if loop is only run if the type is not ilr
      #Since the ilr coordinates are based on TPM values and are transformed, calculating InfRV based on the ilr coordinates doesn't make sense
    
    curr_ab_col <- paste0("Sample", curr_samp, type)
    #curr_ab_col <- paste0("Sample", curr_samp, "TPM")
    
    SampVars <- aggregate(curr_dat[,curr_ab_col], by = list(tx_id = curr_dat$tx_id), FUN = var, na.rm = T)
    colnames(SampVars) <- c("tx_id", "SampVars")
    
    SampMeans <- aggregate(curr_dat[,curr_ab_col], by = list(tx_id = curr_dat$tx_id), FUN = mean, na.rm = T)
    colnames(SampMeans) <- c("tx_id", "SampMeans")
    
    SampsMeansAndVars <- merge(SampMeans, SampVars, by = "tx_id")
    
    #numerator <- SampsMeansAndVars[,"SampVars", drop = FALSE] - SampsMeansAndVars[,"SampMeans", drop = FALSE]
    #numerator[numerator < 0] <- 0
    #denom <- SampsMeansAndVars[,"SampMeans", drop = FALSE] + 5
    
    
    numerator <- SampsMeansAndVars$SampVars- SampsMeansAndVars$SampMeans
    numerator[numerator < 0] <- 0
    denom <- SampsMeansAndVars$SampMeans + 5
    
    
    InfRV <- (numerator/ denom) + 0.01
    #colnames(InfRV) <- "InfRV"
    SampsMeansAndVars$InfRV <- InfRV
    
    SampsMeansAndVars2T <- merge(SampsMeansAndVars, tx2gene, by = "tx_id")
    
    SampsMeansAndVars2 <- SampsMeansAndVars2T[order(SampsMeansAndVars2T$gene_id, SampsMeansAndVars2T$tx_id),]
    
    SampsMeansAndVarsGeneMax <- aggregate(SampsMeansAndVars2$InfRV, by = list(gene_id = SampsMeansAndVars2$gene_id), FUN = max, na.rm = T)
    colnames(SampsMeansAndVarsGeneMax) <- c("gene_id", "GeneMaxInfRV")
    
    
    SampsMeansAndVarsGeneMedian <- aggregate(SampsMeansAndVars2$InfRV, by = list(gene_id = SampsMeansAndVars2$gene_id), FUN = median, na.rm = T)
    colnames(SampsMeansAndVarsGeneMedian) <- c("gene_id", "GeneMedianInfRV")
    
    SampsMeansAndVarsGeneMean <- aggregate(SampsMeansAndVars2$InfRV, by = list(gene_id = SampsMeansAndVars2$gene_id), FUN = mean, na.rm = T)
    colnames(SampsMeansAndVarsGeneMean) <- c("gene_id", "GeneMeanInfRV")
    
    
    t1 <- merge(SampsMeansAndVarsGeneMax, SampsMeansAndVarsGeneMedian, by = "gene_id")
    t2 <- merge(t1, SampsMeansAndVarsGeneMean, by = "gene_id")
    
    FinalInfRV <- t2
    FinalInfRV$Sample <- curr_samp
    res[[i]] <- FinalInfRV
    cttime <- (proc.time() - sts)[3]
    print(paste0("Total time for Sample ", i, " is ", cttime , " seconds"))
    
    rm(curr_dat)
    gc()
  }
  
  if(type=="ilr"){
    InfRepVar <- data.frame(data.table::rbindlist(res, fill = TRUE))
    
    temp1 <- aggregate(InfRepVar$GeneMaxVarInfRepilrScale, by = list(gene_id = InfRepVar$gene_id), FUN = mean, na.rm = T)
    colnames(temp1) <- c("gene_id", "MeanGeneMaxVarInfRepilrScale")
    temp2 <- aggregate(InfRepVar$GeneMedianVarInfRepilrScale, by = list(gene_id = InfRepVar$gene_id), FUN = mean, na.rm = T)
    colnames(temp2) <- c("gene_id", "MeanGeneMedianVarInfRepilrScale")
    temp3 <- aggregate(InfRepVar$GeneMeanVarInfRepilrScale, by = list(gene_id = InfRepVar$gene_id), FUN = mean, na.rm = T)
    colnames(temp3) <- c("gene_id", "MeanGeneMeanVarInfRepilrScale")
    temp4 <- aggregate(InfRepVar$GeneMaxVarPropScale, by = list(gene_id = InfRepVar$gene_id), FUN = mean, na.rm = T)
    colnames(temp4) <- c("gene_id", "MeanGeneMaxVarPropScale")
    temp5 <- aggregate(InfRepVar$GeneMeanVarPropScale, by = list(gene_id = InfRepVar$gene_id), FUN = mean, na.rm = T)
    colnames(temp5) <- c("gene_id", "MeanGeneMeanVarPropScale")
    temp6 <- aggregate(InfRepVar$GeneMedianVarPropScale, by = list(gene_id = InfRepVar$gene_id), FUN = mean, na.rm = T)
    colnames(temp6) <- c("gene_id", "MeanGeneMedianVarPropScale")
    
    
    InfRepVarRes <- Reduce(function(x, y) {merge(x, y, by = "gene_id")}, list(temp1, temp2, temp3, temp4, temp5, temp6))
    #m1 <- merge(temp1, temp2, by  = "gene_id")
    #m2 <- merge(m1, temp3, by  = "gene_id")
    
    
    #InfRepVarRes <- m2
    return(InfRepVarRes)
    
  }else{
    InfRV <- data.frame(data.table::rbindlist(res, fill = TRUE))
    
    temp1 <- aggregate(InfRV$GeneMaxInfRV, by = list(gene_id = InfRV$gene_id), FUN = mean, na.rm = T)
    colnames(temp1) <- c("gene_id", "MeanGeneMaxInfRV")
    temp2 <- aggregate(InfRV$GeneMedianInfRV, by = list(gene_id = InfRV$gene_id), FUN = mean, na.rm = T)
    colnames(temp2) <- c("gene_id", "MeanGeneMedianInfRV")
    temp3 <- aggregate(InfRV$GeneMeanInfRV, by = list(gene_id = InfRV$gene_id), FUN = mean, na.rm = T)
    colnames(temp3) <- c("gene_id", "MeanGeneMeanInfRV")
    
    
    m1 <- merge(temp1, temp2, by  = "gene_id")
    m2 <- merge(m1, temp3, by  = "gene_id")
    
    
    InfRVRes <- m2
    return(InfRVRes)
  }
}




GeneratePlots <- function(infReps = "Boot", DRIMSeqFiltering = TRUE, SuppPlots = FALSE, UseRealSensRes = FALSE, FilterGenes = FALSE, FilterCutoff = 0.50,
                          changestouse, PlotTypes, stdal = F, TwentySamplesOnly = FALSE, def_wd, def_wd2, save_dir, abDatasets,
                          AbFromInfReps = NULL, geneGroupings = NULL,
                          saveCompiledPDFs = FALSE){
  
  if(saveCompiledPDFs==TRUE){
    stdal <- TRUE
    print("To save pre-compiled pdf files of the plots, the stdal option has been set to TRUE (even if it was input to function as FALSE)")
  }
  
  if(is.null(geneGroupings)){
    stop("geneGroupings (ex sets of genes to stratify plots based on) must be specified")
  }
  types <- geneGroupings
  
  if(is.null(AbFromInfReps)){
    AbFromInfReps <-- FALSE
  }
  
  
  if(DRIMSeqFiltering==TRUE){
    load("abDatasetsNoOtherGroupsFiltered.RData")
    abDatasets <- abDatasetsFiltered
  }else if(DRIMSeqFiltering==FALSE){
    load("abDatasets.RData")
  }
  
  if(TwentySamplesOnly==TRUE){
    dirr_modifier <- "TwentySamples"
  }else{
    dirr_modifier <- ""
  }
  
  print(paste0("Current change is ", changestouse))
  
    if(infReps=="Gibbs" & DRIMSeqFiltering==FALSE){
      if(AbFromInfReps == TRUE){
        dirr <- paste0("/ROCAllResAbFromInfReps", dirr_modifier, "/")
      }else if(SuppPlots==TRUE){
        dirr <- paste0("/ROCAllResSupp", dirr_modifier, "/")
      }else if(SuppPlots==FALSE){
        dirr <- paste0("/ROCAllRes", dirr_modifier, "/")
      }
      
      if(FilterGenes==TRUE){
        direc <- paste0(def_wd, "GEUV1PowerRes/", "AllResultsFilter", FilterCutoff, dirr_modifier, "/")
      }else{
        direc <- paste0(def_wd, "GEUV1PowerRes/", "AllResults", dirr_modifier, "/")
      }
    }else if(infReps=="Boot" & DRIMSeqFiltering==FALSE){
      if(AbFromInfReps == TRUE){
        dirr <- paste0("/ROCAllResBootAbFromInfReps", dirr_modifier, "/")
      }else if(SuppPlots==TRUE){
        dirr <- paste0("/ROCAllResBootSupp", dirr_modifier, "/")
      }else if(SuppPlots==FALSE){
        dirr <- paste0("/ROCAllResBoot", dirr_modifier, "/")
      }
      
      if(FilterGenes==TRUE){
        direc <- paste0(def_wd, "GEUV1PowerResBoot/", "AllResultsFilter", FilterCutoff, dirr_modifier, "/")
      }else{
        direc <- paste0(def_wd, "GEUV1PowerResBoot/", "AllResults", dirr_modifier, "/")
      }
    }else if(infReps=="Boot" & DRIMSeqFiltering==TRUE){
      if(AbFromInfReps == TRUE){
        dirr <- paste0("/ROCAllResBootDRIMSeqFilteringAbFromInfReps", dirr_modifier, "/")
      }else if(SuppPlots==TRUE){
        dirr <- paste0("/ROCAllResBootDRIMSeqFilteringSupp", dirr_modifier, "/")
      }else if(SuppPlots==FALSE){
        dirr <- paste0("/ROCAllResBootDRIMSeqFiltering", dirr_modifier, "/")
      }
      
      if(FilterGenes==TRUE){
        direc <- paste0(def_wd, "GEUV1PowerResBootDRIMSeqFiltering/", "AllResultsFilter", FilterCutoff, dirr_modifier, "/")
      }else{
        direc <- paste0(def_wd, "GEUV1PowerResBootDRIMSeqFiltering/", "AllResults", dirr_modifier, "/")
      }
    }
    
  
  
  print(paste0("dirr is ", dirr))
  fil <- list.files(direc, pattern = ".RData", full.names = TRUE)
  for (i in 1:length(fil)){
    load(fil[i])
  }
  
  
  
  for (i in 1:length(changestouse)){
    curr_change <- changestouse[i]
    fi <- ls(pattern = paste0("ResListChange", curr_change))
    for(j in 1:length(fi)){
      curr_list <- get(fi[j])
      list2env(curr_list, envir = .GlobalEnv)
    }
  }
  
  #genesGenes2Columns <- genes2Columns
  #genesGenesWithNoOtherGroups <- genesWithNoOtherGroups
  #genesGenesWithNoOtherGroupsAnd2Columns <- genesWithNoOtherGroupsAnd2Columns
  
  #genesGenesWithNoOtherGroupsHighTGE <- genesWithNoOtherGroupsHighTGE
  #genesGenesWithNoOtherGroupsAnd2ColumnsHighTGE <- genesWithNoOtherGroupsAnd2ColumnsHighTGE
  
  
  
  #Total number of genes used in this analysis (Will differ based on whether FilterGenes is turned on or not)
  ngenestotal <- length(genestouse)
  genesFull <- genestouse
  
  #Set number of digits to use with fr function
  ndig <- 3
  
  # Force correct number of digits and round properly for in text references and table columns too
  # And insert <0.0001 for small p-values with number of zeros after decimal controlled by digits option
  # Updated 3-6-17
  fr <- function(x, n = getOption("digits")) {
    y <- rep(0, length(x))
    for (i in 1:length(x)) {
      if (x[i] >=0 & x[i]< 10^-(n) & !is.na(x[i])){
        if(n==0){n <- 1}
        y[i] <- paste("$<$", ".", paste(rep("0", n-1), collapse=""), "1", sep = "")
      }
      if (is.numeric(x[i]) & !is.na(x[i]) &  (x[i] >= 10^-(n) | x[i] < 0)) {
        y[i] <- format(round(x[i], digits = n), nsmall = n, scientific = FALSE)
      }
      # Only truly missing values should ever show up as 0 in y created above, so change them back to missing
      #Note that this is not the same as having a value of 0 in the input vector (x), as these will correctly
      #Just be formatted to show as <0.0001, etc
      if(y[i] == "0") {
        y[i] <- " "
      }
    }
    return(y)
  }
  
  if(SuppPlots==T){
    options("tikzDocumentDeclaration" = "\\documentclass[10pt]{article}\n")
  }else{
    options("tikzDocumentDeclaration" = "\\documentclass[10pt]{article}\n")
  }
  
  options("tikzMetricsDictionary" = paste(def_wd2, "/TikzTemp", sep = ""))
  
  if(AbFromInfReps==TRUE){
    
    methdstouse <- c("CompDTU", "CompDTUAbRowMeanInfRepBoot", "CompDTUme")
    methdnamestouse <- c("CompDTU", "CompDTUAbMeanBoot", "CompDTUme")
    colrs <- c("blue", "darkgreen", "red")
    

    lins <- c(1,1,1,1,1,1,1)
  }else if(SuppPlots==TRUE){
    methdstouse <- c("CompDTU", "CompDTUme", "DRIMSeqAdd1ToEveryCount","RATs", "RATsInfReps", "CompMI", "CompMICombineCoefs")
    methdnamestouse <- c("CompDTU", "CompDTUme", "DRIMSeq","RATsNoBoot", "RATsBoot", "CompMICombinePvals", "CompMICombineCoefs")
    colrs <- c("blue", "red", "black", "purple", "orange", "green", "pink", "darkgreen")
    lins <- c(1,1,1,1,1,1,1,1)
  }else if(SuppPlots==FALSE){
    methdstouse <- c("CompDTU", "CompDTUme", "DRIMSeqAdd1ToEveryCount", "RATs", "RATsInfReps")
    methdnamestouse <- c("CompDTU", "CompDTUme", "DRIMSeq", "RATsNoBoot", "RATsBoot")
    
    
    colrs <- c("blue", "red", "black", "purple", "orange", "green", "pink", "darkgreen")
    lins <- c(1,1,1,1,1,1,1,1)
  }
  
  RATsmethds <- c("RATs", "RATsInfReps")
  
  
  #Currently:
  #type 1 is the usual ROC curve (fpr vs sens),
  #type 2 is pvalcutoff vs FPR
  #type 3 is pvalcutoff vs sens
  #type 4 is fdr vs sens (modified to use eFDR adjusted pvalues for now)
  #type 5 is pvalcutoff vs fdr
  #type 6 is pvalue cutoff vs accuracy
  #type 7 is  fdr vs accuracy
  #type 8 is fdr adjusted pvalue cutoff (eFDR) vs accuracy
  #type 9 is fdr adjusted pvalue cutoff (eFDR) vs sens (so similar to type 4)
  #type 10 is fpr sens based on fdr adjusted pvalues (so similar to usual roc curve (type 1))
  #type 11 is fdr adjusted pvalue cutoff (eFDR) vs true FDR- so, should be a straight line if everything is
  #working properly
  #type 12 is fdr adjusted pvalue cutoff (eFDR) vs FPR- so the same as type 2 but using adjusted pvalues this time
  #This type might not really make sense to use under the null since type I error definition is just that
  #We would expect only 5 percent of rejections to be false (with no mention or need to adjust the pvals)
  #type 13 is fdr adjusted pvalue cutoff (eFDR) vs true negative rate (ie 1 - fpr)
  #type 14 is true FDR vs Sensitivity
  #type 15 is eFDR vs "Real Sensitivity"
  #plotType <- 1
  #changestouse <- changes
  #changestouse <- changes[changes==0.01 | changes==2]
  #changestouse <- changes[changes==4]
  #changestouse <- 1
  #1:14
  for(z in PlotTypes){
    if(UseRealSensRes==TRUE & z < 15){
      stop("Don't use RealSensRes with z <=14")
    }
    #Set the RATs pch type to 18 because all RATs results correspond to an eFDR cutoff of 0.05 and a diamond corresponds to what is done
    #for the other methods
    #Don't set this for plot 2 since eFDR doesn't really make sense with the unadjusted pvalues used to make this plot
    #And don't use for type 11 either since eFDR is on the x axis already so a simple circle will do fine
    RATspch <- 18
    #Make diamond larger to match size of the one in other figures
    RATscex <- 2.5
    
    plotType <- z
    print(paste0("Current plotType is ", plotType))
    
    for(j in 1:length(types)){
      type <- types[j]
      for (i in 1:length(changestouse)){
        change <- changestouse[i]
        upperx <- 0.25
        #upperx <- 0.30
        xlim1 <- c(0,upperx)
        #if(change==1)
        if(change == 1){
          ylim1 <- c(0,0.25)
        }else{
          if(plotType==1 | plotType==3 | plotType==4 | plotType==6 | plotType==7 | plotType==8 | plotType==9 | plotType==10 | plotType==13 | plotType==14 | plotType==15){
            #ylim1 <- c(0.75, 1)
            ylim1 <- c(0,1)
          }else{
            ylim1 <- c(0,0.25)
          }
        }
        
        if(plotType==8 | plotType==9){
          eFDRvals <- list()
        }
        
        for (k in 1:length(methdstouse)){
          curr_methd <- methdstouse[k]
          
          if(type=="Full"){
            if(UseRealSensRes==TRUE){
              assign(paste0(curr_methd, "Res", change, type), get(paste0("SensFpr", curr_methd, change, "UpdatedToCalcRealSens")))
            }else{
              assign(paste0(curr_methd, "Res", change, type), get(paste0("SensFpr", curr_methd, change)))
            }
            
          }else{
            if(UseRealSensRes==TRUE){
              assign(paste0(curr_methd, "Res", change, type), get(paste0("SensFpr", curr_methd, change, type, "UpdatedToCalcRealSens")))
            }else{
              assign(paste0(curr_methd, "Res", change, type), get(paste0("SensFpr", curr_methd, change, type)))
            }
            
          }
          CurrResT <- get(paste0(curr_methd, "Res", change, type))
          CurrRes <- CurrResT[order(CurrResT$pvalcutoff),]
          assign(paste0("p", curr_methd, "used"), fr((CurrRes$TotalNumCombosUsed / CurrRes$MaxPossibleCombos)[1], ndig))
        } # end k loop
        
        #if(change==1)
        #if(change < 0
        if(SuppPlots==TRUE & change==1 & plotType==2){
          legendloc <- "topleft"
          #legendloc <- "bottomright"
        }else if(plotType==11){
          legendloc <- "topleft"
        }else if(change==1 & plotType==2 & type=="LowerThirdOverlap"){
          legendloc <- "bottomright"
        }else if(SuppPlots==FALSE & change==1 & plotType==2){
          legendloc <- "topleft"
        }else if(SuppPlots==FALSE & change==2 & plotType==10 & type=="genesMeanGeneMaxVarilrScaleGr90" & TwentySamplesOnly==TRUE){
          legendloc <- "topleft"
        }else{
          legendloc <- "bottomright"
        }
        #legendloc=="top" | legendloc=="topleft")
        if(SuppPlots==T & change==1 & plotType==2){
          Legendpt.cexToUse <- 1.25
          LegendcexToUse <- 1.25
        }else if (SuppPlots==F & change==1 & plotType==2){
          Legendpt.cexToUse <- 1.75
          LegendcexToUse <- 1.75
        }else if(SuppPlots==F & change==2 & plotType==11) {
          Legendpt.cexToUse <- 1.75
          LegendcexToUse <- 1.75
        }else{
          Legendpt.cexToUse <- 2
          LegendcexToUse <- 2
        }
        
        #stdal <- F
        
        curr_save_dir <- paste0(def_wd2, dirr)
        if(!dir.exists(curr_save_dir)){dir.create(curr_save_dir, recursive = TRUE)}
        
        fil_piece <- paste0("ROCChange", change, type, "PlotNumber", plotType)
        fil_name <- paste0(curr_save_dir, fil_piece , ".tex")
        
        if(change==1){                                                                                  #height=4,width=6.5
          if(SuppPlots==TRUE){
            tikz(file = fil_name, height=6.5,width=6.5, standAlone = stdal, sanitize = F)
          }else{
            tikz(file = fil_name, height=6.5,width=6.5, standAlone = stdal, sanitize = F)
          }
        }else{
          if(SuppPlots==TRUE){
            tikz(file = fil_name, height=6.5,width=6.5, standAlone = stdal, sanitize = F)
          }else{
            tikz(file = fil_name, height=6.5,width=6.5, standAlone = stdal, sanitize = F)
          }
        }
        
        for (k in 1:length(methdstouse)){
          curr_methd <- methdstouse[k]
          #CurrRes <- get(paste0(curr_methd, "Res", change, type))
          CurrResT <- get(paste0(curr_methd, "Res", change, type))
          CurrRes <- CurrResT[order(CurrResT$pvalcutoff),]
          #if(change==1)
          if(change < 0){
            # plot(CurrRes$pvalcutoff, CurrRes$fpr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
            #      col = colrs[k], pch = 19, type = "l", lty=lins[k],  lwd = 4)
            plot(CurrRes$pvalcutoff, CurrRes$fpr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                 col = colrs[k], pch = 19, type = "l", lty=lins[k],  lwd = 4)
            abline(0, 1)
          }else{
            #Increase plot margins to ensure text does not get cut off
            par(mar=c(5.1, 4.1, 4.1, 2.1)+c(0,.8,.8,0))
            
            if(plotType==1){
              pvalcutoff_0.01 <- CurrRes[CurrRes$pvalcutoff==0.01,]
              if(nrow(pvalcutoff_0.01)==0){stop()}
              pvalcutoff_0.05 <- CurrRes[CurrRes$pvalcutoff==0.05,]
              pvalcutoff_0.10 <- CurrRes[CurrRes$pvalcutoff==0.10,]
              
              plot(CurrRes$fpr, CurrRes$sens, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                   col = colrs[k], pch = 19, type = "l", lty=lins[k],  lwd = 4)
              # plot(CurrRes$fpr, CurrRes$sens, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
              #      col = colrs[k], pch = 19, type = "s", lty=lins[k],  lwd = 4)
              par(new = TRUE)
              plot(pvalcutoff_0.01$fpr_padjfdr, pvalcutoff_0.01$sens_padjfdr, col = colrs[k], pch = 17, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE, cex = 1.5)
              par(new = TRUE)
              plot(pvalcutoff_0.05$fpr_padjfdr, pvalcutoff_0.05$sens_padjfdr, col = colrs[k], pch = 18, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE, cex = 2)
              par(new = TRUE)
              plot(pvalcutoff_0.10$fpr_padjfdr, pvalcutoff_0.10$sens_padjfdr, col = colrs[k], pch = 15, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE, cex = 1.5)
              
              
            }else if(plotType==2){
              if(SuppPlots==TRUE){
                CurrRes <- CurrRes[CurrRes$pvalcutoff>0.0001 & CurrRes$pvalcutoff < xlim1[2],]
              }else{
                #CurrRes <- CurrRes[CurrRes$pvalcutoff>0.0001 & CurrRes$pvalcutoff < 0.20,]
                CurrRes <- CurrRes[CurrRes$pvalcutoff>0.0001 & CurrRes$pvalcutoff < xlim1[2],]
              }
              
              if(curr_methd %in% RATsmethds){
                CurrRes2 <- subset(CurrRes, CurrRes$pvalcutoff==0.05)
                plot(CurrRes2$pvalcutoff, CurrRes2$fpr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                     col = colrs[k], pch = 19, lty=lins[k], cex = 2)
              }else{
                plot(CurrRes$pvalcutoff, CurrRes$fpr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                     col = colrs[k], pch = 19, type = "l", lty=lins[k],  lwd = 4)
              }
              
              #abline(0, 1)
              abline(0, 1, lty=2, col = 'gray')
            }else if(plotType==3){
              plot(CurrRes$pvalcutoff, CurrRes$sens, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                   col = colrs[k], pch = 19, type = "l", lty=lins[k],  lwd = 4)
            }else if(plotType==4){
              #FDR vs Sens (using fdr adjusted pvalues for both)
              #Restrict for the plot of fdr vs sens to be only pvalue cutoffs that are reasonably large
              #This is because fdr can be relatively high with a low pvalue cutoff if there are a low number
              #of positives and a lot of them happen to reject- these cutoffs aren't really realistic so don't consider them
              pvalcutoff_0.01 <- CurrRes[CurrRes$pvalcutoff==0.01,]
              pvalcutoff_0.05 <- CurrRes[CurrRes$pvalcutoff==0.05,]
              pvalcutoff_0.10 <- CurrRes[CurrRes$pvalcutoff==0.10,]
              
              CurrRes2 <- subset(CurrRes, CurrRes$pvalcutoff > 1e-4)
              
              
              if(curr_methd %in% RATsmethds){
                plot(CurrRes2$fdr_padjfdr, CurrRes2$sens_padjfdr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                     col = colrs[k], pch = RATspch, cex = RATscex, lty=lins[k])
                #par(new = TRUE)
                #plot(pvalcutoff_0.05$fdr_padjfdr, pvalcutoff_0.05$sens_padjfdr, col = colrs[k], pch = 17, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE)
              }else{
                plot(CurrRes2$fdr_padjfdr, CurrRes2$sens_padjfdr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                     col = colrs[k], pch = 19, type = "l", lty=lins[k],  lwd = 4)
                par(new = TRUE)
                #triangles, diamonds, squares are plotting the points where eFDR values are 0.01, 0.05, and 0.10
                plot(pvalcutoff_0.01$fdr_padjfdr, pvalcutoff_0.01$sens_padjfdr, col = colrs[k], pch = 17, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE, cex = 1.5)
                par(new = TRUE)
                plot(pvalcutoff_0.05$fdr_padjfdr, pvalcutoff_0.05$sens_padjfdr, col = colrs[k], pch = 18, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE, cex = 2)
                par(new = TRUE)
                plot(pvalcutoff_0.10$fdr_padjfdr, pvalcutoff_0.10$sens_padjfdr, col = colrs[k], pch = 15, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE, cex = 1.5)
              }
              
            }else if(plotType==5){
              #Restrict for the plot of fdr vs sens to be only pvalue cutoffs that are reasonably large
              #This is because fdr can be relatively high with a low pvalue cutoff if there are a low number
              #of positives and a lot of them happen to reject- these cutoffs aren't really realistic so don't consider them
              CurrRes2 <- subset(CurrRes, CurrRes$pvalcutoff > 1e-5)
              
              plot(CurrRes2$pvalcutoff, CurrRes2$fdr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                   col = colrs[k], pch = 19, type = "l", lty=lins[k],  lwd = 4)
              
              #abline(0, 1)
            }else if(plotType==6){
              plot(CurrRes$pvalcutoff, CurrRes$acc, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                   col = colrs[k], pch = 19, type = "l", lty=lins[k],  lwd = 4)
            }else if(plotType==7){
              plot(CurrRes$fpr, CurrRes$acc, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                   col = colrs[k], pch = 19, type = "l", lty=lins[k],  lwd = 4)
            }else if(plotType==8){ #eFDR vs accuracy
              
              pvalcutoff_0.05 <- CurrRes[CurrRes$pvalcutoff==0.05,]
              CurrRes2 <- subset(CurrRes, CurrRes$pvalcutoff > 1e-4)
              t1 <- abs(CurrRes2$fdr_padjfdr - 0.01)
              fdr_0.01 <- CurrRes2[which(t1==min(t1)),]
              if(nrow(fdr_0.01) > 1){
                fdr_0.01 <- fdr_0.01[1,]
              }
              
              t2 <- abs(CurrRes2$fdr_padjfdr - 0.05)
              fdr_0.05 <- CurrRes2[which(t2==min(t2)),]
              if(nrow(fdr_0.05) > 1){
                fdr_0.05 <- fdr_0.05[1,]
              }
              
              
              t3 <- abs(CurrRes2$fdr_padjfdr - 0.10)
              fdr_0.10 <- CurrRes2[which(t3==min(t3)),]
              if(nrow(fdr_0.10) > 1){
                fdr_0.10 <- fdr_0.10[1,]
              }
              
              # fdr_0.01 <- CurrRes[420,]
              # fdr_0.05 <- CurrRes[CurrRes$pvalcutoff==0.05,]
              # fdr_0.10 <- CurrRes[CurrRes$pvalcutoff==0.10,]
              
              if(curr_methd %in% RATsmethds){
                CurrRes3 <- subset(CurrRes2, CurrRes2$pvalcutoff==0.05)
                plot(CurrRes3$pvalcutoff, CurrRes3$acc_padjfdr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                     col = colrs[k], pch = RATspch, cex = RATscex, lty=lins[k])
                #par(new = TRUE)
                #plot(fdr_0.05$pvalcutoff, fdr_0.05$acc_padjfdr, col = colrs[k], pch = 18, cex = 1.5, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE)
                
              }else{
                plot(CurrRes2$pvalcutoff, CurrRes2$acc_padjfdr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                     col = colrs[k], pch = 19, type = "l", lty=lins[k],  lwd = 4)
                #triangles, diamonds, squares are plotting points where true FDR values are 0.01, 0.05, and 0.10
                par(new = TRUE)
                plot(fdr_0.01$pvalcutoff, fdr_0.01$acc_padjfdr, col = colrs[k], pch = 17, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE, cex = 1.5)
                par(new = TRUE)
                plot(fdr_0.05$pvalcutoff, fdr_0.05$acc_padjfdr, col = colrs[k], pch = 18, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE, cex = 2)
                par(new = TRUE)
                plot(fdr_0.10$pvalcutoff, fdr_0.10$acc_padjfdr, col = colrs[k], pch = 15, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE, cex = 1.5)
                
              }
              
              
              eFDRvals[[k]] <- fr(pvalcutoff_0.05$fdr_padjfdr, 2)
            }else if(plotType==9){ #eFDR vs sens based on fdr adjusted pvalues
              pvalcutoff_0.01 <- CurrRes[CurrRes$pvalcutoff==0.01,]
              pvalcutoff_0.05 <- CurrRes[CurrRes$pvalcutoff==0.05,]
              pvalcutoff_0.10 <- CurrRes[CurrRes$pvalcutoff==0.10,]
              
              CurrRes2 <- subset(CurrRes, CurrRes$pvalcutoff > 5e-3)
              t1 <- abs(CurrRes2$fdr_padjfdr - 0.01)
              fdr_0.01 <- CurrRes2[which(t1==min(t1)),]
              if(nrow(fdr_0.01) > 1){
                fdr_0.01 <- fdr_0.01[1,]
              }
              
              t2 <- abs(CurrRes2$fdr_padjfdr - 0.05)
              fdr_0.05 <- CurrRes2[which(t2==min(t2)),]
              if(nrow(fdr_0.05) > 1){
                fdr_0.05 <- fdr_0.05[1,]
              }
              
              
              t3 <- abs(CurrRes2$fdr_padjfdr - 0.10)
              fdr_0.10 <- CurrRes2[which(t3==min(t3)),]
              if(nrow(fdr_0.10) > 1){
                fdr_0.10 <- fdr_0.10[1,]
              }
              
              #print(paste0("Pval Cutoff that gives a true FDR (based on adjusted pvalues) of 0.01 is ", fdr_0.01$pvalcutoff, " for ", curr_methd))
              #print(paste0("Pval Cutoff that gives a true FDR (based on adjusted pvalues) of 0.05 is ", fdr_0.05$pvalcutoff, " for ", curr_methd))
              #print(paste0("Pval Cutoff that gives a true FDR (based on adjusted pvalues) of 0.10 is ", fdr_0.10$pvalcutoff, " for ", curr_methd))
              if(curr_methd %in% RATsmethds){
                CurrRes2 <- subset(CurrRes, CurrRes$pvalcutoff==0.05)
                plot(CurrRes2$pvalcutoff, CurrRes2$sens_padjfdr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                     col = colrs[k], pch = RATspch, cex = RATscex, lty=lins[k])
                #par(new = TRUE)
                #plot(pvalcutoff_0.05$pvalcutoff, pvalcutoff_0.05$sens_padjfdr, col = colrs[k], pch = 17, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE)
              }else{
                CurrRes2 <- subset(CurrRes, CurrRes$pvalcutoff > 1e-6)
                plot(CurrRes2$pvalcutoff, CurrRes2$sens_padjfdr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                     col = colrs[k], pch = 19, type = "l", lty=lins[k],  lwd = 4)
                #triangles, diamonds, squares are plotting points where true FDR values are 0.01, 0.05, and 0.10
                par(new = TRUE)
                plot(fdr_0.01$pvalcutoff, fdr_0.01$sens_padjfdr, col = colrs[k], pch = 17, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE, cex = 1.5)
                par(new = TRUE)
                plot(fdr_0.05$pvalcutoff, fdr_0.05$sens_padjfdr, col = colrs[k], pch = 18, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE, cex = 2)
                par(new = TRUE)
                plot(fdr_0.10$pvalcutoff, fdr_0.10$sens_padjfdr, col = colrs[k], pch = 15, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE, cex = 1.5)
              }
              
              #plot(CurrRes2$pvalcutoff, CurrRes2$sens_padjfdr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
              #     col = colrs[k], pch = 19, type = "l", lty=lins[k],  lwd = 4)
              
              eFDRvals[[k]] <- fr(pvalcutoff_0.05$fdr_padjfdr, 2)
            }else if(plotType==10){ #FPR(based on fdr adjusted pvalues) vs sens (also based on FDR adjusted pvalues)
              #Find points at which the eFDR is 0.01, 0.05, or 0.10
              #Note that this is eFDR even though pvalcutoff is not "adjusted" because the plotted points are based on
              #eFDR adjusted pvalues
              pvalcutoff_0.01 <- CurrRes[CurrRes$pvalcutoff==0.01,]
              if(nrow(pvalcutoff_0.01)==0){stop()}
              pvalcutoff_0.05 <- CurrRes[CurrRes$pvalcutoff==0.05,]
              pvalcutoff_0.10 <- CurrRes[CurrRes$pvalcutoff==0.10,]
              
              CurrRes2 <- subset(CurrRes, CurrRes$pvalcutoff > 1e-6)
              #CurrRes2 <- CurrRes2[order(CurrRes2$pvalcutoff),]
              
              if(curr_methd %in% RATsmethds){
                plot(CurrRes2$fpr_padjfdr, CurrRes2$sens_padjfdr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                     col = colrs[k], pch = RATspch, cex = RATscex, lty=lins[k])
              }else{
                plot(CurrRes2$fpr_padjfdr, CurrRes2$sens_padjfdr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                     col = colrs[k], pch = 19, type = "l", lty=lins[k],  lwd = 4)
                #par(new = TRUE)
                #plot(pvalcutoff_0.05$fpr_padjfdr, pvalcutoff_0.05$sens_padjfdr, col = colrs[k], pch = 17, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE)
                #triangles, diamonds, squares are plotting points where eFDR values are 0.01, 0.05, and 0.10
                par(new = TRUE)
                plot(pvalcutoff_0.01$fpr_padjfdr, pvalcutoff_0.01$sens_padjfdr, col = colrs[k], pch = 17, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE, cex = 2)
                par(new = TRUE)
                plot(pvalcutoff_0.05$fpr_padjfdr, pvalcutoff_0.05$sens_padjfdr, col = colrs[k], pch = 18, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE, cex = 2.5)
                par(new = TRUE)
                plot(pvalcutoff_0.10$fpr_padjfdr, pvalcutoff_0.10$sens_padjfdr, col = colrs[k], pch = 15, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE, cex = 2)
              }
              
              
            }else if(plotType==11){
              #CurrRes2 <- CurrRes[CurrRes$pvalcutoff>0.00075,]
              #CurrRes2 <- CurrRes[CurrRes$pvalcutoff>0.001 & CurrRes$pvalcutoff < 0.20,]
              CurrRes2 <- CurrRes[CurrRes$pvalcutoff>0.005 & CurrRes$pvalcutoff < xlim1[2],]
              if(curr_methd %in% RATsmethds){
                CurrRes2 <- subset(CurrRes2, CurrRes2$pvalcutoff==0.05)
                plot(CurrRes2$pvalcutoff, CurrRes2$fdr_padjfdr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                     col = colrs[k], pch = 19, lty=lins[k], cex = 2)
              }else{
                plot(CurrRes2$pvalcutoff, CurrRes2$fdr_padjfdr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                     col = colrs[k], pch = 19, type = "l", lty=lins[k],  lwd = 4)
              }
              #abline(0, 1)
              abline(0, 1, lty=2, col = 'gray')
            }else if(plotType==12){
              CurrRes <- subset(CurrRes, CurrRes$pvalcutoff>0.001 & CurrRes$pvalcutoff < 0.50)
              if(curr_methd %in% RATsmethds){
                CurrRes2 <- subset(CurrRes, CurrRes$pvalcutoff==0.05)
                plot(CurrRes2$pvalcutoff, CurrRes2$fpr_padjfdr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                     col = colrs[k], pch = RATspch, cex = RATscex, lty=lins[k])
              }else{
                plot(CurrRes$pvalcutoff, CurrRes$fpr_padjfdr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                     col = colrs[k], pch = 19, type = "l", lty=lins[k],  lwd = 4)
              }
              
              #abline(0, 1)
              abline(0, 1, lty=2, col = 'gray')
            }else if(plotType==13){
              CurrRes <- subset(CurrRes, CurrRes$pvalcutoff>0.001 & CurrRes$pvalcutoff < 0.50)
              if(curr_methd %in% RATsmethds){
                CurrRes2 <- subset(CurrRes, CurrRes$pvalcutoff==0.05)
                plot(CurrRes2$pvalcutoff, 1-CurrRes2$fpr_padjfdr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                     col = colrs[k], pch = RATspch, cex = RATscex, lty=lins[k])
              }else{
                plot(CurrRes$pvalcutoff, 1-CurrRes$fpr_padjfdr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                     col = colrs[k], pch = 19, type = "l", lty=lins[k],  lwd = 4)
              }
              
            }else if(plotType==14){
              #CurrRes <- subset(CurrRes, CurrRes$pvalcutoff>0.001 & CurrRes$pvalcutoff < 0.50)
              CurrRes2 <- subset(CurrRes, CurrRes$pvalcutoff > 5e-3)
              t1 <- abs(CurrRes2$pvalcutoff - 0.01)
              efdr_0.01 <- CurrRes2[which(t1==min(t1)),]
              if(nrow(efdr_0.01) > 1){
                efdr_0.01 <- efdr_0.01[1,]
              }
              
              t2 <- abs(CurrRes2$pvalcutoff - 0.05)
              efdr_0.05 <- CurrRes2[which(t2==min(t2)),]
              if(nrow(efdr_0.05) > 1){
                efdr_0.05 <- efdr_0.05[1,]
              }
              
              
              t3 <- abs(CurrRes2$pvalcutoff - 0.10)
              efdr_0.10 <- CurrRes2[which(t3==min(t3)),]
              if(nrow(efdr_0.10) > 1){
                efdr_0.10 <- efdr_0.10[1,]
              }
              
              if(curr_methd %in% RATsmethds){
                CurrRes2 <- subset(CurrRes, CurrRes$pvalcutoff==0.05)
                plot(CurrRes2$fdr_padjfdr, CurrRes2$sens_padjfdr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                     col = colrs[k], pch = RATspch, cex = RATscex, lty=lins[k])
              }else{
                plot(CurrRes$fdr_padjfdr, CurrRes$sens_padjfdr, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                     col = colrs[k], pch = 19, type = "l", lty=lins[k],  lwd = 4)
                
                #triangles, diamonds, squares are plotting points where eFDR values (ie adjusted pvalue cutoff values) are 0.01, 0.05, and 0.10
                par(new = TRUE)
                plot(efdr_0.01$fdr_padjfdr, efdr_0.01$sens_padjfdr, col = colrs[k], pch = 17, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE, cex = 1.5)
                par(new = TRUE)
                plot(efdr_0.05$fdr_padjfdr, efdr_0.05$sens_padjfdr, col = colrs[k], pch = 18, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE, cex = 2)
                par(new = TRUE)
                plot(efdr_0.10$fdr_padjfdr, efdr_0.10$sens_padjfdr, col = colrs[k], pch = 15, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE, cex = 1.5)
                
              }
              
            }else if(plotType==15){
              real_sens <- CurrRes$sens_padjfdr_real
              for (p in 2:length(real_sens)){
                if(real_sens[p] < real_sens[p-1]){
                  real_sens[p] <- real_sens[p-1]
                }
              }          
              plot(CurrRes$pvalcutoff, real_sens, xlim = xlim1, ylim = ylim1, ann=FALSE, axes=FALSE,
                   col = colrs[k], pch = 19, type = "l", lty=lins[k],  lwd = 4)
            }
            
            
          }
          
          if(k!=length(methdstouse)){
            par(new=TRUE)
          }
          
        } #end k loop
        for (k in 1:length(methdstouse)){
          curr_methd <- methdstouse[k]
          curr_methd_name <- methdnamestouse[k]
          curr_used <- get(paste0("p", curr_methd, "used"))
          #plotType==8 | plotType==9
          if(plotType==3.14159){
            val1 <- eFDRvals[[k]]
            assign(paste0("legendnum", k), paste0(curr_methd_name, " (", curr_used, ") ", "(trueFDR is " , val1, " at eFDR of 0.05", ")"))
          }else{
            #assign(paste0("legendnum", k), paste0(curr_methd_name, " (", curr_used, ")"))
            assign(paste0("legendnum", k), paste0(curr_methd_name))
          }
          
        }
        # l1 <- paste0("Compositional (",pcompused, ")")
        # l2 <- paste0("DRIMSeq (", pdrimused, ")")
        # l3 <- paste0("CompositionalMI (", pcompmiused, ")")
        # l4 <- paste0("Comp. No Other Groups (",pcompnootherused, ")")
        if(SuppPlots==T & change==1 & plotType==2){
          legend(legendloc, legend = mget(ls(pattern="legendnum")), col = colrs, pch = 19, pt.cex = Legendpt.cexToUse, cex = LegendcexToUse)
        }else{
          legend(legendloc, legend = mget(ls(pattern="legendnum")), col = colrs, pch = 19, pt.cex = Legendpt.cexToUse, cex = LegendcexToUse)
        }
        
        
        
        
        box()
        axis(side = 1, at = seq(0, upperx, 0.05), cex.axis=2)
        #if(change < 0)
        if(change==1){
          axis(side = 2, at = seq(0, upperx, 0.05), cex.axis=2)
        }else{
          if(plotType==1 | plotType==3 | plotType==4 | plotType==6 | plotType==7 | plotType==8 | plotType==9 | plotType==10 | plotType==13 | plotType==14 | plotType==15){
            #axis(side = 2, at = c(0.50, 0.60, 0.70, 0.80, 0.90, 1), labels = c("0.50", "0.60", "0.70", "0.80", "0.90", "1"))
            axis(side = 2, at = c(0, 0.20, 0.40, 0.60, 0.80, 1), labels = c("0", "0.20", "0.40", "0.60", "0.80", "1"), cex.axis=2)
            #axis(side = 2, at = c(0, 0.05, 0.10, 0.15, 0.20, 0.25), labels = c("0", "0.05", "0.10", "0.15", "0.20", "0.25"))
          }else{
            axis(side = 2, at = seq(0, upperx, 0.05), cex.axis=2)
          }
        }
        #if(change==1)
        if(change < 0){
          mtext("P Value Cutoff", side = 1, line = 2, cex = 2)
          mtext('FPR', side = 2, line = 3, cex = 2)
        }else{
          if(plotType==1){
            mtext('FPR', side = 1, line = 3, cex = 2)
            mtext("Sensitivity", side = 2, line = 3, cex = 2)
          }else if(plotType==2){
            mtext("$p$-value Significance Threshold", side = 1, line = 3, cex = 2)
            mtext('FPR', side = 2, line = 3, cex = 2)
          }else if (plotType==3){
            mtext("$p$-value Significance Threshold", side = 1, line = 3, cex = 2)
            mtext("Sensitivity", side = 2, line = 3, cex = 2)
          }else if(plotType==4){
            mtext('True FDR', side = 1, line = 3, cex = 2)
            mtext("Sensitivity", side = 2, line = 3, cex = 2)
          }else if(plotType==5){
            mtext("$p$-value Significance Threshold", side = 1, line = 3, cex = 2)
            mtext('FDR', side = 2, line = 3, cex = 2)
          }else if(plotType==6){
            mtext("$p$-value Significance Threshold", side = 1, line = 3, cex = 2)
            mtext('Accuracy', side = 2, line = 3, cex = 2)
          }else if(plotType==7){
            mtext("FPR", side = 1, line = 3, cex = 2)
            mtext('Accuracy', side = 2, line = 3, cex = 2)
          }else if(plotType==8){
            mtext("eFDR", side = 1, line = 3, cex = 2)
            mtext('Accuracy', side = 2, line = 3, cex = 2)
          }else if(plotType==9){
            mtext("eFDR", side = 1, line = 3, cex = 2)
            mtext('Sensitivity', side = 2, line = 3, cex = 2)
          }else if(plotType==10){
            #mtext("FPR (Based on FDR Adj Pvalues)", side = 1, line = 2)
            mtext("FPR", side = 1, line = 3, cex = 2)
            mtext('Sensitivity', side = 2, line = 3, cex = 2)
          }else if(plotType==11){
            mtext("eFDR", side = 1, line = 3, cex = 2)
            mtext('True FDR', side = 2, line = 3, cex = 2)
          }else if(plotType==12){
            mtext("eFDR", side = 1, line = 3, cex = 2)
            mtext('FPR', side = 2, line = 3, cex = 2)
          }else if(plotType==13){
            mtext("eFDR", side = 1, line = 3, cex = 2)
            mtext('True Negative Rate (1-FPR)', side = 2, line = 3, cex = 2)
          }else if(plotType==14){
            mtext("True FDR", side = 1, line = 3, cex = 2)
            mtext('Sensitivity', side = 2, line = 3, cex = 2)
          }else if(plotType==15){
            mtext("eFDR", side = 1, line = 3, cex = 2)
            mtext('True Sensitivity', side = 2, line = 3, cex = 2)
          }
        }
        
        
        if(type=="Full"){
          mess <- "All Genes"
        }else if(type=="genesMeanGeneMaxVarilrScaleGr90"){
          mess <- paste0("Genes in the Top 10\\% of \nInferential Variability")
        }else if(type=="genesMeanGeneMaxVarilrScaleGr80"){
          mess <- paste0("Genes in the Top 20\\% of \nInferential Variability")
        }else if(type=="genesMeanGeneMaxVarilrScaleHigh"){
          mess <- paste0("Genes in the Top 33.3\\% of \nInferential Variability")
        }else if(type=="genesMeanGeneMaxVarilrScaleLow"){
          mess <- paste0("Genes in the Lower 33.3\\% of \nInferential Variability")
        }else if(type=="genesMeanGeneMaxVarilrScaleMedium"){
          mess <- paste0("Genes in the Middle 33.3\\% of \nInferential Variability")
        }else if(type=="genesMeanGeneMaxInfRVTPMGrEighty"){
          #mess <- paste0("Genes in the Upper 20\\% of Maximum \nInfRV Level Across its Transcripts")
          mess <- paste0("Genes in the Top 20\\% of \nInferential Variability")
        }else if(type=="genesMeanGeneMaxInfRVTPMGrNinety"){
          #mess <- paste0("Genes in the Upper 10\\% of Maximum \nInfRV Level Across its Transcripts")
          mess <- paste0("Genes in the Top 10\\% of \nInferential Variability")
        }else if(type=="LowerThirdOverlap"){
          mess <- paste0("Genes in the Lowest Tertile of Overlap \nBetween the Sequences of its Transcripts")
        }else if(type=="HighThirdOverlap"){
          mess <- paste0("Genes in the Highest Tertile of Overlap \nBetween the Sequences of its Transcripts")
        }else{
          mess <- paste0("Genes with ", type)
        }
        #tempp <- get(paste0("genes", type))
        tempp <- NULL
        mess2 <- NULL
        if(!is.null(tempp)){
          propused <- length(tempp)/ngenestotal
          propused2 <- fr(propused*100, 0)
          #mess2 <- paste0(" (", length(tempp), " Genes Used Out of ", ngenestotal, ")")
          mess2 <- paste0(" (", propused2, "\\%", " of Genes Used)")
        }
        if(is.null(mess2)){
          #if(plotType==1 & change!=1)
          if(plotType==1){
            mtext(paste0('ROC Curve For Change Value ', change, " for ", mess), side = 3)
          }else{
            #mtext(paste0('Curve For Change Value ', change, " for ", mess), side = 3)
            mtext(mess, side = 3, cex = 2, line = 1.25)
          }
        }else{
          #if(plotType==1 & change!=1)
          if(plotType==1){
            mtext(paste0('ROC Curve For Change Value ', change, " for ", mess, mess2), side = 3)
          }else{
            #mtext(paste0('Curve For Change Value ', change, " for ", mess, mess2), side = 3)
            mtext(paste0('Change Value ', change, " for ", mess, mess2), side = 3)
            #mtext(paste0(mess, mess2), side = 3)
          }
        }
        
        rm(mess2)
        rm(propused)
        dev.off()
        
        system(paste0("pdflatex ","-output-directory=", curr_save_dir, " ", fil_name), ignore.stdout = TRUE)
        system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".aux"))
        system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".log"))
        system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".tex"))
        rm(fil_name)
        rm(fil_piece)
        
      }
    }
  }
}






