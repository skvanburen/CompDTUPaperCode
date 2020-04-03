#Code to calculate the amount of "overlap" of transcripts of a specific gene
#This is defined as the ratio of number of exon base pairs that are on multiple transcripts 
#within a specific gene to the ones that are not on multiple transcripts
#Procedure to calculate this:
#1. Create the Transcript database from the .gtf GENCODE annotation (using v27 currently)
#2. Extract the exon ranges (no introns) for each transcript
#3. Calculate the total length of exon base pairs in a specific transcript that are not found in any other transcript and divide that number by the total length of of all exon sequences within that transcript. This will give each gene an overlap value between 0 and 1, with 1 being full/perfect overlap and 0 being no overlap at all.


library(parallel)
onCluster <- TRUE #Change to True if running on the cluster, False if running on Macbook

#Set FilterGenes to FALSE to calculate the overlap for all genes
FilterGenes <- FALSE
#FilterCutoff <- 0.50

if(onCluster==TRUE){
  def_wd <- "~/res/GEUV1Data/"
  setwd(def_wd)
  clust <- makeCluster(4) # Change this to match whatever number of cores requested on cluster or are on macbook
  func_loc <-  "~/code/CompFunctions.R"
  #load(paste0(def_wd,"Salmon/SalmonData.RData"))
  source("~/.Rprofile")
}else{
  def_wd <- "/Users/Scott/Documents/Dissertation Data/SQCCData/"
  setwd(def_wd)
  clust <- makeCluster(4) # Change this to match whatever number of cores requested on cluster or are on macbook
  #load("SalmonData.RData")
  source("/Users/Scott/Documents/Dissertation/res/SQCCData/SQCCFunctions.R")
}
source(func_loc)
library(GenomicFeatures)
library(Biostrings)
#library(BSgenome.Hsapiens.NCBI.GRCh38)

#Load the necessary packages to the cluster since they are not in the default location
if(onCluster==TRUE){
  clusterEvalQ(clust, source("~/.Rprofile"))
}
#clusterEvalQ(clust, library(BSgenome.Hsapiens.NCBI.GRCh38))
#clusterEvalQ(clust, source(func_loc))
clusterEvalQ(clust, library(GenomicFeatures))
clusterEvalQ(clust, library(Biostrings))

#Use gene names from abDatasets to get the set of genes that is used in the analysis to avoid unnecessary computation
load("abGene.RData")

if(FilterGenes==TRUE){
  #This will give a vector 
  genestouset <- unique(abGene$gene_id[abGene$MeanTGE > FilterCutoff & abGene$NTrans > 1])
  genestouse <- genestouset[genestouset %in% genenames]
}else{
  genestouse <- genenames
}


#Matches transcripts to their associated gene
load("tx2gene.RData")




#Rename the seqnames to add chr in front of 1, 2, 3, etc so the BSgenome names match the gencode gtf annotation
#Check on this, not sure why this is necessary since both and BSgenome and annotation are NCBI GRCh38 and not UCSC
#But, as of now it seems necessary
# seqnames(Hsapiens)[seqnames(Hsapiens) %in% as.character(1:22)] <- paste0("chr", seqnames(Hsapiens)[1:22])
# seqnames(Hsapiens)[seqnames(Hsapiens)=="MT"] <- "chrM"
# seqnames(Hsapiens)[seqnames(Hsapiens)=="X"] <- "chrX"
# seqnames(Hsapiens)[seqnames(Hsapiens)=="Y"] <- "chrY"


#Create transcript database from the gencode v27 annotation directly
txdb <- makeTxDbFromGFF("~/gencode.v27.annotation.gtf.gz")

#Extract list of exons by transcript
exon <- exonsBy(txdb, by = "tx", use.names=TRUE)
#exons <- exons(txdb, columns = c("exon_id", "exon_name", "tx_id", "tx_name", "gene_id"))

exon_ranges <- ranges(exon)

#Extract the transcript sequences corresponding to the exons only
#trans_seq_exon_only <- extractTranscriptSeqs(Hsapiens, exon)
start <- proc.time()
overlapres <- data.frame(t(data.frame(parLapply(clust, genestouse, calcOverlap2, tx2gene = tx2gene, exon_ranges = exon_ranges))))
colnames(overlapres) <- "overlap"
rownames(overlapres) <- genestouse
overlaprestime <- proc.time() - start


save(overlapres, overlaprestime, file = "OverlapRes.RData")


#Code below is older, probably can delete this soon



# for(i in 1:length(names(curr_exon_ranges))){
#   curr_tns <- names(curr_exon_ranges)[1]
#   other_tns <- names(curr_exon_ranges)[-1]
#   
#   uniq[[i]] <- setdiff(curr_exon_ranges[curr_tns], curr_exon_ranges[other_tns])
# }
# 
# 
# x <- "ENSG00000000003.14" #sample, just for testing
# 
# 
# 
# 
# #Sequence of the exons from each transcript
# #Just given to be able to compare the exons of a specific transcript to the final exon-only transcript sequences
#   #Which should be all the exon sequences pasted together
# exon_seqs <- getSeq(Hsapiens, exon)
# curr_trans <- exon_seqs[[1]]
# 
# 
# #This is the protein coding sequence, not what is of interest now
# #cds_seqs <- extractTranscriptSeqs(Hsapiens, cdsBy(txdb, by="tx", use.names=TRUE)
# 
# 
# 
# 
# 
# #Now, calculate the overlap among exon regions of the transcripts within a specific gene
#   #Where the overlap score is discussed at the top of the file
# calcOverlap <- function(x, tx2gene, trans_seq_exon_only, testingCode = FALSE){
#   #x <- "ENSG00000000003.14" #sample, just for testing
#   
#   if(testingCode==TRUE){
#     #Test list
#     curr_trans <- trans_seq_exon_only
#     trans_widths <- nchar(curr_trans)
#     
#   }else{
#     #Restrict analysis to a specific gene and its transcripts
#     sub <- subset(tx2gene, tx2gene$gene_id==x)
#     
#     #Get sequences only from the transcripts associated with the current gene
#     curr_trans <- trans_seq_exon_only[names(trans_seq_exon_only) %in% sub$tx_id]
#     trans_widths <- width(curr_trans)
#     
#   }
#   
#   
#   maxwidth <- max(trans_widths)
#   sum_lengths <- sum(trans_widths)
#   
#   ntrans <- length(curr_trans)
# 
#   
#   #This has been changed for now #Start scores at 1 and change any invalid positions to 0 after the code below is done
#   scores <- matrix(0, nrow = ntrans, ncol = maxwidth)
#   mainseq <- 1:ntrans
#   for(i in mainseq){
#     for(j in 1:maxwidth){
#       if(j > trans_widths[i]){
#         next
#       }else{
#         for (k in mainseq[!(mainseq %in% i)]){
#           if(j > trans_widths[k]){
#             next
#           }else{
#             if(subseq(curr_trans[[i]], j, j)==subseq(curr_trans[[k]], j, j)){
#               scores[i, j] <- scores[i, j] + 1
#             }
#           }
#         }
#       }
#     }
#   }
#   
#   #Now, change any scores corresponding to a spot that is beyond the width of that transcript to 0
#   for(i in mainseq){
#     for(j in 1:maxwidth){
#       if(j > trans_widths[i]){
#         scores[i,j] <- 0
#       }
#     }
#   }
#   
#   
#   #Now, get final overlap measure
#   
#   
#   #Sum the total score matrix and divide by the "maximum possible score", which is the ntrans * sum of trans lengths
#   #As of now maxscore could never be attainted unless all transcripts are the same length?  Is this good?- I'd say yes
#   #Maybe maxscore should be defined differently?
#   #maxscore <- (ntrans*sum_lengths) - (ntrans*maxwidth)
#   maxscore <- ((ntrans-1)*sum_lengths)
#   
#   #Score is a normalized measure between 0 and 1 that will be 0 with no overlap and 1 with perfect overlap
#   #gene_score <- (sum(scores)-(ntrans*maxwidth))/maxscore
#   gene_score <- (sum(scores))/maxscore
#   return(list(gene_score = gene_score, scores = scores))
# }
# 
# 
# 
# 
#   #First, restrict to only genes used in the compositional analysis
# genestouse <- genenames
# 
# #Test
# tr <- list("GTC","GA" ,"GA" ,"AG")
# tr <- list("GG", "GG", "GGG")
# tr <- list("AG","TA" ,"CT" ,"GC")
# tr <- list("AG","AG" ,"AG" ,"AG")
# tr <- list("AG", "TCT")
# tr <- list("AG","AG" ,"AG" ,"AGGTGTGTGTGTGTGTGTGTGTGTGTGTGTG")
# tr <- list("ACTGAC", "CTGAC")
# calcOverlap(tx2gene = tx2gene, trans_seq_exon_only = tr, test = TRUE)
# start <- proc.time()
# overlaprestemp <- parLapply(clust, genestouse, calcOverlap, tx2gene = tx2gene, trans_seq_exon_only = trans_seq_exon_only)
# overlapres <- t(data.frame(overlaprestemp))
# rownames(overlapres) <- genestouse
# colnames(overlapres) <- "overlap"
# overlaprestime <- proc.time() - start

# save(overlapres, overlaprestime, file = "OverlapRes.RData")










#OLD code below
# trans_seqs2 <- extractTranscriptSeqs(Hsapiens, txdb, use.names=TRUE)
# trans_seqs <- readDNAStringSet(paste0(def_wd, "gencode.v27.transcripts.fa"))
# transcripts <- transcripts(txdb, columns = c("tx_id", "tx_name", "gene_id"))
# tx_names <- transcripts$tx_name
# exons <- exons(txdb, columns = c("exon_id", "exon_name", "tx_id", "tx_name", "gene_id"))
# 
# 
# n <- as.character(lapply(strsplit(names(trans_seqs), split = "\\|"), function(x){x[1]}))
# 
# #List of transcripts
# trans <- DataFrame(transcripts(txdb, columns = c("tx_id", "tx_name", "gene_id")))
# trans2 <- transcripts(txdb, columns = c("tx_id", "tx_name", "gene_id"))
# 
# #Lists of genes and exons
# genes <- genes(txdb)
# 
# 
# overlapPercs <- data.frame(as.numeric(lapply(fullgenenames, calcOverlapPerc, exons = exons)))
# rownames(overlapPercs) <- fullgenenames
# colnames(overlapPercs) <- c("overlapperc")
# 
# #This calculates and returns the percentage of exons within a gene that are present in multiple transcripts
# calcOverlapPerc <- function(x, exons){
#   #x <- "ENSG00000000003.14"
#   #sub2 <- subset(exons, exons$gene_id=="ENSG00000000003.14")
#   sub2 <- subset(exons, exons$gene_id==x)
#   vals <- sub2$tx_id
#   overlapperc <- mean(as.numeric(lapply(vals, function(x) {length(x)>1})))
#   return(overlapperc)
# }
# 
# overlapPercs$gene_id <- rownames(overlapPercs)
# 
# save(overlapPercs, file = "overlapPercs.RData")
# 
# genelengths <- DataFrame(ranges(genes)@width)
# rownames(genelengths) <- genes$gene_id
# colnames(genelengths) <- c("gene_length")
# 
# sub1 <- subset(trans, trans$gene_id=="ENSG00000000003.14")
# 
# 
# exp <- exonicParts(txdb, linked.to.single.gene.only = TRUE)
# sub3 <- subset(exp, exp$gene_id=="ENSG00000000003.14")
# 
# rang <- ranges(sub1)

