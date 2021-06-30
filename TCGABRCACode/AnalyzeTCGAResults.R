#Analyze TCGA BRCA Results
library(readr)
library(compositions)
library(tikzDevice)
library(xtable)

fr <- function(x, n = getOption("digits"), forLatex = TRUE, NLS = FALSE) {
  y <- rep(0, length(x))
  for (i in 1:length(x)) {
    if (x[i] >=0 & x[i]< 10^-(n) & !is.na(x[i])){
      
      if(NLS==TRUE){
        y[i] <- paste("0", ".", paste(rep("0", n), collapse=""), sep = "")
      }else if(forLatex==TRUE){
        y[i] <- paste("$<$", ".", paste(rep("0", n-1), collapse=""), "1", sep = "")
      }else{
        y[i] <- paste("<", ".", paste(rep("0", n-1), collapse=""), "1", sep = "")
      }
      
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

"Interestingly, GSEA of differentially expressed isoforms demonstrated enrichment
for genes involved in apoptosis, and thus differential splicing may contribute more broadly
than CASP2 to differences in regulation of apoptosis between breast cancer subtypes"

gtf <- rtracklayer::import("/Users/Scott/Documents/Dissertation Data/GEUV1Data/gencode.v27.annotation.gtf.gz")

#Load v32 annotation for comparison, byt v27 was used for the analysis so don't use this one
#gtf_v32 <- rtracklayer::import("/Users/Scott/Documents/Dissertation Data/SingleCellProject/gencode_v32_files/gencode.v32.annotation.gtf.gz")
#gtf_v32_CASP2_AllTrans <- subset(gtf, gtf$gene_name=="CASP2"& gtf$type=="transcript")
#gtf2 <- gtf[,c("gene_id", "gene_name")]

gtf2 <- data.frame(gtf$gene_id, gtf$gene_name, stringsAsFactors = F)
colnames(gtf2) <- c("gene_id", "gene_name")

gtf3 <- unique(gtf2)

full_trans <- subset(gtf, gtf$type=="transcript")

gtf_transT <- data.frame(gtf$gene_id, gtf$gene_name, gtf$transcript_id, gtf$transcript_name, stringsAsFactors = F)
colnames(gtf_transT) <- c("gene_id", "gene_name", "transcript_id", "transcript_name")

gtf_transT2 <- unique(gtf_transT)
load("/Users/Scott/Documents/Dissertation Data/TCGABRCAAnalysis/abDatasetsNoOtherGroupsFiltered.RData")
load("/Users/Scott/Documents/Dissertation Data/TCGABRCAAnalysis/abGeneFiltered.RData")

load("/Users/Scott/Documents/Dissertation Data/TCGABRCAAnalysis/key.RData")

#Genes that make up the pam50 subtype that determines brest cancer subtypes
pam50 = c(
  "ACTR3B",
  "ANLN",
  "BAG1",
  "BCL2",
  "BIRC5",
  "BLVRA",
  "CCNB1",
  "CCNE1",
  "CDC20",
  "CDC6",
  "CDH3",
  "CENPF",
  "CEP55",
  "CXXC5",
  "EGFR",
  "ERBB2",
  "ESR1",
  "EXO1",
  "FGFR4",
  "FOXA1",
  "FOXC1",
  "GPR160",
  "GRB7",
  "KIF2C",
  "KRT14",
  "KRT17",
  "KRT5",
  "MAPT",
  "MDM2",
  "MELK",
  "MIA",
  "MKI67",
  "MLPH",
  "MMP11",
  "MYBL2",
  "MYC",
  "NAT1",
  "NDC80",
  "NUF2",
  "ORC6L",
  "PGR",
  "PHGDH",
  "PTTG1",
  "RRM2",
  "SFRP1",
  "SLC39A6",
  "TMEM45B",
  "TYMS",
  "UBE2C",
  "UBE2T"
)



BasalVsNonBasal <- FALSE

if(BasalVsNonBasal==TRUE){
  load("/Users/Scott/Documents/Dissertation/res/TCGABRCAAnalysis/CompDTURes/CompDTUObsResBasalVsNonBasal.RData")
  load("/Users/Scott/Documents/Dissertation/res/TCGABRCAAnalysis/DRIMSeqResultsBasalVsNonBasal/DRIMSeqResDRIMSeqFiltering.RData")
  load("/Users/Scott/Documents/Dissertation/res/TCGABRCAAnalysis/RATsResultsBasalVsNonBasal/RATsObsRes.RData")
  load("/Users/Scott/Documents/Dissertation/res/TCGABRCAAnalysis/RATsResultsBasalVsNonBasal/RATsObsResInfReps.RData")
  key$group <- relevel(as.factor(ifelse(key$Condition=="Basal", "A", "B")), ref = "A")
}else{
  load("/Users/Scott/Documents/Dissertation/res/TCGABRCAAnalysis/CompDTURes/CompDTUObsRes.RData")
  load("/Users/Scott/Documents/Dissertation/res/TCGABRCAAnalysis/DRIMSeqResults/DRIMSeqResDRIMSeqFiltering.RData")
  key$group <- key$Condition
}

load("/Users/Scott/Documents/Dissertation Data/TCGABRCAAnalysis/InfRVResAllForBootSamps.RData")
InfRVToUse <- InfRVResAllBootSamps[,c("gene_id", "MeanGeneMaxInfRVCnt")]
colnames(InfRVToUse) <- c("gene_id", "InfRV")

CompDTU <- CompDTUResults[,c("gene_id", "pval_CompDTU")]
CompDTU$adj_pval_CompDTU <- p.adjust(CompDTU$pval_CompDTU, method = "fdr")

CompDTUme <- CompDTUmeResults[,c("gene_id", "pval_CompDTUme")]
CompDTUme$adj_pval_CompDTUme <- p.adjust(CompDTUme$pval_CompDTUme, method = "fdr")

DRIMPvals$gene_id <- rownames(DRIMPvals)
DRIMSeq <- DRIMPvals[,c("gene_id", "pvalue")]
colnames(DRIMSeq) <- c("gene_id", "pval_DRIMSeq")
DRIMSeq$adj_pval_DRIMSeq <- p.adjust(DRIMSeq$pval_DRIMSeq, method = "fdr")
if(BasalVsNonBasal==TRUE){
  RATs <- res_reg
  #Drop these adj pvalues and recalculate them yourself just to confirm how they are calculated
  RATs$adj_pvalue <- NULL
  RATs$adj_pvalue_RATs <- p.adjust(RATs$pval, method = "fdr")
  RATs$pval <- NULL
  
  RATsInfRep <- res_infRep
  RATsInfRep$adj_pvalue_RATsInfRep <- p.adjust(RATsInfRep$pval_infReps, method = "fdr")
  RATsInfRep$pval_infReps <- NULL
  
  final_resT <- Reduce(function(x, y) merge(x, y, all=TRUE), list(CompDTU, CompDTUme, DRIMSeq, RATs, RATsInfRep))
  all_methods_resT <- merge(final_resT, gtf3, by = "gene_id")
  all_methods_res <- merge(all_methods_resT, InfRVToUse, by = "gene_id")
  write.csv(all_methods_res, file = "all_methods_res_basal_non_basal.csv")
}else{
  final_resT <- Reduce(function(x, y) merge(x, y, all=TRUE), list(CompDTU, CompDTUme, DRIMSeq))
  all_methods_resT <- merge(final_resT, gtf3, by = "gene_id")
  all_methods_res <- merge(all_methods_resT, InfRVToUse, by = "gene_id")
  write.csv(all_methods_res, file = "all_methods_res.csv")
}

#Computation Time Table Values
sum(CompDTUResults$ComputationTime)
length(!is.na(CompDTUResults$pval_CompDTU))
sum(CompDTUResults$ComputationTime)/length(!is.na(CompDTUResults$pval_CompDTU))
sum(CompDTUmeResults$ComputationTime)
length(!is.na(CompDTUmeResults$pval_CompDTUme))
sum(CompDTUmeResults$ComputationTime)/length(!is.na(CompDTUmeResults$pval_CompDTU))
DRIMCompTime[3]
length(!is.na(DRIMSeq$pval_DRIMSeq))
DRIMCompTime[3]/length(!is.na(DRIMSeq$pval_DRIMSeq))


#Lists of genes of interest - a list of genes associated with BRCA from genecards.org
GeneCardsBRCA <- data.frame(read_csv("/Users/Scott/Documents/Dissertation/res/TCGABRCAAnalysis/GeneCardsBRCA.csv"))
#genes_of_interest <- c("BRCA1", "BRCA2", "BRCA3", "BRCA1P1", "RAD51", "BRIP1", "BARD1", "BAP1")

#These genes come from "Robust stratification of breast cancer subtypes using differential patterns of transcript isoform expression"
genes_of_interest <- c("TP53", "BRCA1", "PTEN", "CD44", "ECHDC1", "PPP2R5D", "CASP2", "YBX1", "YBX2", "MAGOH", "MAGOHB", "PCBP2")

#This list is from "Widespread alternative exon usage in clinically distinct subtypes of Invasive Ductal Carcinoma"
  #https://www.nature.com/articles/s41598-017-05537-0
genes_of_interest2 <- c("EPB41L1", "TPD52", "ACOX2", "MYO6", "IQCG")

fun1 <- function(gene, returnAvgProps = TRUE){
  library(compositions)
  gene_idd <- gtf3$gene_id[gtf3$gene_name==gene]
  lol <- abDatasetsFiltered[[gene_idd]]
  n_trans <- ncol(lol)
  lol2 <- data.frame(as.matrix(unclass(ccomp(lol, total = 1))))
  lol2$group <- key$group
  
  lol$group <- key$group
  
  if(returnAvgProps==FALSE){
    t1 <-aggregate(lol[,1:n_trans], by = list(group = lol$group), FUN = mean)
    t1$gene_id <- gene_idd
    return(t1)
  }else{
    t1 <- aggregate(lol2[,1:n_trans], by = list(group = lol2$group), FUN = mean)
    t1$gene_id <- gene_idd
    return(t1)
  }

}

genes_of_interest3 <- c("PIK3CA", "ERBB2", "CDKN2A")

genes_of_interest4 <- c("CTNND1", "PRICKLE1")

#These genes come from "Alternative Splicing in Breast Cancer and the Potential Development of Therapeutic Tools"
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5664086/
genes_of_interest6 <- c("BRCA1", "DMTF1", "HER2", "FGFR", "KLF6", "BIRC5", "TP53")

#These genes were identified to have major isoform switching events in breast cancer
  #https://link.springer.com/article/10.1186/s12864-016-2521-9
genes_of_interest7 <- c("CTNND1","PRICKLE1")

#These genes came from Naim (with 3 out of the 4 contained in the PAM50 cluster)
genes_of_interest8 <- c("ERBB2", "ESR1", "AURKA","PGR")

#These genes are discussed specifically in the paper
genes_of_interest9 <- c("MYO6", "EPB41L1", "TPD52","IQCG", "ACOX2", "PRICKLE1", "CTNND1", "CASP2", "KIF2C")

sub_res_GeneCards <- subset(all_methods_res, all_methods_res$gene_name %in% GeneCardsBRCA$Gene.Symbol)

sub_res <- subset(all_methods_res, all_methods_res$gene_name %in% genes_of_interest)
sub_res2 <- subset(all_methods_res, all_methods_res$gene_name %in% genes_of_interest2)
sub_res3 <- subset(all_methods_res, all_methods_res$gene_name %in% genes_of_interest3)
sub_res4 <- subset(all_methods_res, all_methods_res$gene_name %in% genes_of_interest4)

sub_res5 <- subset(all_methods_res, all_methods_res$gene_name %in% pam50)
sum(sub_res5$adj_pval_CompDTU < 0.05, na.rm=T)
sum(sub_res5$adj_pval_CompDTUme < 0.05, na.rm=T)
sum(sub_res5$adj_pval_DRIMSeq < 0.05, na.rm=T)
sub_res6 <- subset(all_methods_res, all_methods_res$gene_name %in% genes_of_interest6)
sub_res7 <- subset(all_methods_res, all_methods_res$gene_name %in% genes_of_interest7)
sub_res8 <- subset(all_methods_res, all_methods_res$gene_name %in% genes_of_interest8)

sub_res9 <- subset(all_methods_res, all_methods_res$gene_name %in% genes_of_interest9)

fun1("PTTG1")
fun1("PTTG1", returnAvgProps = F)

CASP2_props <- fun1("CASP2")
fun1("CASP2", returnAvgProps = F)

#Now, get 95% CI Values for each proportion
Prop_CI_fn <- function(prop){
  z_val <- qnorm(.975)
  p1 <- prop * (1 - prop)
  lower <- prop - (z_val * sqrt(p1/740))
  upper <- prop + (z_val * sqrt(p1/740))
  return(paste0("(", fr(lower, 3), ", ", fr(upper, 3), ")"))
}
Prop_CI_fn(CASP2_props[,2])
Prop_CI_fn(CASP2_props[,3])

fun1("KIF2C")
fun1("KIF2C", returnAvgProps = F)

#These transcripts correspond to CASP2, the one ending on 447.9 (name CASP2-201) is the longer one and the one ending in 992.4 is the shorter one (name CASP2-207)
#See here for more information http://apr2018.archive.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000106144;r=7:143288215-143307696;t=ENST00000310447
gtf_transT3 <- subset(gtf_transT2, gtf_transT2$transcript_id%in%c("ENST00000310447.9", "ENST00000619992.4"))
gtf_transT4 <- subset(gtf, gtf$transcript_name%in%gtf_transT3$transcript_name & gtf$type=="transcript")

trans1 <- subset(gtf_transT2, gtf_transT2$gene_name=="CASP2")[,"transcript_id"]
gtf_trans_CASP2_AllTrans <- subset(gtf, gtf$transcript_id%in%trans1& gtf$type=="transcript")




#Now, generate plot that compares CompDTU vs CompDTUme results for the genes from the PAM50 cluster that show improvement from uncertainty
#Use these lines to split the filtered genes into parts to figure out which parts the example genes below come from and thus which data parts need to be loaded
nparts <- 100
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

filteredgenenames_split <- chunk2(sort(names(abDatasetsFiltered)), nparts)

#These genes correspond to the KIF2C and PTTG1 genes respectively
plot_genes_to_use <- c("ENSG00000142945.12", "ENSG00000164611.12")



for(zz in 1:length(filteredgenenames_split)){
  
  if(plot_genes_to_use[1] %in% filteredgenenames_split[[zz]]){
    KIF2CPartNum <- zz
  }
  
  if(plot_genes_to_use[2]%in% filteredgenenames_split[[zz]]){
    PTTG1PartNum <- zz
  }
}

#After looking, KIF2C is part 50 and PTTG1 is part 66
load("/Users/Scott/Documents/Dissertation Data/TCGABRCAAnalysis/abDatasetsBootPart50.RData")
load("/Users/Scott/Documents/Dissertation Data/TCGABRCAAnalysis/abDatasetsBootPart66.RData")

KIF2CMajTrans <- attr(abDatasetsFiltered[["ENSG00000142945.12"]], "MajorTrans")
PTTG1MajTrans <- attr(abDatasetsFiltered[["ENSG00000164611.12"]], "MajorTrans")
t1 <- subset(gtf, gtf$transcript_id==KIF2CMajTrans)
KIF2CMajTransName <- unique(t1$transcript_name)
t2 <- subset(gtf, gtf$transcript_id==PTTG1MajTrans)
PTTG1MajTransName <- unique(t2$transcript_name)

KIF2CDat <- as.data.frame(unclass(ccomp(abDatasetsFiltered[["ENSG00000142945.12"]], total = 1)))
KIF2CDat <- KIF2CDat[,KIF2CMajTrans, drop = FALSE]
colnames(KIF2CDat) <- "Prop"

PTTG1Dat <- as.data.frame(unclass(ccomp(abDatasetsFiltered[["ENSG00000164611.12"]], total = 1)))
PTTG1Dat <- PTTG1Dat[,PTTG1MajTrans, drop = FALSE]
colnames(PTTG1Dat) <- "Prop"
  

KIF2CFullBoot <- as.data.frame(unclass(ccomp(abDatasetsBootPart50[["ENSG00000142945.12"]], total = 1)))
KIF2CFullBoot <- KIF2CFullBoot[,KIF2CMajTrans, drop = FALSE]
colnames(KIF2CFullBoot) <- "Prop"


PTTG1FullBoot <- as.data.frame(unclass(ccomp(abDatasetsBootPart66[["ENSG00000164611.12"]], total = 1)))
PTTG1FullBoot <- PTTG1FullBoot[,PTTG1MajTrans, drop = FALSE]
colnames(PTTG1FullBoot) <- "Prop"

KIF2CDat$group <- key$group
PTTG1Dat$group <- key$group
KIF2CDat$BootSamps <- 0
PTTG1Dat$BootSamps <- 0

KIF2CFullBoot$group <- rep(key$group, 100)
PTTG1FullBoot$group <- rep(key$group, 100)

KIF2CFullBoot$BootSamps <- 1
PTTG1FullBoot$BootSamps <- 1

KIF2CAllDat <- rbind(KIF2CDat, KIF2CFullBoot)
PTTG1AllDat <- rbind(PTTG1Dat, PTTG1FullBoot)
library(ggplot2)
my_xlab <- c("Regular Abundance Estimates", "All Bootstrap Replicates")
curr_save_dir <- "/Users/Scott/Documents/Dissertation/res/TCGABRCAAnalysis/Plots/"

if(!dir.exists(curr_save_dir)){dir.create(curr_save_dir)}

options("tikzDocumentDeclaration" = "\\documentclass[12pt]{article}\n")
options("tikzMetricsDictionary" = paste(curr_save_dir, "/TikzTemp", sep = ""))
stdal <- T

fil_piece <- "KIF2CPlot"
fil_name <- paste0(curr_save_dir, fil_piece,".tex")
#width=6.5
tikz(file = fil_name, height=6.5, width=13, standAlone = stdal, sanitize = F)
pl <- ggplot(data=KIF2CAllDat,
       aes(x=as.factor(BootSamps),y=Prop,fill=as.factor(group)))+geom_boxplot()+
  ggtitle(paste0("Boxplots of the RTA for ", KIF2CMajTransName, " (Major Transcript of ", "KIF2C", ")"))+
  xlab("")+ylab(paste0("RTA of ", KIF2CMajTransName))+
  scale_x_discrete(labels=my_xlab)+
  scale_fill_discrete(name="Subtype")+
  theme(legend.text = element_text(size=rel(2), face = "bold"),legend.title = element_text(size=rel(2), face = "bold"),axis.text.y=element_text(size=rel(2), face = "bold"),
        axis.text.x=element_text(size=rel(2.5), face = "bold"), plot.title=element_text(size=rel(2), hjust=.5,face = "bold"), axis.title=element_text(size=rel(2), face = "bold"),
        legend.key.size = unit(3, "line"))

print(pl)



dev.off()
system(paste0("pdflatex ","-interaction=nonstopmode ", "-output-directory=", curr_save_dir, " ", fil_name), ignore.stdout = TRUE)
system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".aux"))
system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".log"))
system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".tex"))
rm(fil_name)
rm(fil_piece)


fil_piece <- "PTTG1Plot"
fil_name <- paste0(curr_save_dir, fil_piece,".tex")
tikz(file = fil_name, height=6.5, width=13, standAlone = stdal, sanitize = F)

pl2 <- ggplot(data=PTTG1AllDat,
       aes(x=as.factor(BootSamps),y=Prop,fill=as.factor(group)))+geom_boxplot()+
  ggtitle(paste0("Boxplots of the RTA for ", PTTG1MajTransName, " (Major Transcript of ", "PTTG1", ")"))+
  xlab("")+ylab(paste0("RTA of ", PTTG1MajTransName))+
  scale_x_discrete(labels=my_xlab)+
  scale_fill_discrete(name="Subtype")+
  theme(legend.text = element_text(size=rel(2), face = "bold"),legend.title = element_text(size=rel(2), face = "bold"),axis.text.y=element_text(size=rel(2), face = "bold"),
        axis.text.x=element_text(size=rel(2.5), face = "bold"), plot.title=element_text(size=rel(2), hjust=.5,face = "bold"), axis.title=element_text(size=rel(2), face = "bold"),
        legend.key.size = unit(3, "line"))

print(pl2)

dev.off()
system(paste0("pdflatex ","-interaction=nonstopmode ", "-output-directory=", curr_save_dir, " ", fil_name), ignore.stdout = TRUE)
system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".aux"))
system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".log"))
system(paste0("rm ", paste0(curr_save_dir, fil_piece), ".tex"))
rm(fil_name)
rm(fil_piece)


Pam50DTUSignT1 <- sub_res5
Pam50DTUSignT2 <- Pam50DTUSignT1[order(Pam50DTUSignT1$gene_name),]
Pam50DTUSignT3 <- Pam50DTUSignT2[,c("gene_name", "adj_pval_CompDTU", "adj_pval_CompDTUme", "adj_pval_DRIMSeq")]
colnames(Pam50DTUSignT3) <- c("Gene Name", "CompDTU", "CompDTUme", "DRIMSeq")

Pam50DTUSignT4 <- Pam50DTUSignT3
Pam50DTUSignT4[,2] <- fr(Pam50DTUSignT3[,2], n=6)
Pam50DTUSignT4[,3] <- fr(Pam50DTUSignT3[,3], n=6)
Pam50DTUSignT4[,4] <- fr(Pam50DTUSignT3[,4], n=6)

curr_table_dir <- "/Users/Scott/Documents/Dissertation/res/TCGABRCAAnalysis/Tables/"

if(!dir.exists(curr_table_dir)){dir.create(curr_table_dir)}
TCGABRCAPAM50Table <- xtable(Pam50DTUSignT4, align=rep("c", ncol(Pam50DTUSignT4)+1),caption = c("FDR adjusted $p$-values for the genes from the PAM50 classifier that pass filtering for the DTU analysis."), label = "table:TCGABRCAPAM50")

print(TCGABRCAPAM50Table, caption.placement="top", include.rownames = FALSE, sanitize.text.function=function(x){x}, scalebox = 1, table.placement = "!t", hline.after=c(-1,-1,0,nrow(TCGABRCAPAM50Table),nrow(TCGABRCAPAM50Table)), file = paste0(curr_table_dir, "TCGABRCAPAM50Table.tex"))


