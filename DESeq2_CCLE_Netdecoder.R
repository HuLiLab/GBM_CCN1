## Developed by Cristina Correia 2020 ##############
## Updated March 2021 ##############################
## Script for Deseq2 ###############################
####################################################
#=============================================================
#REF. https://bioconductor.riken.jp/packages/3.7/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#countmat
#=============================================================

# setpath ====================================================
setwd("/Users/m001542/Documents/PROJECTs/Li_Lab/DileepDM/GB_Netdecoder/CCLE_DESeq2")
## Results: in CCLE_DESeq2

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

#====================================================================
# Load Packages 
#===================================================================
#rm(list=ls())

BiocManager::install("rlang") # you many not need this!
BiocManager::install(c("gplots","ggplot2", "reshape", "reshape2","ggthemes","ggrepel")) # "grid"
BiocManager::install(c("devtools", "dplyr","CePa")) #,"mygene", "GenomicFeatures"))
BiocManager::install(c("RColorBrewer","ComplexHeatmap", "circlize"))
BiocManager::install("DESeq2")
BiocManager::install("biomaRt")

library(rlang)
library(biomaRt)
library(gplots)
library(ggplot2)
library(reshape)
library(reshape2)
#library(grid)
library(ComplexHeatmap)
library(circlize)
library(ggthemes) # custom plotting themes for ggplot2
library(ggrepel) # prevent overlap of datapoints
library(dplyr)
library(CePa) # read a gct format file
#library(GenomicFeatures)
#library(mygene) # Convert geneIDs with mygene
library(DESeq2)


#=================================================================================
# LOAD data CCLE data for gene or interest CCN1=CYR61
#=================================================================================
# Define cell line labels for DF analysis e.g. CCN1_low and CCN1_high

CCLE.df <- read.csv("./CCLE_mRNACYR61.csv") # 

# Define aggregate statistics to generate labels
sd(CCLE.df$CYR61)
mean(CCLE.df$CYR61)
median(CCLE.df$CYR61)

CCLE.df$Dx <- c("GBM")
CCLE.df$CCN1_Status <- c("High")

# Define cell lines labels
CCLE.df$CCN1_Status <- ifelse(CCLE.df$CYR61 < 7.39, c("Low"), c("High"))

# Determine if the labels classes are balanced
low <- CCLE.df[(CCLE.df$CYR61 > 7.39), ]
length(rownames(low))
write.csv(CCLE.df, "CCLE_RNAseq_GBM_samples.csv") # labels assigned

# Box plot to check CCN1 range
pdf(file="Boxplot_CCN1_Mar2021.pdf", width=4, height=8);
p2 <- ggplot(CCLE.df, aes(x=Dx, y=CYR61, color=c("orange"))) +
  geom_boxplot(fill='white', color="darkblue",outlier.colour="black")+
  geom_jitter(position='identity',size=.5) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3)+
  theme_bw() + xlab("GBM") + ylab("Expression") + theme(legend.position="bottom") +
  theme(axis.text.x = element_text(size=rel(1), hjust=1))
plot(p2)
dev.off()
#=================================================================================


#=============================================================================
# DESeq2 Analysis 
#=============================================================================

# Retrieve CCLE RNAseq count data
eset <- read.gct("/Users/m001542/Desktop/DATA/CCLE_RNAseq_genes_counts_20180929.gct") #use CePa package
head(eset)
eset.df <- as.data.frame(eset) # convert to df
eset.df.lymph <- eset.df[,grepl("CENTRAL_NERVOUS_SYSTEM",colnames(eset.df))] # select only BREAST cell lines,
eset.df.lymph.2 <- eset.df.lymph
colnames(eset.df.lymph.2) <- gsub("(_CENTRAL_NERVOUS_SYSTEM)", "", colnames(eset.df.lymph.2))
write.csv(eset.df.lymph.2, file='CCLE_RNAseq_Counts_GBM_Mar2021.csv')

GBM.df <- eset.df.lymph.2

#GBM.df <- read.csv("./CCLE_RNAseq_Counts_GBM_Mar2021.csv")
#rownames(GBM.df) <- GBM.df$X
#GBM.df <- GBM.df[,-1] # n=65

#write.csv(CCLE.df, "CCLE_RNAseq_GBM_samples.csv") # labels assigned
samples <- read.csv("./Scheme_CCLE_RNAseq_GBM_samples.csv") # n=66
rownames(samples) <- samples$cell #samples$labels
#View(CCLE.df)

a <- intersect(colnames(GBM.df),samples$cell) # common cell lines
b <- setdiff(colnames(GBM.df),samples$cell) # what is different

GBM.df1 <- GBM.df[, intersect(rownames(samples), colnames(GBM.df))]
samples1 <- samples[intersect(colnames(GBM.df1), rownames(samples)),]
GBM.df1 <- GBM.df1[, c(samples1$cell)] # order DF
all(colnames(GBM.df1) == rownames(samples1)) # check that order is correct

#=============================================================================
# DESeq2 design, define default comparisons  
#=============================================================================
coldata <- samples1
coldata <- samples1[,c("condition","group","tx","cell")] # labels
coldata$condition <- factor(coldata$condition) # generate factors for comparisons 
coldata$cell <- factor(coldata$cell) # generate factors for comparisons
coldata$group <- factor(coldata$group)
coldata$tx <- factor(coldata$tx)
#=============================================================================

# Note the expression DF needs to have the order as coldata
GBM.df2 <- GBM.df1[,c(rownames(coldata))] # order the DF as coldata

# Full count matrix # DESeq2 input is a matrix
# DESeq2 input is a matrix of un_normalized counts
GBM.exp <- as.matrix(GBM.df2) #convert to all count matrix for DESeq2
# head(GBM.exp)

# Build a DESeq2 object =====================================================
#colData = coldata.sub, condition
ddsFull <- DESeqDataSetFromMatrix(countData = GBM.exp,
                                colData = coldata,
                               design = ~ tx)
ddsFull # confirm number of samples and names
as.data.frame( colData(ddsFull)) # check colData info
View(counts(ddsFull)) # check matrix
#===============================================================================

# prefilter step
keep <- rowSums(counts(ddsFull)) >= 1 # remove rows with a lot of zeros
ddsFull <- ddsFull[keep,]

set.seed(123)
ddsFull <- estimateSizeFactors(ddsFull) # pre-estimate the size factors so these do no change wth desgn 
ddsFull <- DESeq(ddsFull)
#res <- results(ddsFull) 
#resultsNames(ddsFull)
sf <- sizeFactors(ddsFull)
write.csv(sf, "CCLE_GBM_sizeFactors_Mar_2021.csv")


# adjust to dds as needed )see drop levels)
colData(ddsFull)
View(counts(ddsFull))
dds <- DESeq(ddsFull)
res <- results(ddsFull)
resultsNames(ddsFull)

#=================================================================================================================
# Define comparisons 
#contrast = c('factorName','numeratorLevel','denominatorLevel')
#==================================================================================================================

GBM.res.all <- results(dds, contrast = c("tx", "CCN1_high", "CCN1_low"))#, alpha = 0.05) # alpha for pval.
#TLRB1.res <- results(dds, contrast = c("tx", "T1", "ctr"), alpha = 0.005)
GBM.Sig <- subset(GBM.res.all, padj < 0.05) #pvalue < 0.05)
GBM.Sig <- as.data.frame(GBM.Sig)
write.csv(GBM.res.all, file="./data/GBM_high_lowCCN1_counts.csv")
write.csv(GBM.Sig, file="./data/GBM_high_lowCCN1_p005counts.csv")

# Extract normalized counts ============================================================================================
#Variance stabilizing transformations (VST) 
vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 3)
vsd.df <- assay(vsd) # extract the DESeq2 object with assay

write.csv(vsd.df, file="./data/GBM_high_lowCCN1_normalizedcounts_March2021.csv")
#=================================================================================================================

#=================================================================================================================
## Convert gene IDs - I need to remove the "." in ENSG00000237613.2 - does not allow for ID conversion
#=================================================================================================================
# Remove the "." in ENSG00000237613.2

GBMfinal.df  <- vsd.df
GBMfinal.df  <- read.csv("./data/GBM_high_lowCCN1_normalizedcounts_March2021.csv")
rownames(GBMfinal.df) <- GBMfinal.df$X
GBMfinal.df$ID <- as.character(rownames(GBMfinal.df))

GBMfinal.df$gene2 <- gsub('(.*).\\w+', '\\1', GBMfinal.df$ID) # use gsub to remove all after the "." and keep what's first
# note. note all "." were removed  ENSG00000078808.13  to ENSG00000078808.
#gene_test <- c("ENSG00000078808.")

genes.comma <- as.vector(unlist(lapply(strsplit(as.character(GBMfinal.df$gene2), split="\\."), "[", 1))) 
unique.genes <- make.unique(genes.comma, sep='')
GBMfinal.df$gene3 <- as.vector(unique.genes)
rownames(GBMfinal.df) <- as.vector(GBMfinal.df$gene3) 

# Convert gene IDs with bioMart
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- GBMfinal.df$gene3
GBMfinal.df$id <- NA
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                          "hgnc_symbol", "description"),values=genes,mart= mart)

df.final <- merge(GBMfinal.df,G_list,by.x="gene3",by.y="ensembl_gene_id") # merge dfs
rownames(df.final) <- df.final$hgnc_symbol # rename df.final

#duplicate_genes <- c("GOLGA8M", "ITFG2-AS1", "LINC01238", "RN7SL274P", "RNU6-318P", "TUBB7P")
df.final.2 <- df.final[!duplicated(df.final$hgnc_symbol),] # find and remove duplicate IDs
rownames(df.final.2) <- df.final.2$hgnc_symbol # rename rownames
colnames(df.final.2)
df.final.3 <- df.final.2[, c(3:65)]
write.csv(df.final.3 , file='CCLE_RNAseq_CountsNormalized_GBM_Mar2021.csv')


# ConverT the DF for Significant genes 
GBM.Sig.df <- read.csv("./data/GBM_high_lowCCN1_p005counts.csv")
genesID.df <- df.final.2[, c("gene3","X","ID","gene2","hgnc_symbol", "description")]
write.csv(genesID.df  , file='GBM_GeneConversion_DF.csv')

GBM.Sig.df.2 <-  merge(GBM.Sig.df,genesID.df,by.x="X",by.y="X")
#View(GBM.Sig.df.2)
rownames(GBM.Sig.df.2) <- GBM.Sig.df.2$hgnc_symbol
write.csv(df.final.3 , file='GBM_high_lowCCN1_p005countsGenes.csv')
#==================================================================
