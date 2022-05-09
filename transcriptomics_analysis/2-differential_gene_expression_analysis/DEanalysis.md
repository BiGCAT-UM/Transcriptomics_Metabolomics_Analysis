## Introduction

In this workflow, differential gene expression analysis will be
performed on host-transcriptomics data.

## R environment setup

First, we need to make sure all required packages

``` r
# check if libraries are already installed > otherwise install it
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager",repos = "http://cran.us.r-project.org")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"baySeq" %in% installed.packages()) BiocManager::install("baySeq")
if(!"DESeq2" %in% installed.packages()) BiocManager::install("DESeq2")
if(!"edgeR" %in% installed.packages()) BiocManager::install("edgeR")
if(!"bioDist" %in% installed.packages()) BiocManager::install("bioDist")
if(!"biomaRt" %in% installed.packages()) BiocManager::install("biomaRt")
if(!"dplyr" %in% installed.packages()) BiocManager::install("dplyr")
if(!"magrittr" %in% installed.packages()) BiocManager::install("magrittr")

#load packages
library(rstudioapi)
library(baySeq)
library(DESeq2)
library(edgeR)
library(bioDist)
library(biomaRt)
library(dplyr)
library(ggplot2)
library(limma)
library(gplots)
library(magrittr)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(VennDiagram)

# set your working environment to the location where your current source file is saved into.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#include some functions adapted from ArrayAnalysis.org scripts
source("functions_ArrayAnalysis_v2.R")
WORK.DIR <- getwd()
```

## Data Preparations

The following section will prepar input data to be used in the analysis

``` r
#read  the input data file
htxCount <- read.table(file = 'data/htxCount', sep = '\t', header = TRUE)

#checking which samples have all zero values across all genes
#these sample should be removed otherwise there will be a problem when calculating estimate size factors
idx <- which(colSums(htxCount) == 0)
#CSMDRVXI MSM719ME  are samples who has all zero values so we remove them
htxCount <- htxCount[ , -idx]

#read metadata file sample labels
sampleLabels <- read.table(file = "data/sampleLabels", row.names=1,sep = '\t', stringsAsFactors = TRUE)
#apply same filtering to the metadata remove samples from sample labels metadata
sampleLabels <- sampleLabels[-idx , ]

#add column names to the metadata file
colnames(sampleLabels) <- c( "sampleID", "biopsy_location","disease")
#check whether sample names are in same order
#all(colnames(htxCount) == rownames(sampleLabels))

#select only biopsy_location and disease columns
sampleLabels<-sampleLabels[, c(2,3)]
sampleLabels$disease <- relevel(sampleLabels$disease,ref="nonIBD")
#add an experimental group variable to sampleLabels
sampleLabels$group <- as.factor(paste(sampleLabels$disease,sampleLabels$biopsy_location,sep="_"))
```

## Filtering Steps

We will apply some filtering process to filter out genes in the input
data

``` r
#remove genes which has all zero values for all samples then start DE analysis
nonzero <- rowSums(htxCount) > 0
htxCount %<>% .[nonzero,]

#############################CPM FILTERING#############################################
#aveLogCPM function compute average log2 counts-per-million for each row of counts.
#the below function is similar to log2(rowMeans(cpm(y, ...)))
mean_log_cpm = aveLogCPM(htxCount)

# We plot the distribution of average log2 CPM values to verify that our chosen presence threshold is appropriate. The distribution is expected to be bimodal, with a low-abundance peak representing non-expressed genes and a high-abundance peak representing expressed genes. The chosen threshold should separate the two peaks of the bimodal distribution. 
filter_threshold <- -1# we can try different threshold values
#jpeg(file="avgLogCpmDist.jpeg")#if you want to save the histogram uncomment the following command  
ggplot() + aes(x=mean_log_cpm) +
    geom_histogram(binwidth=0.2) +
    geom_vline(xintercept=filter_threshold) +
    ggtitle("Histogram of mean expression values")
```

![](DEanalysis_files/figure-markdown_github/filtering-1.png)

``` r
#dev.off()#to save the plot to the file
#Having chosen our threshold, lets pick the subset of genes whose average expression passes that threshold.
keep_genes <- mean_log_cpm >= filter_threshold 
htxCount <- htxCount[keep_genes,]
#dim(htxCount)#to check dimension of the data
###############################################################################
```

## Differential Gene Expression Analysis

In the following section differential gene expression analysis will be
performed using DESeq2 package

``` r
#we wil firstly create a DESeqDataSet object
#(non-intercept) statistical model based on the disease and biopsy_location, group column represent both of them 
dds <- DESeqDataSetFromMatrix(countData = htxCount, colData=sampleLabels, design= ~0 + group)

#estimate the size factors
#To perform the median of ratios method of normalization, DESeq2 has a single estimateSizeFactors() function that will generate size factors for us.
dds <- estimateSizeFactors(dds)
#normalise the data (here for Quality Control(QC) plotting)
#QC plotting is optional
norm <- counts(dds,normalize=TRUE)
#create a 2logged data for original object (here for QC plotting)
datlog <- log(htxCount+1,2)
#create a 2logged norm object (here for QC plotting)
normlog <- log(norm+1,2)
#for QC remove genes that have not been measured in any sample in the experiment
datlogQC <- datlog[rowSums(datlog)!=0,]
normlogQC <- normlog[rowSums(normlog)!=0,]
#create QC plots for raw data, colored by different variables
factors <- c("disease","biopsy_location","group")
if(!dir.exists("QCraw")) dir.create("QCraw")
setwd(paste(WORK.DIR,"QCraw",sep="/"))
png("sizefactors.png")
plot(sizeFactors(dds),type='h',lwd=5,ylim=c(0,max(sizeFactors(dds))),col="darkblue")
dev.off()
```

    ## png 
    ##   2

``` r
createQCPlots(datlogQC, factors, Table=sampleLabels, normMeth="", postfix="")
setwd("..")
#create QC plots for normalized data coloured by different variables
if(!dir.exists("QCnorm")) dir.create("QCnorm")
setwd(paste(WORK.DIR,"QCnorm",sep="/"))
createQCPlots(normlogQC, factors, Table=sampleLabels, normMeth="DESeq", postfix="")
setwd("..")
#sample MSM719M9 is an outlier remove it from dataset
#sample HSM5FZAZ is an outlier remove it from dataset
htxCount <- htxCount[,-match(c("MSM719M9","HSM5FZAZ"),colnames(htxCount))]
sampleLabels <- sampleLabels[-match(c("MSM719M9","HSM5FZAZ"),rownames(sampleLabels)),]
#doublecheck whether the order of the samples in samplekey and dat still match
#sum(rownames(sampleLabels) == colnames(htxCount))==dim(sampleLabels)[1]

###REDO ALL STEPS AS ABOVE FOR THE NEW DATA SET
dds <- DESeqDataSetFromMatrix(countData = htxCount, colData=DataFrame(sampleLabels), design= ~0 + group)
dds <- estimateSizeFactors(dds)
norm <- counts(dds,normalize=TRUE)
datlog <- log(htxCount+1,2)
normlog <- log(norm+1,2)
datlogQC <- datlog[rowSums(datlog)!=0,]
normlogQC <- normlog[rowSums(normlog)!=0,]
#create QC plots for raw data, coloured by different variables
if(!dir.exists("QCraw2")) dir.create("QCraw2")
setwd(paste(WORK.DIR,"QCraw2",sep="/"))
png("sizefactors2.png")
plot(sizeFactors(dds),type='h',lwd=5,ylim=c(0,max(sizeFactors(dds))),col="darkblue")
dev.off()
```

    ## png 
    ##   2

``` r
createQCPlots(datlogQC, factors, Table=sampleLabels, normMeth="", postfix="")
setwd("..")
#create QC plots for normalized data coloured by different variables
if(!dir.exists("QCnorm2")) dir.create("QCnorm2")
setwd(paste(WORK.DIR,"QCnorm2",sep="/"))
createQCPlots(normlogQC, factors, Table=sampleLabels, normMeth="DESeq", postfix="")
setwd("..")

################################################################################
#######################statistical modelling####################################
################################################################################
#set directory for stat output
if(!dir.exists("statsmodel")) dir.create("statsmodel")
setwd(paste(WORK.DIR,"statsmodel",sep="/"))
#run differential analysis
dds <- DESeq(dds)
#results(dds) #runs the default DESeq2 comparison    (last factor level versus first factor level)
cont.matrix <- makeContrasts(
  #CD disease on ileum and rectum 
  CD_Ileum_vs_nonIBD_Ileum   = groupCD_Ileum - groupnonIBD_Ileum,
  CD_Rectum_vs_nonIBD_Rectum = groupCD_Rectum - groupnonIBD_Rectum,
  #UC disease on ileum and rectum
  UC_Ileum_vs_nonIBD_Ileum   = groupUC_Ileum - groupnonIBD_Ileum,
  UC_Rectum_vs_nonIBD_Rectum = groupUC_Rectum - groupnonIBD_Rectum,
  #UC and CD disease comparison
  UC_Ileum_vs_CD_Illeum     = groupUC_Ileum - groupCD_Ileum,
  UC_Rectum_vs_CD_Rectum    = groupUC_Rectum - groupCD_Rectum,
  #biopsy location comparisons
  CD_Rectum_vs_CD_Ileum     = (groupCD_Rectum - groupnonIBD_Rectum) - (groupCD_Ileum - groupnonIBD_Ileum),
  UC_Rectum_vs_UC_Ileum     = (groupUC_Rectum - groupnonIBD_Rectum) - (groupUC_Ileum - groupnonIBD_Ileum),
  levels = resultsNames(dds)
)
#extract resulting contrasts based on the model, and save those in a table; also save some graphical representations
#the function results() is called from within the saveStatOutputDESeq2 function to compute the contrasts
files <- saveStatOutputDESeq2(cont.matrix,dds,postfix="",annotation=NULL)
```

    ## --[[ Saving table for coefficient    CD_Ileum_vs_nonIBD_Ileum     ]]--
    ## ----[[  0  probes with NA estimates removed ]]

    ## --[[ Saving table for coefficient    CD_Rectum_vs_nonIBD_Rectum   ]]--
    ## ----[[  0  probes with NA estimates removed ]]

    ## --[[ Saving table for coefficient    UC_Ileum_vs_nonIBD_Ileum     ]]--
    ## ----[[  0  probes with NA estimates removed ]]

    ## --[[ Saving table for coefficient    UC_Rectum_vs_nonIBD_Rectum   ]]--
    ## ----[[  0  probes with NA estimates removed ]]

    ## --[[ Saving table for coefficient    UC_Ileum_vs_CD_Illeum    ]]--
    ## ----[[  0  probes with NA estimates removed ]]

    ## --[[ Saving table for coefficient    UC_Rectum_vs_CD_Rectum   ]]--
    ## ----[[  0  probes with NA estimates removed ]]

    ## --[[ Saving table for coefficient    CD_Rectum_vs_CD_Ileum    ]]--
    ## ----[[  0  probes with NA estimates removed ]]

    ## --[[ Saving table for coefficient    UC_Rectum_vs_UC_Ileum    ]]--
    ## ----[[  0  probes with NA estimates removed ]]

``` r
#create summary table of the contrast results
#up and down for p-val and adj p-val are decided by for example log2fc>=0.58 and  log2fc<0.58
createPvalTab(files,postfix="",namePVal="pvalue",nameAdjPVal="padj",nameFC="FoldChange",nameLogFC="log2FoldChange",html=TRUE)
```

    ## --[[ Saving Summary_tables.tab ]]--
    ## --[[ Saving Summary_tables.html ]]--

``` r
setwd("..")
```

## Enhanced Volcano Plots Visualization

Differential gene expression analysis results will be visualized by
volcano plots

``` r
#go to the result folder named statsmodel
setwd(paste(WORK.DIR,"statsmodel",sep="/"))
readFilePath <- paste(WORK.DIR,"statsmodel",sep="/")
# create an empty list
plot_list = list()
for(i in 1:length(files))
{
  #read file
  splitted <- strsplit(files[i], split = "_")[[1]]
  title <- splitted [2:6]
  title <- paste (splitted[2],splitted[3],splitted[4],splitted[5],splitted[6],sep=" ") 
  title <- (strsplit(title, split = "\\.")[[1]])[1]
  tab <- read.delim(paste(readFilePath, files[i],sep="/"),header=TRUE,as.is=TRUE)  
  p<-EnhancedVolcano(tab , lab = tab$X, labSize = 3, title = title, 
                     x = 'log2FoldChange', y = 'pvalue', pCutoff = 0.05, FCcutoff = 0.58)
  plot_list[[i]] = p
}
#path of the outout folder
outFolder <- paste(WORK.DIR,"volcano_plots",sep="/")
if(!dir.exists(outFolder)) dir.create(outFolder)

for(i in 1:length(files)) {
    file_name = paste(outFolder,"/",files[i],".png",sep="")
    png(file_name)
    print(plot_list[[i]])
    dev.off()
}  
```
