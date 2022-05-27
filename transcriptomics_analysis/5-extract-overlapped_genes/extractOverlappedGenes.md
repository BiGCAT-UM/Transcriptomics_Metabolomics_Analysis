## Introduction

In this workflow, extracting overlapped differential expressed genes
between Chron\`s disease and ulcerative colitis will be performed

## R environment setup

``` r
# check if libraries are already installed > otherwise install it
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"org.Hs.eg.db" %in% installed.packages()) BiocManager::install("org.Hs.eg.db")  
if(!"dplyr" %in% installed.packages()) BiocManager::install("dplyr")

#loading installed libraries
library(rstudioapi) 
library(org.Hs.eg.db) 
library (dplyr)

# set your working environment to the location where your current source file is saved into.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
```

## Importing dataset and filtering

``` r
#Obtain data from step 2:
setwd('..')
## Select a disorder to analyse (options; CD or UC)
location <- "rectum"

#we will first perform process for ileum and then rectum 
#read dataset to be processed for ileum biopsy location
if (location == "ileum") {
dataset.CD <- read.delim("2-differential_gene_expression_analysis/statsmodel/table_CD_Ileum_vs_nonIBD_Ileum.tab")
dataset.UC <- read.delim("2-differential_gene_expression_analysis/statsmodel/table_UC_Ileum_vs_nonIBD_Ileum.tab")
 print("Selected location is ileum")
}else if(location == "rectum"){ 
#read dataset for for rectum biopsy location
dataset.CD <- read.delim("2-differential_gene_expression_analysis/statsmodel/table_CD_Rectum_vs_nonIBD_Rectum.tab")
dataset.UC <- read.delim("2-differential_gene_expression_analysis/statsmodel/table_UC_Rectum_vs_nonIBD_Rectum.tab")
 print("Selected location is rectum")
}
```

    ## [1] "Selected location is rectum"

``` r
#filter out unused columns
dataset.CD <- subset( dataset.CD, select = c(1,3,7) )
dataset.UC <- subset( dataset.UC, select = c(1,3,7) )
#change column name
colnames(dataset.CD)[1] = "gene_symbol"
colnames(dataset.UC)[1] = "gene_symbol"

#return back to root folder
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
```

## Getting ENTREZ IDs gor gene symbols

``` r
hs <- org.Hs.eg.db
entrezID <- AnnotationDbi::select(hs, 
                                  keys = dataset.CD$gene_symbol,
                                  columns = c("ENTREZID", "SYMBOL"),
                                  keytype = "SYMBOL")
#filter out double gene symbols
entrezID<- entrezID %>% distinct(entrezID$SYMBOL, .keep_all = TRUE)

# add entrezIDs to each dataset
dataset.CD <- cbind(entrezID$ENTREZID,dataset.CD)
dataset.UC <- cbind(entrezID$ENTREZID,dataset.UC)
#change column names 
colnames(dataset.CD)[1] = "ENTREZ.ID"
colnames(dataset.UC)[1] = "ENTREZ.ID"
```

## Getting differential expressed genes for each disease

``` r
#list of all deg from CD
deg.CD  <-unique(dataset.CD[!is.na(dataset.CD$pvalue) & dataset.CD$pvalue < 0.05 & abs(dataset.CD$log2FoldChange) > 0.58,])
CD.up   <-unique(dataset.CD[!is.na(dataset.CD$pvalue) & dataset.CD$pvalue < 0.05 & dataset.CD$log2FoldChange > 0.58,])
CD.down <-unique(dataset.CD[!is.na(dataset.CD$pvalue) & dataset.CD$pvalue < 0.05 & dataset.CD$log2FoldChange < -0.58,])
#list of all deg from UC 
deg.UC  <-unique(dataset.UC[!is.na(dataset.UC$pvalue) & dataset.UC$pvalue < 0.05 & abs(dataset.UC$log2FoldChange) > 0.58,])
UC.up   <-unique(dataset.UC[!is.na(dataset.UC$pvalue) & dataset.UC$pvalue < 0.05 & dataset.UC$log2FoldChange > 0.58,])
UC.down <-unique(dataset.UC[!is.na(dataset.UC$pvalue) & dataset.UC$pvalue < 0.05 & dataset.UC$log2FoldChange < -0.58,])
```

## Finding overlapped genes between diseases on ileum biopsy location

``` r
######################################FOR ILEUM biopsy location#######################################
# overlap genes between CD down and UC down
overlap.genes1   <- CD.down [CD.down$gene_symbol %in% intersect(CD.down$gene_symbol,UC.down$gene_symbol), c(1:2)]
CD.down.filtered <- CD.down [CD.down$gene_symbol %in% overlap.genes1$gene_symbol,]
UC.down.filtered <- UC.down [UC.down$gene_symbol %in% overlap.genes1$gene_symbol,]
merged.DEG.1 <- cbind (CD.down.filtered, UC.down.filtered)
merged.DEG.1 <- merged.DEG.1 [,c(1,2,3,7)] 
colnames(merged.DEG.1) <- c ("ENTREZ", "SYMBOL", "log2FC_CD", "log2FC_UC")

# overlap genes between CD up and UC down
overlap.genes2 <- CD.up [CD.up$gene_symbol %in% intersect(CD.up$gene_symbol,UC.down$gene_symbol), c(1:3)]
CD.up.filtered <- CD.up [CD.up$gene_symbol %in% overlap.genes2$gene_symbol,]
UC.down.filtered <- UC.down [UC.down$gene_symbol %in% overlap.genes2$gene_symbol,]
merged.DEG.2 <- cbind (CD.up.filtered, UC.down.filtered)
merged.DEG.2 <- merged.DEG.2 [,c(1,2,3,7)] 
colnames(merged.DEG.2) <- c ("ENTREZ", "SYMBOL", "log2FC_CD", "log2FC_UC")

# overlap genes between CD up and UC up
overlap.genes3 <- CD.up [CD.up$gene_symbol %in% intersect(CD.up$gene_symbol,UC.up$gene_symbol), c(1,2)]
CD.up.filtered <- CD.up [CD.up$gene_symbol %in% overlap.genes3$gene_symbol,]
UC.up.filtered  <- UC.up [UC.up$gene_symbol %in% overlap.genes3$gene_symbol,]
merged.DEG.3 <- cbind (CD.up.filtered , UC.up.filtered)
merged.DEG.3 <- merged.DEG.3 [,c(1,2,3,7)] 
colnames(merged.DEG.3) <- c ("ENTREZ", "SYMBOL", "log2FC_CD", "log2FC_UC")

#merge all DEG with corresponding logFC for each disease
DEG.overlapped <- rbind(merged.DEG.1,merged.DEG.2, merged.DEG.3)
if(!dir.exists("output")) dir.create("output")
write.table(DEG.overlapped ,"output/DEG.overlapped_ileum",row.names=FALSE,col.names = TRUE,quote= FALSE, sep = "\t")
##############################################################################################
```

## Finding overlapped genes between diseases on rectum bioopsy location

``` r
######################################FOR RECTUM biopsy location#######################################
# overlap genes between CD down and UC down
overlap.genes1   <- CD.down [CD.down$gene_symbol %in% intersect(CD.down$gene_symbol,UC.down$gene_symbol), c(1:2)]
CD.down.filtered <- CD.down [CD.down$gene_symbol %in% overlap.genes1$gene_symbol,]
UC.down.filtered <- UC.down [UC.down$gene_symbol %in% overlap.genes1$gene_symbol,]
merged.DEG.1 <- cbind (CD.down.filtered, UC.down.filtered)
merged.DEG.1 <- merged.DEG.1 [,c(1,2,3,7)] 
colnames(merged.DEG.1) <- c ("ENTREZ", "SYMBOL", "log2FC_CD", "log2FC_UC")

# overlap genes between CD down and UC up
overlap.genes2 <- CD.down [CD.down$gene_symbol %in% intersect(CD.down$gene_symbol,UC.up$gene_symbol), c(1:3)]
CD.down.filtered <- CD.down [CD.down$gene_symbol %in% overlap.genes2$gene_symbol,]
UC.up.filtered <- UC.up [UC.up$gene_symbol %in% overlap.genes2$gene_symbol,]
merged.DEG.2 <- cbind (CD.down.filtered, UC.up.filtered)
merged.DEG.2 <- merged.DEG.2 [,c(1,2,3,7)] 
colnames(merged.DEG.2) <- c ("ENTREZ", "SYMBOL", "log2FC_CD", "log2FC_UC")

# overlap genes between CD up and UC up
overlap.genes3 <- CD.up [CD.up$gene_symbol %in% intersect(CD.up$gene_symbol,UC.up$gene_symbol), c(1,2)]
CD.up.filtered <- CD.up [CD.up$gene_symbol %in% overlap.genes3$gene_symbol,]
UC.up.filtered  <- UC.up [UC.up$gene_symbol %in% overlap.genes3$gene_symbol,]
merged.DEG.3 <- cbind (CD.up.filtered , UC.up.filtered)
merged.DEG.3 <- merged.DEG.3 [,c(1,2,3,7)] 
colnames(merged.DEG.3) <- c ("ENTREZ", "SYMBOL", "log2FC_CD", "log2FC_UC")

#merge all DEG with corresponding logFC for each disease
DEG.overlapped <- rbind(merged.DEG.1,merged.DEG.2, merged.DEG.3)
if(!dir.exists("output")) dir.create("output")
write.table(DEG.overlapped ,"output/DEG.overlapped_rectum",row.names=FALSE,col.names = TRUE,quote= FALSE, sep = "\t")
##############################################################################################
```

## Last, we create a Jupyter notebook from this script

``` r
#Jupyter Notebook file
if(!"devtools" %in% installed.packages()) BiocManager::install("devtools")
devtools::install_github("mkearney/rmd2jupyter", force=TRUE)
```

    ## 
    ## * checking for file 'C:\Users\dedePC\AppData\Local\Temp\RtmpUrvcyd\remotes40704de21e0e\mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION' ... OK
    ## * preparing 'rmd2jupyter':
    ## * checking DESCRIPTION meta-information ... OK
    ## * checking for LF line-endings in source and make files and shell scripts
    ## * checking for empty or unneeded directories
    ## Omitted 'LazyData' from DESCRIPTION
    ## * building 'rmd2jupyter_0.1.0.tar.gz'
    ## 

``` r
library(devtools)
library(rmd2jupyter)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
rmd2jupyter("extractOverlappedGenes.Rmd")
```
