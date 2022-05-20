## Introduction

In this workflow, we will perform pathway enrichment analysis to
differential expressed genes for both diseases on two biopsy locations
ileum and rectum. We will be performing pathway enrichment analysis on a
differential gene expression dataset. The dataset compares the
expression of transcripts in inflammatory bowel disease (CD and UC)
biopsy locations (ileum and rectum) versus healthy people. Differential
expression analysis has already been performed (see step 2 of this
workflow), generating log2foldchange and (adjusted) p-values data for
each gene. \## R environment setup

``` r
# check if libraries are already installed > otherwise install it
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"org.Hs.eg.db" %in% installed.packages()) BiocManager::install("org.Hs.eg.db")  
if(!"AnnotationDbi" %in% installed.packages()) BiocManager::install("AnnotationDbi")
if(!"dplyr" %in% installed.packages()) BiocManager::install("dplyr") #for using %>% function
if(!"rWikiPathways" %in% installed.packages()) BiocManager::install("rWikiPathways")
if(!"clusterProfiler" %in% installed.packages()) BiocManager::install("clusterProfiler") 
if(!"ggplot2" %in% installed.packages()) BiocManager::install("ggplot2") 
if(!"pheatmap" %in% installed.packages()) BiocManager::install("pheatmap")
if(!"RColorBrewer" %in% installed.packages()) BiocManager::install("RColorBrewer")

#loading installed libraries
library(rstudioapi) # interface for interacting with RStudio IDE with R code.
library(org.Hs.eg.db) #This is the organism annotation package ("org") for Homo sapiens ("Hs"), organized as an AnnotationDbi   package ("db"), using Entrez Gene IDs ("eg") as primary key.
library(AnnotationDbi) # for connecting and querying annotation databases
library(dplyr)#for data manipulation
library(rWikiPathways) # for programmatic access to WikiPathways content
library(clusterProfiler) # for implementing methods to analyze and visualize functional profiles of genomic data
library(ggplot2) # for creating graphics
library (pheatmap) # for creating heatmap
library (RColorBrewer) # for managing colors 

# set your working environment to the location where your current source file is saved into.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
```

## Importing dataset

The data will be read for the disease on two biopsy locations

``` r
##TODO: add script to include UC or CD disorder data, depending on User input.

##Obtain data from step 2:
setwd('..')
work_DIR <- getwd()
#we have two datasets from different biopsy locations
dataset1 <- read.delim("2-differential_gene_expression_analysis/statsmodel/table_UC_Ileum_vs_nonIBD_Ileum.tab")
dataset2 <- read.delim("2-differential_gene_expression_analysis/statsmodel/table_UC_Rectum_vs_nonIBD_Rectum.tab")

# Set Working Directory back to current folder
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
work_DIR <- getwd()

#filter out  unused columns, we select geneSymbol, log2FC and pvalue
dataset1<- subset( dataset1, select = c(1,3,7))
dataset2<- subset( dataset2, select = c(1,3,7))

#merge two dataset into one data 
dataset <- merge(dataset1, dataset2,by.x="X", by.y="X",sort = TRUE, all.x = TRUE, all.y = TRUE)
#change column names
colnames(dataset) <- c("GeneSymbol","log2FC_ileum","pvalue_ileum","log2FC_rectum","pvalue_rectum")
```

## Converting gene symbols to the corresponding entrez IDs

``` r
#converting gene symbols to entrez ID since the enrichR function needs entrez IDs for each gene symbol
hs <- org.Hs.eg.db #This object is a simple mapping of Entrez Gene identifier
entrezID <- AnnotationDbi::select(hs, keys = dataset$GeneSymbol,
            columns = c("ENTREZID", "SYMBOL"),
            keytype = "SYMBOL")
#filter out double gene symbols
entrezID <- entrezID %>% distinct (entrezID$SYMBOL, .keep_all = TRUE)
# add entrezIDs for each gene symbol in the dataset
dataset <- cbind(entrezID$ENTREZID,dataset)
#change column name
colnames(dataset)[1] = "ENTREZ.ID"
#filter out genes that has NA value for entrezID
dataset<- dataset %>% tidyr::drop_na(ENTREZ.ID)
```

## Getting differentially expressed genes including up-regulated and down-regulated genes using the criteria

``` r
#we will use selection criteria as Fold change=1.5,log2FC=0.58 and p.value < 0.05
#for ileum location
up.genes.ileum   <- dataset[dataset$log2FC_ileum >= 0.58 & dataset$pvalue_ileum < 0.05, 1] 
down.genes.ileum <- dataset[dataset$log2FC_ileum <= -0.58 & dataset$pvalue_ileum < 0.05, 1] 
deg.ileum <- unique(dataset[!is.na(dataset$pvalue_ileum) & dataset$pvalue_ileum < 0.05 
                            & abs(dataset$log2FC_ileum) > 0.58,])
if(!dir.exists("results")) dir.create("results")
write.table(deg.ileum, file="results/DEGs_UC_ileum",sep="\t", quote=FALSE, row.names = FALSE)

#for rectum location
up.genes.rectum   <- dataset[dataset$log2FC_rectum >= 0.58 & dataset$pvalue_rectum < 0.05, 1] 
down.genes.rectum <- dataset[dataset$log2FC_rectum <= -0.58 & dataset$pvalue_rectum < 0.05, 1] 
deg.rectum <- unique(dataset[!is.na(dataset$pvalue_rectum) & dataset$pvalue_rectum < 0.05 & abs(dataset$log2FC_rectum) > 0.58,])
write.table(deg.rectum, file="results/DEGs_UC_rectum",sep="\t", quote=FALSE, row.names = FALSE)

#background genes to be used in enrichment analysis
bkgd.genes <- unique(dataset[,c(1,2)])
```

## Pathway Enrichment Analysis

In this section, we will perform pathway enrichment analysis. The
clusterProfiler R-package is used to perform over-representation
analysis (ORA). The function can be easily replaced to use other
enrichment methods (GSEA / rSEA / etc).

``` r
#below code should be performed first to handle the ssl certificate error while uploading pathways 
options(RCurlOptions = list(cainfo = paste0( tempdir() , "/cacert.pem" ), ssl.verifypeer = FALSE))
#downloading latest pathway gmt files for human 
wp.hs.gmt <- rWikiPathways::downloadPathwayArchive(organism="Homo sapiens", format = "gmt")

#Now that we have got the latest GMT file for human pathways, 
#we can process it to generate the two dataframes we need for enricher
wp2gene   <- rWikiPathways::readPathwayGMT(wp.hs.gmt)
wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid,name) #TERM2NAME

##################ILEUM location#######################
# The clusterProfiler R-package is used to perform over-representation analysis (ORA)
# The function can be easily replaced to use other enrichment methods (GSEA / rSEA / etc). 
ewp.ileum <- clusterProfiler::enricher(
  deg.ileum$ENTREZ.ID,# a vector of gene IDs
  universe = bkgd.genes$ENTREZ.ID,
  pAdjustMethod = "fdr",#you can change this to hochberg, bonferronni or none etc.
  pvalueCutoff = 1, #adjusted pvalue cutoff on enrichment tests to report, we set it a wider criteria then we will filter out
  #results based on padjust and qvalue in next section which is enrichment result visualization
  qvalueCutoff = 1, #qvalue cutoff on enrichment tests to report as significant, 
                       # multiple hypothesis testing
  TERM2GENE = wpid2gene, #user input annotation of TERM TO GENE mapping
  TERM2NAME = wpid2name) #user input of TERM TO NAME mapping

ewp.ileum.res <- as.data.frame(ewp.ileum) 

#interpretation of the output: 
#     BgRatio   = (number of genes measured in the current pathway) / (number of genes measured in all pathways)
#     geneRatio = (number of DEGs in the current pathway) / (total number of DEGs in all pathways)

##TODO: print some next before numbers, so output is more user friendly.
# number of genes measured in all pathways
length(ewp.ileum@universe)
```

    ## [1] 6030

``` r
# number of DEGs in all pathways
length(deg.ileum$ENTREZ.ID[deg.ileum$ENTREZ.ID %in% unique(wp2gene$gene)])
```

    ## [1] 185

``` r
#number of enriched pathways
num.pathways.ileum <- dim(ewp.ileum.res)[1]

#exporting results to the file
write.table(ewp.ileum.res, file=paste("results/enrichResults_","UC_ileum",sep = ""),
            sep = "\t" ,quote = FALSE, row.names = FALSE)

##################RECTUM location#######################
ewp.rectum <- clusterProfiler::enricher(
  deg.rectum$ENTREZ.ID,
  universe = bkgd.genes$ENTREZ.ID,
  pAdjustMethod = "fdr",#you can change it as BH, bonferronni etc.
  pvalueCutoff = 1, #padjust cutoff
  qvalueCutoff = 1, #q value cutoff 
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)
ewp.rectum.res <- as.data.frame(ewp.rectum) 

# number of genes measured in all pathways
length(ewp.rectum@universe)
```

    ## [1] 6030

``` r
# number of DEGs in all pathways
length(deg.rectum$ENTREZ.ID[deg.rectum$ENTREZ.ID %in% unique(wp2gene$gene)])
```

    ## [1] 1485

``` r
#number of enriched pathways
num.pathways.rectum <- dim(ewp.rectum.res)[1]

#exporting results to the file
write.table(ewp.rectum.res, file=paste("results/enrichResults_","UC_rectum",sep = ""),
            sep = "\t" ,quote = FALSE, row.names = FALSE)
```

### Last, we create a Jupyter notebook from this script

``` r
#Jupyter Notebook file
if(!"devtools" %in% installed.packages()) BiocManager::install("devtools")
devtools::install_github("mkearney/rmd2jupyter", force=TRUE)
```

    ## 
    ## * checking for file ‘/tmp/Rtmpqu2h1G/remotes17f352a4aa1f/mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION’ ... OK
    ## * preparing ‘rmd2jupyter’:
    ## * checking DESCRIPTION meta-information ... OK
    ## * checking for LF line-endings in source and make files and shell scripts
    ## * checking for empty or unneeded directories
    ## Omitted ‘LazyData’ from DESCRIPTION
    ## * building ‘rmd2jupyter_0.1.0.tar.gz’

``` r
library(devtools)
library(rmd2jupyter)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
rmd2jupyter("pathway_analysis.Rmd")
```
