## Introduction

In this script, to perform visualization of multi-omics data in
CytoScape we need to prepare data to be imported in CytoScape To do so,
gene and metabolites data are mapped to have ensembelIDs and ChEBI IDS
and corresponding log2FC

## Setup

``` r
# check if libraries are already installed > otherwise install it
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"dplyr" %in% installed.packages()) BiocManager::install("dplyr")
if(!"org.Hs.eg.db" %in% installed.packages()) BiocManager::install("org.Hs.eg.db")
if(!"BridgeDbR" %in% installed.packages()) BiocManager::install("BridgeDbR")
install.packages("rJava",repos = "http://cran.us.r-project.org" )
```

    ## package 'rJava' successfully unpacked and MD5 sums checked
    ## 
    ## The downloaded binary packages are in
    ##  C:\Users\dedePC\AppData\Local\Temp\RtmpamTzyW\downloaded_packages

``` r
 ## See https://www.r-bloggers.com/2018/02/installing-rjava-on-ubuntu/ if you have issues with this package on Ubuntu.

#load libraries
library(rstudioapi)
library(dplyr)
library(org.Hs.eg.db)
library(BridgeDbR)
library(rJava)

# set your working environment to the location where your current source file is saved into.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
```

## Process gene data

``` r
#read data
data <- read.csv("data/table_CD_Ileum_vs_nonIBD_Ileum.tab", sep= "\t")#input file may change depend on the analysis, different data in the data folder 
#filter out unused columns
data <- data [,c(1,3,7)]

#transform Gene Symbols to ensembl IDs
#get ENSEMBL ids
hs <- org.Hs.eg.db
ensemblID <- AnnotationDbi::select(hs, 
            keys = data$X,
            columns = c("ENSEMBL", "SYMBOL"),
            keytype = "SYMBOL")

#filter out double gene symbols, because there are one-to-many relationship between symbol and ensembl IDs
ensemblID <- ensemblID %>% distinct (ensemblID$SYMBOL, .keep_all = TRUE)
#add ensemblIDs to the data
data <- cbind(ensemblID$ENSEMBL,data)
#change column names
colnames(data) <- c("ENSEMBL.ID","SYMBOL", "log2FC_gene", "pvalue_gene")
#filter out genes that has NA value for ensemblID
data<- data %>% tidyr::drop_na(ENSEMBL.ID)
#add type column
data ["data.type"] <- "transcriptomics" 
#write data to file
write.table(data, "output/transcriptomics",row.names = FALSE, quote = FALSE, sep="\t")
```

## Process metabolite data

``` r
#read data
mbxData <- read.csv("data/mbxDataCD_nonIBD.csv")
#filter out unused columns
mbxData <- mbxData [,c(1,3,4)]

##Transform hmdb IDs to CHebi IDs
##Download the Metabolite mapping file:
fileUrl <- "https://ndownloader.figshare.com/files/26001794?accessType=DOWNLOAD"
require(downloader)
download(fileUrl, "data/metabolites.bridge", mode = "wb")

mapper <- BridgeDbR ::loadDatabase('data/metabolites.bridge')

## Obtain the System codes for the databases HMDB (source database of dataset) and ChEBI (intended output database)
code = getOrganismCode("Homo sapiens")
code_mappingFrom <- getSystemCode("HMDB")
code_mappingTo   <- getSystemCode("ChEBI")

## Create a data frame with the mappings and the correct SystemCode
input = data.frame(
    source = rep(code_mappingFrom, length(mbxData$X)),
    identifier = mbxData$X)
#Obtain all mappings from HMDB to ChEBI
MultiMappings = BridgeDbR::maps(mapper, input, code_mappingTo)
#remove all rows in the mapped data which do not include the prefix "CHEBI"
MultiMappings <- MultiMappings %>% filter(grepl("CHEBI",mapping, fixed = TRUE))
#filter out double identifiers because there are one-to-many relationship between hmdb and chebi IDs
MultiMappings <- MultiMappings %>% distinct (MultiMappings$identifier, .keep_all = TRUE)
MultiMappings <- MultiMappings [,c(2,4)]

merged.data<- merge(MultiMappings, mbxData,by.x="identifier", by.y="X",sort = TRUE, all.x = TRUE, all.y = TRUE)
#filter out metabolites that has NA value for CHEBI
merged.data<- merged.data %>% tidyr::drop_na(mapping)
#change column names
colnames(merged.data) <- c("HMDBID","CHEBI", "log2FC_met", "pvalue_met")
#add type column 
merged.data ["data.type"] <- "metabolomics" 
#write data to the file
write.table(merged.data, "output/metabolomics",row.names = FALSE, quote = FALSE, sep="\t")
```
