## Introduction

In this script, filtering options will be applied for transcriptomics
data to be prepared for the analysis

## Setup

``` r
# check if libraries are already installed > otherwise install it
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager",repos = "http://cran.us.r-project.org")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"readxl" %in% installed.packages()) BiocManager::install("readxl")
if(!"dplyr" %in% installed.packages()) BiocManager::install("dplyr")

#load libraries
library(rstudioapi)
library(readxl)
library(dplyr)

# set working environment to the location where current source file is saved into.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
```

## Read and filter out metadata

``` r
##Download metadata, extract transcriptomics sample IDs, location and disorders.
if(file.exists("data/hmp2_metadata.csv")){print("Metadata already downloaded")}else{
fileUrl <- "https://ibdmdb.org/tunnel/products/HMP2/Metadata/hmp2_metadata.csv?accessType=DOWNLOAD"
require(downloader)
download(fileUrl, "data/hmp2_metadata.csv", mode = "wb")
}
```

    ## [1] "Metadata already downloaded"

``` r
#Read metadata
htxMeta <- read.csv("data/hmp2_metadata.csv")
#filter out by data type as host-transcriptomics
htxMeta <- htxMeta  %>% filter(htxMeta$data_type == "host_transcriptomics")

#filter out data by biopsy location, include CD, UC and nonIBD samples from ileum and rectum location 
htxMeta <-htxMeta  %>% filter(  (htxMeta$diagnosis == "CD" & htxMeta$biopsy_location=="Ileum") 
                              | (htxMeta$diagnosis == "CD" & htxMeta$biopsy_location=="Rectum")
                              | (htxMeta$diagnosis == "UC" & htxMeta$biopsy_location=="Ileum") 
                              | (htxMeta$diagnosis == "UC" & htxMeta$biopsy_location=="Rectum") 
                              | (htxMeta$diagnosis == "nonIBD" & htxMeta$biopsy_location=="Rectum") 
                              | (htxMeta$diagnosis == "nonIBD" & htxMeta$biopsy_location=="Ileum") 
)
#filter out samples by visit_num=1
htxMeta <-htxMeta  %>% filter(htxMeta$visit_num == "1")

#filter out unused columns
htxMeta <- htxMeta %>% dplyr::select(External.ID,Participant.ID,biopsy_location,diagnosis)

#Order htxMeta data based on external ID to match samples with htx count correctly
htxMeta<- htxMeta[order(htxMeta$External.ID),]#order htxMeta by external ID
```

# Filter out host transcriptomics (htx) count data based on sample names obtained from htx meta data, save data in output file.

``` r
#transcript count (htx count) original file is read
htxOrj <- read.csv("host_tx_counts.tsv",sep = "\t")

#Convert sample names to upper (some of them are in lower case)
colnames(htxOrj)<-toupper(colnames(htxOrj))

#htx count data are filtered based on col names in htxMeta
names.use <- names(htxOrj)[(names(htxOrj) %in% htxMeta$External.ID)]
#filter out htxOrj based on names.use and create htxCount
htxCount <- htxOrj[, names.use]
#htxCount data are ordered based on column names to match samples between htxCount and sampleLabels
htxCount <- htxCount[,order(names(htxCount))]
#check whether they are in same order
#colnames(htxCount) == htxMeta[,"External.ID"]
#Write all the generated data into the related output files 
write.table(htxCount, "output/htxCount.csv", sep="\t",quote=FALSE, row.names = TRUE )
write.table(htxMeta, "output/sampleLabels.csv", sep="\t",quote=FALSE,row.names = FALSE, col.names = FALSE)
```

### Last, we create a Jupyter notebook from this script

``` r
#Jupyter Notebook file
if(!"devtools" %in% installed.packages()) BiocManager::install("devtools")
devtools::install_github("mkearney/rmd2jupyter", force=TRUE)
```

    ## 
    ## * checking for file ‘/tmp/RtmpKwm5dh/remotes7ca0109cad0/mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION’ ... OK
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
rmd2jupyter("dataPrep.Rmd")
```
