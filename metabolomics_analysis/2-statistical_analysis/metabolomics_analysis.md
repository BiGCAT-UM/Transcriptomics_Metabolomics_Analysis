## Introduction

In this workflow, we will apply statis one factoral analysis on
metabolomics data using MetaboAnalyst software R packages.
MetaboAnalystR 3 contains the R functions and libraries underlying the
popular MetaboAnalyst web server, including metabolomic data analysis,
visualization, and functional interpretation. The package is
synchronized with the MetaboAnalyst web server. After installing and
loading the package, users will be able to reproduce the same results
from their local computers using the corresponding R command history
downloaded from MetaboAnalyst web site, thereby achieving maximum
flexibility and reproducibility.
source:<https://www.metaboanalyst.ca/MetaboAnalyst/docs/RTutorial.xhtml>

## Setup

Installing and loading required libraries

``` r
# check if libraries are already installed > otherwise install it
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"MetaboAnalystR" %in% installed.packages()) BiocManager::install("MetaboAnalystR")

#load installed librariers
library(BiocManager)
library(rstudioapi)
library(MetaboAnalystR)
```

# Data import and preprocessing

``` r
# set your working environment to the location where your current source file is saved into.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
work_DIR <- getwd()

#initialize data object
mSet<-InitDataObjects("pktable", "stat", FALSE)
```

    ## Starting Rserve...
    ##  "C:\Users\dedePC\DOCUME~1\R\WIN-LI~1\4.1\Rserve\libs\x64\Rserve.exe" --no-save 
    ## [1] "MetaboAnalyst R objects initialized ..."

``` r
#read data
mSet<-Read.TextData(mSet, "data/mbxDataCD_nonIBD.csv", "colu", "disc");
mSet<-SanityCheckData(mSet)
```

    ##  [1] "Successfully passed sanity check!"                                                                                
    ##  [2] "Samples are not paired."                                                                                          
    ##  [3] "2 groups were detected in samples."                                                                               
    ##  [4] "Only English letters, numbers, underscore, hyphen and forward slash (/) are allowed."                             
    ##  [5] "<font color=\"orange\">Other special characters or punctuations (if any) will be stripped off.</font>"            
    ##  [6] "All data values are numeric."                                                                                     
    ##  [7] "<font color=\"red\"> 3 features with a constant or single value across samples were found and deleted.</font>"    
    ##  [8] "A total of 1721 (7.6%) missing values were detected."                                                             
    ##  [9] "<u>By default, missing values will be replaced by 1/5 of min positive values of their corresponding variables</u>"
    ## [10] "Click the <b>Proceed</b> button if you accept the default practice;"                                              
    ## [11] "Or click the <b>Missing Values</b> button to use other methods."

``` r
#missing value removal by filtering out features has more than 50% missing
mSet<-RemoveMissingPercent(mSet, percent=0.5)
#set the rest features as column mean
mSet<-ImputeMissingVar(mSet, method="mean")
mSet<-SanityCheckData(mSet)
```

    ##  [1] "Successfully passed sanity check!"                                                                                
    ##  [2] "Samples are not paired."                                                                                          
    ##  [3] "2 groups were detected in samples."                                                                               
    ##  [4] "Only English letters, numbers, underscore, hyphen and forward slash (/) are allowed."                             
    ##  [5] "<font color=\"orange\">Other special characters or punctuations (if any) will be stripped off.</font>"            
    ##  [6] "All data values are numeric."                                                                                     
    ##  [7] "<font color=\"red\"> 3 features with a constant or single value across samples were found and deleted.</font>"    
    ##  [8] "A total of 1721 (7.6%) missing values were detected."                                                             
    ##  [9] "<u>By default, missing values will be replaced by 1/5 of min positive values of their corresponding variables</u>"
    ## [10] "Click the <b>Proceed</b> button if you accept the default practice;"                                              
    ## [11] "Or click the <b>Missing Values</b> button to use other methods."

``` r
#no filter applied
mSet<-FilterVariable(mSet, "none", "F", 25)
```

    ## [1] " No filtering was applied"
    ## [1] "  No filtering was applied"

``` r
# Saving output
if(dir.exists("output"))#if the output folder already exist
  unlink("output", recursive=TRUE)#first delete the existing one
dir.create("output")#create a new output folder
#setwd(paste0(work_DIR,"/output"))

#normalization
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "CrNorm", "NULL", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
```

## Volcano plot

``` r
mSet<-Volcano.Anal(mSet, FALSE, 1.0, 0, F, 1.0, TRUE, "raw")
```

    ## [1] "Performing regular t-tests ...."
    ## [1] "A total of 436 significant features were found."

``` r
mSet<-PlotVolcano(mSet, "volcano_1_",1, 0, "png", 72, width=NA)
```
