## Introduction

In this section of the workflow, we will obtain the metabolomics data
and apply filtering options, to create a dataset ready for further
statistical and pathway analysis.

## First, we setup the required libraries to get started.

``` r
# check if libraries are already installed > otherwise install it
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager",repos = "http://cran.us.r-project.org")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"dplyr" %in% installed.packages()) BiocManager::install("dplyr")
#Libraries required for markdown documents:
if(!"markdown" %in% installed.packages()){install.packages("markdown")}
if(!"rmarkdown" %in% installed.packages()){install.packages("rmarkdown")}
if(!"stringr" %in% installed.packages()){install.packages("stringr")}
#load libraries
library(rstudioapi)
library(dplyr)
library(stringr)
# set your working environment to the location where your current source file is saved into.
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
```

## Second, we download the required data, read the metadata and filter out not-relevant data.

``` r
#Library to download data from online files:
if(!"downloader" %in% installed.packages()){install.packages("downloader")}
require(downloader)

##Download metadata, extract metabolomics sample IDs, location and disorders.
if(file.exists("data/hmp2_metadata.csv")){print("Metadata already downloaded")}else{
fileUrl <- "https://ibdmdb.org/tunnel/products/HMP2/Metadata/hmp2_metadata.csv?accessType=DOWNLOAD"
require(downloader)
download(fileUrl, "data/hmp2_metadata.csv", mode = "wb")
}
```

    ## [1] "Metadata already downloaded"

``` r
#read metadata file
metaData <- read.csv("data/hmp2_metadata.csv")
#filter out by data type and week number
metaDataMBX <- subset(metaData, metaData$data_type == "metabolomics" )
#we need to have the samples which has same visit number
metaDataMBX<- subset(metaDataMBX, metaDataMBX$visit_num == 4)
#we should match transcriptomics (htx) samples and metabolomics (mbx) samples with participantID
#but samples are given by their externalID in mbx file so we should keep them both
#select columns which will be used
metaDataMBX <- metaDataMBX %>% dplyr::select(External.ID,Participant.ID,diagnosis)
#rename columns of metaDataMBX
colnames(metaDataMBX) <- c("ExternalID","ParticipantID","disease" )

#download and read metabolomics peak intensity data
if(file.exists("data/metabolomics.csv.gz")){print("Metabolomics zipped data already downloaded")}else{
fileUrl <- "https://ibdmdb.org/tunnel/products/HMP2/Metabolites/1723/HMP2_metabolomics.csv.gz?accessType=DOWNLOAD"
download(fileUrl, "data/metabolomics.csv.gz", mode = "wb")
}
```

    ## [1] "Metabolomics zipped data already downloaded"

``` r
#Note: if the URL download does not work, the zipped file is located on GitHub to continue the rest of this script.
if(file.exists("data/metabolomics.csv")){print("Unzipped Metabolomics data already downloaded")}else{
if(!"R.utils" %in% installed.packages()){install.packages("R.utils")}
library(R.utils)
gunzip("data/metabolomics.csv.gz", remove=FALSE)
}
```

    ## [1] "Unzipped Metabolomics data already downloaded"

``` r
mbxData <- read.csv("data/metabolomics.csv")
#delete not used columns
mbxData = subset(mbxData, select = -c(1,2,3,4,7) )
```

## Third, we perform data extraction, and process the data

``` r
### row (metabolite) filtering ###
#delete metabolite or row if it has NA or empty value for hmdbID
mbxData<- mbxData[!(is.na(mbxData$HMDB...Representative.ID.) | mbxData$HMDB...Representative.ID.=="") , ]
#remove rows which has hmdb as "redundant ion"
mbxData<- mbxData[!(mbxData$HMDB...Representative.ID.=="redundant ion") , ]
#remove character (asterisk) in some hmdb column values
mbxData$HMDB...Representative.ID.<- stringr::str_replace(mbxData$HMDB...Representative.ID., '\\*', '')
#Update HMDB IDs to new data structure
mbxData$HMDB...Representative.ID.<- stringr::str_replace(mbxData$HMDB...Representative.ID., 'HMDB', 'HMDB00')
#back up original mbxdata
mbxData.b <- mbxData

### modify mbxData based on sample names given in metaData file (created with the criteria visit_num=4 )###
#filter out mbxData columns (samples) based metaDataMBX externalIDs
names.use <- names(mbxData)[ names(mbxData) %in% metaDataMBX$ExternalID]
#update mbx data with used names
mbxData <- mbxData [ ,names.use]
#order data based on col names
mbxData <- mbxData[ , order(names(mbxData))]

#order metadata based on externalID
metaDataMBX <- metaDataMBX[order(metaDataMBX$ExternalID),]

#add HMDBID and Compound Name column to the mbx data
mbxData <- cbind(mbxData.b$HMDB...Representative.ID., mbxData.b$Metabolite,mbxData)
colnames(mbxData)[1] <- "HMDB.ID"
colnames(mbxData)[2] <- "Compound.Name"

#add disease labels to the mbx data
diseaseLabels <- metaDataMBX$disease
##Add two empty strings to macth with additional column data.
diseaseLabels <- append(diseaseLabels, "NA",after = 0)
diseaseLabels <- append(diseaseLabels, "NA",after = 0)
mbxData <- rbind(diseaseLabels, mbxData)
```

## Fourth, we split up the data for UC and CD, include the control data nonIBD, and save this data to an output folder.

``` r
#write only UC versus nonIBD comparison
mbxDataUC <- mbxData[ ,(mbxData[1, ] == "UC" | mbxData[1, ] == "nonIBD")]
#add hmdb id again
mbxDataUC <- cbind(mbxData[,1:2],mbxDataUC)
colnames(mbxDataUC)[1]="HMBDB.ID"
colnames(mbxDataUC)[2] <- "Compound.Name"
write.table(mbxDataUC, "output/mbxDataUC_nonIBD.csv", sep =",", row.names = FALSE)

#write only CD_healthy comparison
mbxDataCD <- mbxData[ ,(mbxData[1, ] == "CD" | mbxData[1, ] == "nonIBD")]
mbxDataCD <- cbind(mbxData[,1:2],mbxDataCD)
colnames(mbxDataCD)[1]="HMBDB.ID"
colnames(mbxDataCD)[2] <- "Compound.Name"
write.table(mbxDataCD, "output/mbxDataCD_nonIBD.csv", sep =",", row.names = FALSE)
```

### Last, we create a Jupyter notebook and markdown file from this script

``` r
#Jupyter Notebook file
if(!"devtools" %in% installed.packages()) BiocManager::install("devtools")
devtools::install_github("mkearney/rmd2jupyter", force=TRUE)
```

    ## 
    ## ── R CMD build ─────────────────────────────────────────────────────────────────
    ##       ✔  checking for file 'C:\Users\duygu\AppData\Local\Temp\RtmpgVap3m\remotes422049c37151\mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION'
    ##       ─  preparing 'rmd2jupyter':
    ##    checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
    ##       ─  checking for LF line-endings in source and make files and shell scripts
    ##   ─  checking for empty or unneeded directories
    ##    Omitted 'LazyData' from DESCRIPTION
    ##       ─  building 'rmd2jupyter_0.1.0.tar.gz'
    ##      
    ## 

``` r
library(devtools)
library(rmd2jupyter)
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
rmd2jupyter("dataPreprocessing.Rmd")

#markdown_file <- "dataPreprocessing.md"
#if (file.exists(markdown_file)) {
#   unlink(markdown_file, recursive=TRUE)#firs delete the existing one
# }
#Still trows an error, now build with knittr manual selection.
#rmarkdown::render("dataPreprocessing.Rmd", "md_document")

##Clean up R-studio environment
remove(diseaseLabels, fileUrl, names.use, mbxData, mbxData.b, mbxDataCD, mbxDataUC, metaDataMBX, metaData)
```
