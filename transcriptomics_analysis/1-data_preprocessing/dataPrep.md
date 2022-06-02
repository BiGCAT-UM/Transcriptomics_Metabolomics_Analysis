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
htxOrj <- read.csv("data/host_tx_counts.tsv",sep = "\t")

#Convert sample names to upper (some of them are in lower case)
colnames(htxOrj)<-toupper(colnames(htxOrj))

#htx count data is filtered based on column names in htxMeta
names.use <- names(htxOrj)[(names(htxOrj) %in% htxMeta$External.ID)]
#filter out htxOrj based on names.use and create a new htxCount
htxCount <- htxOrj[, names.use]
#htxCount data are ordered based on column names to match samples between htxCount and sampleLabels
htxCount <- htxCount[,order(names(htxCount))]

#sample distribution based on biopsy locations
ileum =nrow(htxMeta[htxMeta$biopsy_location=="Ileum",])
rectum = nrow(htxMeta[htxMeta$biopsy_location=="Rectum",])
cat ("Number of samples in ileum:", ileum ,"\nNumber of samples in rectum:",rectum)
```

    ## Number of samples in ileum: 84 
    ## Number of samples in rectum: 91

``` r
#check whether they are in same order
#colnames(htxCount) == htxMeta[,"External.ID"]
#Write all the generated data into the related output files 
write.table(htxCount, "output/htxCount.csv", sep=",",quote=FALSE, row.names = TRUE )
write.table(htxMeta, "output/sampleLabels.csv", sep=",",quote=FALSE,row.names = FALSE, col.names = FALSE)
```

##Print session info:

``` r
##Print session info:
sessionInfo()
```

    ## R version 4.2.0 (2022-04-22)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.6 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=nl_NL.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=nl_NL.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=nl_NL.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=nl_NL.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] dplyr_1.0.9     readxl_1.4.0    rstudioapi_0.13
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] knitr_1.39          magrittr_2.0.3      tidyselect_1.1.2   
    ##  [4] R6_2.5.1            rlang_1.0.2         fastmap_1.1.0      
    ##  [7] fansi_1.0.3         stringr_1.4.0       tools_4.2.0        
    ## [10] xfun_0.31           utf8_1.2.2          DBI_1.1.2          
    ## [13] cli_3.3.0           htmltools_0.5.2     ellipsis_0.3.2     
    ## [16] assertthat_0.2.1    yaml_2.3.5          digest_0.6.29      
    ## [19] tibble_3.1.7        lifecycle_1.0.1     crayon_1.5.1       
    ## [22] purrr_0.3.4         BiocManager_1.30.17 vctrs_0.4.1        
    ## [25] glue_1.6.2          evaluate_0.15       rmarkdown_2.14     
    ## [28] stringi_1.7.6       compiler_4.2.0      pillar_1.7.0       
    ## [31] cellranger_1.1.0    generics_0.1.2      pkgconfig_2.0.3

### Last, we create a Jupyter notebook from this script

``` r
#Jupyter Notebook file
if(!"devtools" %in% installed.packages()) BiocManager::install("devtools")
devtools::install_github("mkearney/rmd2jupyter", force=TRUE)
```

    ## 
    ## * checking for file ‘/tmp/Rtmpv6npOb/remotes3d433037c5b0/mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION’ ... OK
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
