## Introduction

Before we can visualize the multi-omics data in CytoScape, we need to
prepare data to be imported according to the unified database
identifiers available. Therefor, metabolite data is mapped to ChEBI IDS.

## Setup

``` r
# check if libraries are already installed > otherwise install it
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"dplyr" %in% installed.packages()) BiocManager::install("dplyr")
if(!"BridgeDbR" %in% installed.packages()) BiocManager::install("BridgeDbR")
install.packages("rJava",repos = "http://cran.us.r-project.org" )
 ## See https://www.r-bloggers.com/2018/02/installing-rjava-on-ubuntu/ if you have issues with this package on Ubuntu.

#load libraries
library(rstudioapi)
library(dplyr)
library(BridgeDbR)
library(rJava)

# set your working environment to the location where your current source file is saved into.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
```

## Process metabolite data

``` r
setwd('..')

#Obtain data from step 8
mSet_CD <- read.csv("8-significantly_changed_metabolites_analysis/output/mbxData_CD.csv", na.strings=c("", "NA"))
mSet_UC <- read.csv("8-significantly_changed_metabolites_analysis/output/mbxData_UC.csv", na.strings=c("", "NA"))
#filter out unused columns
mSet_CD <- mSet_CD [,c(1:4)]
mSet_UC <- mSet_UC [,c(1:4)]

# Set Working Directory back to current folder
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
work_DIR <- getwd()

##Retain all metabolite IDs from CD and UC for identifier mapping:
mSet_total <- unique(rbind(mSet_CD[,c(1,2)], mSet_UC[,c(1,2)]))

#Download the Metabolite mapping file (if it doesn't exist locally yet):
checkfile <- paste0(getwd(), '/' ,"data/metabolites.bridge")
if (!file.exists(checkfile)) {
  download.file("https://figshare.com/ndownloader/files/36197283", checkfile)
}
##Load the metabolite mapping file:
mapper <- BridgeDbR ::loadDatabase(checkfile)

## Obtain the System codes for the databases HMDB (source database of dataset) and ChEBI (intended output database)
code <- getOrganismCode("Homo sapiens")
code_mappingFrom <- getSystemCode("HMDB")
code_mappingTo   <- getSystemCode("ChEBI")

## Create a data frame with the mappings and the correct SystemCode
input = data.frame(
    source = rep(code_mappingFrom, length(mSet_total$HMDB_ID)),
    identifier = mSet_total$HMDB_ID)
#Obtain all mappings from HMDB to ChEBI
MultiMappings = BridgeDbR::maps(mapper, input, code_mappingTo)
#remove all rows in the mapped data which do not include the prefix "CHEBI"
MultiMappings <- MultiMappings %>% filter(grepl("CHEBI",mapping, fixed = TRUE))
#filter out double identifiers because there are one-to-many relationship between hmdb and chebi IDs
MultiMappings <- MultiMappings %>% distinct (MultiMappings$identifier, .keep_all = TRUE)
MultiMappings <- MultiMappings [,c(2,4)]

merged.data_CD<- merge(MultiMappings, mSet_CD,by.x="identifier", by.y="HMDB_ID",sort = TRUE, all.x = TRUE, all.y = TRUE)
merged.data_UC<- merge(MultiMappings, mSet_UC,by.x="identifier", by.y="HMDB_ID",sort = TRUE, all.x = TRUE, all.y = TRUE)
#filter out metabolites that has NA value for CHEBI
merged.data_CD<- merged.data_CD %>% tidyr::drop_na(mapping)
merged.data_UC<- merged.data_UC %>% tidyr::drop_na(mapping)
#filter out metabolites that has NA value for foldchange_disorder
merged.data_CD<- merged.data_CD %>% tidyr::drop_na(foldchange_disorder)
merged.data_UC<- merged.data_UC %>% tidyr::drop_na(foldchange_disorder)
#change column names
colnames(merged.data_CD) <- c("HMDBID","CHEBI", "label", "log2FC_met", "pvalue_met")
colnames(merged.data_UC) <- c("HMDBID","CHEBI", "label", "log2FC_met", "pvalue_met")
```

## Export the mapped metabolomics data:

``` r
##Save the data file
write.table(merged.data_CD, 'output/mbxPWdata_CD.tsv', sep ="\t", row.names = FALSE)
write.table(merged.data_UC, 'output/mbxPWdata_UC.tsv', sep ="\t", row.names = FALSE)
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
    ##  [1] LC_CTYPE=en_US.UTF-8          LC_NUMERIC=C                 
    ##  [3] LC_TIME=nl_NL.UTF-8           LC_COLLATE=en_US.UTF-8       
    ##  [5] LC_MONETARY=nl_NL.UTF-8       LC_MESSAGES=en_US.UTF-8      
    ##  [7] LC_PAPER=nl_NL.UTF-8          LC_NAME=nl_NL.UTF-8          
    ##  [9] LC_ADDRESS=nl_NL.UTF-8        LC_TELEPHONE=nl_NL.UTF-8     
    ## [11] LC_MEASUREMENT=nl_NL.UTF-8    LC_IDENTIFICATION=nl_NL.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] BridgeDbR_2.7.3 rJava_1.0-6     dplyr_1.0.9     rstudioapi_0.13
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] knitr_1.39       magrittr_2.0.3   tidyselect_1.1.2 R6_2.5.1        
    ##  [5] rlang_1.0.2      fastmap_1.1.0    fansi_1.0.3      stringr_1.4.0   
    ##  [9] tools_4.2.0      xfun_0.31        utf8_1.2.2       DBI_1.1.2       
    ## [13] cli_3.3.0        htmltools_0.5.2  ellipsis_0.3.2   assertthat_0.2.1
    ## [17] yaml_2.3.5       digest_0.6.29    tibble_3.1.7     lifecycle_1.0.1 
    ## [21] crayon_1.5.1     tidyr_1.2.0      purrr_0.3.4      vctrs_0.4.1     
    ## [25] curl_4.3.2       glue_1.6.2       evaluate_0.15    rmarkdown_2.14  
    ## [29] stringi_1.7.6    compiler_4.2.0   pillar_1.7.0     generics_0.1.2  
    ## [33] pkgconfig_2.0.3

### Last, we create a Jupyter notebook and markdown file from this script

``` r
#Jupyter Notebook file
if(!"devtools" %in% installed.packages()) BiocManager::install("devtools")
devtools::install_github("mkearney/rmd2jupyter", force=TRUE)
```

    ## 
    ## * checking for file ‘/tmp/RtmpKKUayn/remotes2422526e7276/mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION’ ... OK
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
rmd2jupyter("metabolomics_identifier_mapping.Rmd")
```
