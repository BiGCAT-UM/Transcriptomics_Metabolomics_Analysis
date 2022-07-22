## Introduction

In this script, visualization of the enriched pathway which include both
altered genes and metabolites is performed.

## Setup

``` r
# check if libraries are already installed > otherwise install it
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"RCy3" %in% installed.packages()) BiocManager::install("RCy3")
if(!"rWikiPathways" %in% installed.packages()) BiocManager::install("rWikiPathways")
if(!"RColorBrewer" %in% installed.packages()) BiocManager::install("RColorBrewer")
if(!"dplyr" %in% installed.packages()) BiocManager::install("dplyr")
#load libraries
library(rstudioapi)
library(RCy3)#connect cytoscape via R
library(rWikiPathways)#to get pathways from WikiPathways
library(RColorBrewer)#to manage colors with R
library(dplyr)
# set your working environment to the location where your current source file is saved into.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
```

##Obtain transcriptomics and metabolomics data

``` r
setwd('../..')
work_DIR <- getwd()

#Set location to download data for transcriptomics:
filelocation_t <- paste0(work_DIR, "/transcriptomics_analysis/3-identifier_mapping/output/")
#Obtain data from step 3
tSet_CD <- read.delim(paste0(filelocation_t, 'IDMapping_CD.tsv'), sep = "\t", na.strings=c("", "NA"))
tSet_UC <- read.delim(paste0(filelocation_t, 'IDMapping_UC.tsv'), sep = "\t", na.strings=c("", "NA"))
##Select the corresponding location to be visualized:
location_transcriptomics <- "ileum" ##Options: ileum or rectum
#filter out unused columns
if(location_transcriptomics == "ileum"){
tSet_CD <- na.omit(tSet_CD [,c(1,4:5)])
tSet_UC <- na.omit(tSet_UC [,c(1,4:5)])}else if(location_transcriptomics == "rectum"){
tSet_CD <- na.omit(tSet_CD [,c(1,6:7)])
tSet_UC <- na.omit(tSet_UC [,c(1,6:7)])}else{print("Location for transcriptomics data not recognised.")}

#Rename columns for merger later
colnames(tSet_CD) <- c('ID','log2FC','pvalues')
colnames(tSet_UC) <- c('ID','log2FC','pvalues')

#Set location to download data for metabolomics:
filelocation_m <- paste0(work_DIR, "/metabolomics_analysis/10-identifier_mapping/output/")
#Obtain data from step 10
mSet_CD <- read.csv(paste0(filelocation_m, 'mbx_mapped_data_CD.tsv'), sep = "\t", na.strings=c("", "NA"))
mSet_UC <- read.csv(paste0(filelocation_m, 'mbx_mapped_data_UC.tsv'), sep = "\t", na.strings=c("", "NA"))
#filter out unused columns
mSet_CD <- na.omit(mSet_CD [,c(2,4:5)])
mSet_UC <- na.omit(mSet_UC [,c(2,4:5)])

#Rename columns for merger later
colnames(mSet_CD) <- c('ID','log2FC','pvalues')
colnames(mSet_UC) <- c('ID','log2FC','pvalues')

# Set Working Directory back to current folder
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
work_DIR <- getwd()
```

##Combine both dataframes

``` r
combined.data_CD <- rbind(tSet_CD, mSet_CD)
combined.data_UC <- rbind(tSet_UC, mSet_UC)
##Select disorder to visualize later on:
disorder <- "CD" ##Options are: CD,UC
if(disorder == "CD"){combined.data <- combined.data_CD}else if(disorder == "UC"){combined.data <- combined.data_UC}else{print("Disorder not recognized.")}
```

## Import pathway

``` r
#make sure to launch cytoscape, if you get CyREST error you need to relaunch cytoscape
cytoscapePing()
#close all opened session before starting
closeSession(FALSE)
#Set up WikiPathways app in Cytoscape, v.3.3.10
if("WikiPathways" %in% commandsHelp("")) print("Success: the WikiPathways app is installed") else print("Warning: WikiPathways app is not installed. Please install the WikiPathways app before proceeding.")
```

    ## [1] "Available namespaces:"
    ## [1] "Warning: WikiPathways app is not installed. Please install the WikiPathways app before proceeding."

``` r
if(!"WikiPathways" %in% commandsHelp("")) installApp("WikiPathways")
```

    ## [1] "Available namespaces:"

``` r
#pathway IDs to be visualized
pathway.id <- "WP4726"# Sphingolipid metabolism: integrated pathway is a relevant and significantly altered pathways regarding metabolomics data.
#import pathways as pathway in cytoscape
RCy3::commandsRun(paste0('wikipathways import-as-pathway id=',pathway.id )) 
```

## Data upload

``` r
#get node table from imported pathway in cytoscape
ID.cols <- getTableColumns(table ="node", columns = c("XrefId","Ensembl", "ChEBI"))
#filter out rows which contain NA value for columns Ensembl and ChEBI
ID.cols <- ID.cols[!with(ID.cols, is.na(Ensembl) & is.na(ChEBI)),]
#if a row value in the Ensembl column is NA then replace it with ChEBI  
ID.cols$Ensembl <- ifelse(is.na(ID.cols$Ensembl), ID.cols$ChEBI, ID.cols$Ensembl)
#use the only one column contains both Ensembl and ChEBI identifiers
ID.cols <- data.frame(ID.cols[,c(1,2)])
#change column name
colnames(ID.cols)[2] <- "omics.ID"

#merge two data frames for adding xrefid to the combined data frame
data <- merge(combined.data, ID.cols, by.x = "ID", by.y = "omics.ID" )
#remove duplicate rows
data <- data %>% distinct(ID, .keep_all = TRUE)
colnames(data)[1] <- "omics.ID"
#load data to the imported pathway in cytoscape by key column as XrefId
loadTableData(table = "node", data = data, data.key.column = "XrefId", table.key.column = "XrefId")
```

    ## [1] "Success: Data loaded in defaultnode table"

## Visualization options

``` r
#new visual style is created
RCy3::copyVisualStyle("default","pathwayStyle")
#set new style as the current style
RCy3::setVisualStyle("pathwayStyle")
```

    ##                 message 
    ## "Visual Style applied."

``` r
#set node dimensions as fixed sizes
RCy3::lockNodeDimensions(TRUE, style.name="pathwayStyle")

#node shape mapping
RCy3::setNodeShapeMapping('Type',c('GeneProduct','Protein', 'Metabolite'),c('ELLIPSE','ELLIPSE','RECTANGLE'), style.name="pathwayStyle")
```

    ## NULL

``` r
#change node height
RCy3::setNodeHeightMapping('Type',c('GeneProduct','Protein', 'Metabolite'), c(23,23,25), mapping.type = "d", style.name = "pathwayStyle")
```

    ## NULL

``` r
#change node width
RCy3::setNodeWidthMapping('Type',c('GeneProduct','Protein', 'Metabolite'), c(60,60,100), mapping.type = "d", style.name = "pathwayStyle")
```

    ## NULL

``` r
#set node color based on log2FC for both genes and metabolites
node.colors <- c(rev(brewer.pal(3, "RdBu")))
setNodeColorMapping("log2FC", c(-1,0,1), node.colors, default.color = "#D3D3D3", style.name = "pathwayStyle")
```

    ## NULL

``` r
#Set node border width and color based on p-value
#First we need to get all p-values from node table
pvalues <- getTableColumns(table = 'node', columns = 'pvalues')
pvalues <- na.omit(pvalues)
#Create a range for all sign. p-values, and one for all not significant.
significant_pvalues <- pvalues[(pvalues < 0.05)]
not.significant_pvalues <- pvalues[(pvalues >= 0.05)]
significant_pvalues.colors <- rep("#2e9d1d", length(significant_pvalues))
not.significant_pvalues.colors <- rep("#FFFFFF", length(not.significant_pvalues))
RCy3::setNodeBorderWidthMapping('pvalues', table.column.values = NULL , c(6,6) , mapping.type = "c", style.name = "pathwayStyle")
```

    ## NULL

``` r
RCy3::setNodeBorderColorMapping('pvalues', c(significant_pvalues,not.significant_pvalues), c(significant_pvalues.colors, not.significant_pvalues.colors), default.color = "#AAAAAA", mapping.type = "d", style.name = "pathwayStyle")
```

    ## NULL

``` r
##Update relevant interactions to directional ones:
RCy3::setEdgeTargetArrowShapeMapping(table.column = 'EndArrow', c('mim-conversion', 'Arrow', 'mim-catalysis'), c('DELTA', 'DELTA', 'OPEN_CIRCLE'), style.name = "pathwayStyle")
```

    ## NULL

``` r
#Save output 
filename_multiomics <- paste0("output/", pathway.id, "_", disorder, "_location_", location_transcriptomics,"_visualization.png")
png.file <- file.path(getwd(), filename_multiomics)
exportImage(png.file, 'PNG', zoom = 500)
```

    ##                                                                                                                                                                     file 
    ## "/home/deniseslenter/Documents/GitHub/Transcriptomics_Metabolomics_Analysis/visualization_multiomics/12-visualization/output/WP4726_CD_location_ileum_visualization.png"

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
    ## [1] dplyr_1.0.9          RColorBrewer_1.1-3   rWikiPathways_1.16.0
    ## [4] RCy3_2.16.0          rstudioapi_0.13     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] pbdZMQ_0.3-7        tidyselect_1.1.2    xfun_0.31          
    ##  [4] repr_1.1.4          purrr_0.3.4         vctrs_0.4.1        
    ##  [7] generics_0.1.3      htmltools_0.5.3     stats4_4.2.0       
    ## [10] yaml_2.3.5          base64enc_0.1-3     utf8_1.2.2         
    ## [13] XML_3.99-0.9        rlang_1.0.4         pillar_1.8.0       
    ## [16] DBI_1.1.3           glue_1.6.2          BiocGenerics_0.42.0
    ## [19] uuid_1.1-0          lifecycle_1.0.1     stringr_1.4.0      
    ## [22] evaluate_0.15       uchardet_1.1.0      knitr_1.39         
    ## [25] fastmap_1.1.0       curl_4.3.2          fansi_1.0.3        
    ## [28] IRdisplay_1.1       backports_1.4.1     BiocManager_1.30.17
    ## [31] IRkernel_1.3        graph_1.74.0        jsonlite_1.8.0     
    ## [34] fs_1.5.2            rjson_0.2.21        digest_0.6.29      
    ## [37] stringi_1.7.8       RJSONIO_1.3-1.6     cli_3.3.0          
    ## [40] tools_4.2.0         bitops_1.0-7        magrittr_2.0.3     
    ## [43] base64url_1.4       RCurl_1.98-1.6      tibble_3.1.7       
    ## [46] crayon_1.5.1        tidyr_1.2.0         pkgconfig_2.0.3    
    ## [49] ellipsis_0.3.2      data.table_1.14.2   assertthat_0.2.1   
    ## [52] rmarkdown_2.14      httr_1.4.3          R6_2.5.1           
    ## [55] compiler_4.2.0

## Creating jupyter files

``` r
#Jupyter Notebook file
if(!"devtools" %in% installed.packages()) BiocManager::install("devtools")
devtools::install_github("mkearney/rmd2jupyter", force=TRUE)
```

    ## 
    ## * checking for file ‘/tmp/RtmpqYwyeu/remotes1c3726fc623c/mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION’ ... OK
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
rmd2jupyter("visualization.Rmd")
```
