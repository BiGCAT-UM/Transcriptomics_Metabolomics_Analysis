## Introduction

In this section, identifier (IDs) mapping is performed from the original
data annotation (HGNC symbols) to Entrez Gene and Ensembl IDs, since
tools downstream of this step require different input formats for the
IDs.

## R environment setup

``` r
# check if libraries are already installed > otherwise install it
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"org.Hs.eg.db" %in% installed.packages()) BiocManager::install("org.Hs.eg.db")  
if(!"AnnotationDbi" %in% installed.packages()) BiocManager::install("AnnotationDbi")
#if(!"rWikiPathways" %in% installed.packages()) BiocManager::install("rWikiPathways")
#if(!"clusterProfiler" %in% installed.packages()) BiocManager::install("clusterProfiler") 
if(!"dplyr" %in% installed.packages()){install.packages("dplyr")}

#loading installed libraries
library(rstudioapi) # interface for interacting with RStudio IDE with R code.
library(org.Hs.eg.db) #This is the organism annotation package ("org") for Homo sapiens ("Hs"), organized as an AnnotationDbi   package ("db"), using Entrez Gene IDs ("eg") as primary key.
library(AnnotationDbi) # for connecting and querying annotation databases
#library(rWikiPathways) # for programmatic access to WikiPathways content
#library(clusterProfiler) # for implementing methods to analyze and visualize functional profiles of genomic data
library(dplyr)

# set your working environment to the location where your current source file is saved into.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
```

## Importing dataset

The data will be read for the disease on two biopsy locations

``` r
## Select a disorder to analyse (options; CD or UC)
disorder <- "CD"
##Obtain data from step 2:
setwd('..')
work_DIR <- getwd()
#we have two datasets from different biopsy locations
dataset1 <- read.delim("2-differential_gene_expression_analysis/statsmodel/table_UC_Ileum_vs_nonIBD_Ileum.tab")
dataset2 <- read.delim("2-differential_gene_expression_analysis/statsmodel/table_UC_Rectum_vs_nonIBD_Rectum.tab")
dataset3 <- read.delim("2-differential_gene_expression_analysis/statsmodel/table_CD_Ileum_vs_nonIBD_Ileum.tab")
dataset4 <- read.delim("2-differential_gene_expression_analysis/statsmodel/table_CD_Rectum_vs_nonIBD_Rectum.tab")

# Set Working Directory back to current folder
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
work_DIR <- getwd()

if (disorder == "CD") {
  #filter out  unused columns, we select geneSymbol, log2FC and pvalue
  dataset_ileum<- subset( dataset3, select = c(1,3,7))
  dataset_rectum<- subset( dataset4, select = c(1,3,7))
  print("Selected disorder is Crohn's disease")
}else if(disorder == "UC"){ 
  #filter out  unused columns, we select geneSymbol, log2FC and pvalue
  dataset_ileum<- subset( dataset1, select = c(1,3,7))
  dataset_rectum<- subset( dataset2, select = c(1,3,7))
  print("Selected disorder is Ulcerative Colitis")}else{print("Disorder not Recognised")
}
```

    ## [1] "Selected disorder is Crohn's disease"

``` r
#merge two dataset of two locations into one data 
dataset <- merge(dataset_ileum, dataset_rectum,by.x="X", by.y="X",sort = TRUE, all.x = TRUE, all.y = TRUE)
#change column names
colnames(dataset) <- c("GeneSymbol","log2FC_ileum","pvalue_ileum","log2FC_rectum","pvalue_rectum")
```

## Converting hgnc gene symbols to the corresponding Entrez (NCBI) gene IDs

``` r
#converting gene symbols to entrez ID since these are required for the enrichR function
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
#dataset<- dataset %>% tidyr::drop_na(ENTREZ.ID)
```

## Converting hgnc gene symbols to the corresponding Ensembl IDs

``` r
#converting gene symbols to Ensembl ID since these are required for the Cytoscape multiomics visualization
hs <- org.Hs.eg.db #This object is a simple mapping of Entrez Gene identifier
ensemblID <- AnnotationDbi::select(hs, keys = dataset$GeneSymbol,
            columns = c("ENSEMBL", "SYMBOL"),
            keytype = "SYMBOL")
#filter out double gene symbols
ensemblID <- ensemblID %>% distinct (ensemblID$SYMBOL, .keep_all = TRUE)
# add entrezIDs for each gene symbol in the dataset
dataset <- cbind(ensemblID$ENSEMBL,dataset)
#change column name
colnames(dataset)[1] = "Ensembl.ID"
#filter out genes that has NA value for entrezID
#dataset<- dataset %>% tidyr::drop_na(Ensembl.ID)


##TODO: add NA removal before PW analysis (and network analysis?)
##TODO: print some stats on mapping (issues)
```

##Save data, print session info and remove large datasets:

``` r
##Save data:
#exporting results to the file
write.table(dataset, file=paste0("output/IDMapping_",disorder),
            sep = "\t" ,quote = FALSE, row.names = FALSE)

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
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] dplyr_1.0.9          org.Hs.eg.db_3.15.0  AnnotationDbi_1.58.0
    ## [4] IRanges_2.30.0       S4Vectors_0.34.0     Biobase_2.56.0      
    ## [7] BiocGenerics_0.42.0  rstudioapi_0.13     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] KEGGREST_1.36.2        tidyselect_1.1.2       xfun_0.31             
    ##  [4] purrr_0.3.4            vctrs_0.4.1            generics_0.1.2        
    ##  [7] htmltools_0.5.2        yaml_2.3.5             utf8_1.2.2            
    ## [10] blob_1.2.3             rlang_1.0.2            pillar_1.7.0          
    ## [13] glue_1.6.2             DBI_1.1.2              bit64_4.0.5           
    ## [16] GenomeInfoDbData_1.2.8 lifecycle_1.0.1        stringr_1.4.0         
    ## [19] zlibbioc_1.42.0        Biostrings_2.64.0      memoise_2.0.1         
    ## [22] evaluate_0.15          knitr_1.39             fastmap_1.1.0         
    ## [25] GenomeInfoDb_1.32.2    fansi_1.0.3            Rcpp_1.0.8.3          
    ## [28] BiocManager_1.30.17    cachem_1.0.6           XVector_0.36.0        
    ## [31] bit_4.0.4              png_0.1-7              digest_0.6.29         
    ## [34] stringi_1.7.6          cli_3.3.0              tools_4.2.0           
    ## [37] bitops_1.0-7           magrittr_2.0.3         RCurl_1.98-1.6        
    ## [40] RSQLite_2.2.13         tibble_3.1.7           crayon_1.5.1          
    ## [43] pkgconfig_2.0.3        ellipsis_0.3.2         assertthat_0.2.1      
    ## [46] rmarkdown_2.14         httr_1.4.3             R6_2.5.1              
    ## [49] compiler_4.2.0

``` r
##Remove data objects which are not needed for further processing:
rm(list=setdiff(ls(), c("dataset", "disorder", "work_DIR")))
```

### Last, we create a Jupyter notebook from this script

``` r
#Jupyter Notebook file
if(!"devtools" %in% installed.packages()) BiocManager::install("devtools")
devtools::install_github("mkearney/rmd2jupyter", force=TRUE)
```

    ## 
    ## * checking for file ‘/tmp/RtmpRCEeFx/remotes47da424a42f8/mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION’ ... OK
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
rmd2jupyter("identifier_mapping.Rmd")
```
