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
if(!"dplyr" %in% installed.packages()){install.packages("dplyr")}

#loading installed libraries
library(rstudioapi) # interface for interacting with RStudio IDE with R code.
library(org.Hs.eg.db) #This is the organism annotation package ("org") for Homo sapiens ("Hs"), organized as an AnnotationDbi   package ("db"), using Entrez Gene IDs ("eg") as primary key.
library(AnnotationDbi) # for connecting and querying annotation databases
library(dplyr)

# set your working environment to the location where your current source file is saved into.
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
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
write.table(dataset, file=paste0("output/IDMapping_",disorder, ".tsv"),
            sep = "\t" ,quote = FALSE, row.names = FALSE)

##Print session info:
sessionInfo()
```

    ## R version 4.2.2 (2022-10-31 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19044)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_Netherlands.utf8  LC_CTYPE=English_Netherlands.utf8   
    ## [3] LC_MONETARY=English_Netherlands.utf8 LC_NUMERIC=C                        
    ## [5] LC_TIME=English_Netherlands.utf8    
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] dplyr_1.1.0          org.Hs.eg.db_3.16.0  AnnotationDbi_1.60.0
    ## [4] IRanges_2.32.0       S4Vectors_0.36.2     Biobase_2.58.0      
    ## [7] BiocGenerics_0.44.0  rstudioapi_0.14     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.10            pillar_1.8.1           compiler_4.2.2        
    ##  [4] BiocManager_1.30.20    GenomeInfoDb_1.34.9    XVector_0.38.0        
    ##  [7] bitops_1.0-7           tools_4.2.2            zlibbioc_1.44.0       
    ## [10] digest_0.6.31          bit_4.0.5              tibble_3.1.8          
    ## [13] lifecycle_1.0.3        RSQLite_2.3.0          evaluate_0.20         
    ## [16] memoise_2.0.1          pkgconfig_2.0.3        png_0.1-8             
    ## [19] rlang_1.0.6            DBI_1.1.3              cli_3.6.0             
    ## [22] yaml_2.3.7             xfun_0.37              fastmap_1.1.1         
    ## [25] GenomeInfoDbData_1.2.9 withr_2.5.0            httr_1.4.5            
    ## [28] knitr_1.42             generics_0.1.3         Biostrings_2.66.0     
    ## [31] vctrs_0.5.2            tidyselect_1.2.0       bit64_4.0.5           
    ## [34] glue_1.6.2             R6_2.5.1               fansi_1.0.4           
    ## [37] rmarkdown_2.20         magrittr_2.0.3         blob_1.2.3            
    ## [40] htmltools_0.5.4        KEGGREST_1.38.0        utf8_1.2.3            
    ## [43] RCurl_1.98-1.10        cachem_1.0.7           crayon_1.5.2

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
    ## ── R CMD build ─────────────────────────────────────────────────────────────────
    ##          checking for file 'C:\Users\duygu\AppData\Local\Temp\RtmpgV5PLf\remotes4fb456897fd4\mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION' ...  ✔  checking for file 'C:\Users\duygu\AppData\Local\Temp\RtmpgV5PLf\remotes4fb456897fd4\mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION'
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
rmd2jupyter("identifier_mapping.Rmd")
```
