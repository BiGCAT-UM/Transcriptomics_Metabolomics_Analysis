## Introduction

In this workflow, we will perform pathway enrichment analysis to
differential expressed genes for both diseases on two biopsy locations
ileum and rectum. We will be performing pathway enrichment analysis on a
differential gene expression dataset. The dataset compares the
expression of transcripts in inflammatory bowel disease (CD and UC)
biopsy locations (ileum and rectum) versus healthy people. Differential
expression analysis has already been performed (see step 2 of this
workflow), generating log2foldchange and (adjusted) p-values data for
each gene.

## R environment setup

``` r
# check if libraries are already installed > otherwise install it
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"org.Hs.eg.db" %in% installed.packages()) BiocManager::install("org.Hs.eg.db")  
if(!"AnnotationDbi" %in% installed.packages()) BiocManager::install("AnnotationDbi")
if(!"rWikiPathways" %in% installed.packages()) BiocManager::install("rWikiPathways")
if(!"clusterProfiler" %in% installed.packages()) BiocManager::install("clusterProfiler") 
if(!"dplyr" %in% installed.packages()){install.packages("dplyr")}

#loading installed libraries
library(rstudioapi) # interface for interacting with RStudio IDE with R code.
library(org.Hs.eg.db) #This is the organism annotation package ("org") for Homo sapiens ("Hs"), organized as an AnnotationDbi   package ("db"), using Entrez Gene IDs ("eg") as primary key.
library(AnnotationDbi) # for connecting and querying annotation databases
library(rWikiPathways) # for programmatic access to WikiPathways content
library(clusterProfiler) # for implementing methods to analyze and visualize functional profiles of genomic data
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

if (disorder == "CD") {
  dataset_CD <- read.delim("3-identifier_mapping/output/IDMapping_CD.tsv")
  #filter out  unused columns, we select Entrez.ID, log2FC and pvalue, remove NA values: #background genes to be used in enrichment analysis
  dataset <- na.omit(subset( dataset_CD, select = c(2,4:7)))
  print("Selected disorder is Crohn's disease")
}else if(disorder == "UC"){ 
  dataset_UC <- read.delim("3-identifier_mapping/output/IDMapping_UC.tsv")
  #filter out  unused columns, we select Entrez.ID, log2FC and pvalue
  dataset <- na.omit(subset( dataset_UC, select = c(2,4:7)))
  print("Selected disorder is Ulcerative Colitis")}else{print("Disorder not Recognised")
  }
```

    ## [1] "Selected disorder is Crohn's disease"

``` r
# Set Working Directory back to current folder
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
work_DIR <- getwd()
```

## Getting differentially expressed genes including up-regulated and down-regulated genes using the criteria

``` r
#if the results folder does not exist create it 
if(!dir.exists("output")) dir.create("output")
#we will use selection criteria as Fold change=1.5,log2FC=0.58 and p.value < 0.05
#for ileum location
up.genes.ileum   <- dataset[dataset$log2FC_ileum >= 0.58 & dataset$pvalue_ileum < 0.05, 1] 
down.genes.ileum <- dataset[dataset$log2FC_ileum <= -0.58 & dataset$pvalue_ileum < 0.05, 1] 
deg.ileum <- unique(dataset[!is.na(dataset$ENTREZ.ID) & !is.na(dataset$pvalue_ileum) & dataset$pvalue_ileum < 0.05 & abs(dataset$log2FC_ileum) > 0.58,c(1:3)])
write.table(deg.ileum,file = paste0("output/DEGs_",disorder,"_ileum"),sep="\t", quote=FALSE, row.names = FALSE)

#for rectum location
up.genes.rectum   <- dataset[dataset$log2FC_rectum >= 0.58 & dataset$pvalue_rectum < 0.05, 1] 
down.genes.rectum <- dataset[dataset$log2FC_rectum <= -0.58 & dataset$pvalue_rectum < 0.05, 1] 
deg.rectum <- unique(dataset[!is.na(dataset$ENTREZ.ID) & !is.na(dataset$pvalue_rectum) & dataset$pvalue_rectum < 0.05 & abs(dataset$log2FC_rectum) > 0.58,c(1,4:5)])
write.table(deg.rectum, file=paste0("output/DEGs_",disorder,"_rectum"),sep="\t", quote=FALSE, row.names = FALSE)
```

## Pathway Enrichment Analysis

In this section, we will perform pathway enrichment analysis. The
clusterProfiler R-package is used to perform over-representation
analysis (ORA). The function can be easily replaced to use other
enrichment methods (GSEA / rSEA / etc).

``` r
##Work with local file (for publication), or new download:
pathway_data <- "local" #Options: local, new
if (pathway_data == "local") {
  wp.hs.gmt <-list.files(work_DIR, pattern="wikipathways", full.names=FALSE)
  paste0("Using local file, from: ", wp.hs.gmt )
}else if(pathway_data == "new"){ 
#below code should be performed first to handle the ssl certificate error while downloading pathways 
options(RCurlOptions = list(cainfo = paste0( tempdir() , "/cacert.pem" ), ssl.verifypeer = FALSE))
#downloading latest pathway gmt files for human 
wp.hs.gmt <- rWikiPathways::downloadPathwayArchive(organism="Homo sapiens", format = "gmt")
  paste0("Using new data, from: ", wp.hs.gmt)}else{print("Pathway data type not recognized")
  }
```

    ## [1] "Using local file, from: wikipathways-20220510-gmt-Homo_sapiens.gmt"

``` r
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
  universe = as.character(dataset$ENTREZ.ID), #background genes to be used in enrichment analysis
  pAdjustMethod = "fdr",#you can change this to hochberg, bonferronni or none etc.
  pvalueCutoff = 1, #adjusted pvalue cutoff on enrichment tests to report, we set it a wider criteria then we will filter out
  #results based on padjust and qvalue in next section which is enrichment result visualization
  qvalueCutoff = 1, #qvalue cutoff on enrichment tests to report as significant, 
                    #multiple hypothesis testing
  TERM2GENE = wpid2gene, #user input annotation of TERM TO GENE mapping
  TERM2NAME = wpid2name) #user input of TERM TO NAME mapping

ewp.ileum.res <- as.data.frame(ewp.ileum) 

#interpretation of the output: 
#     BgRatio   = (number of genes measured in the current pathway) / (number of genes measured in all pathways)
#     geneRatio = (number of DEGs in the current pathway) / (total number of DEGs in all pathways)

##Print location:
paste0("Pathways enrichment results for disorder: ", disorder , ", location: ILEUM")
```

    ## [1] "Pathways enrichment results for disorder: CD, location: ILEUM"

``` r
# number of genes measured in all pathways
paste0("The number of genes measured in all pathways is: ", length(ewp.ileum@universe))
```

    ## [1] "The number of genes measured in all pathways is: 6030"

``` r
# number of DEGs in all pathways
paste0("The number of DEGs measured in all pathways is: ", length(deg.ileum$ENTREZ.ID[deg.ileum$ENTREZ.ID %in% unique(wp2gene$gene)]))
```

    ## [1] "The number of DEGs measured in all pathways is: 679"

``` r
#number of enriched pathways
paste0("The number of enriched pathways is: ", num.pathways.ileum <- dim(ewp.ileum.res)[1])
```

    ## [1] "The number of enriched pathways is: 472"

``` r
#exporting results to the file
write.table(ewp.ileum.res, file=paste0("output/enrichResults_ORA_",disorder,"_ileum"),
            sep = "\t" ,quote = FALSE, row.names = FALSE)

##################RECTUM location#######################
ewp.rectum <- clusterProfiler::enricher(
  deg.rectum$ENTREZ.ID,
  universe = as.character(dataset$ENTREZ.ID), #background genes to be used in enrichment analysis
  pAdjustMethod = "fdr",#you can change it as BH, bonferronni etc.
  pvalueCutoff = 1, #padjust cutoff
  qvalueCutoff = 1, #q value cutoff 
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)
ewp.rectum.res <- as.data.frame(ewp.rectum) 

##Print location:
paste0("Pathways enrichment results for disorder: ", disorder , ", location: RECTUM")
```

    ## [1] "Pathways enrichment results for disorder: CD, location: RECTUM"

``` r
# number of genes measured in all pathways
paste0("The number of genes measured in all pathways is: ", length(ewp.rectum@universe))
```

    ## [1] "The number of genes measured in all pathways is: 6030"

``` r
# number of DEGs in all pathways
paste0("The number of DEGs measured in all pathways is: ", length(deg.rectum$ENTREZ.ID[deg.rectum$ENTREZ.ID %in% unique(wp2gene$gene)]))
```

    ## [1] "The number of DEGs measured in all pathways is: 643"

``` r
#number of enriched pathways
paste0("The number of enriched pathways is: ", num.pathways.rectum <- dim(ewp.rectum.res)[1])
```

    ## [1] "The number of enriched pathways is: 484"

``` r
#exporting results to the file
write.table(ewp.rectum.res, file=paste0("output/enrichResults_ORA_",disorder,"_rectum"),
            sep = "\t" ,quote = FALSE, row.names = FALSE)

#if the new data file exist, remove it (so it does not conflict with running the code against the local file) 
if(pathway_data == "new") file.remove(wp.hs.gmt)
```

##Print session info and remove large datasets:

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
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] dplyr_1.0.9           clusterProfiler_4.4.3 rWikiPathways_1.16.0 
    ##  [4] org.Hs.eg.db_3.15.0   AnnotationDbi_1.58.0  IRanges_2.30.0       
    ##  [7] S4Vectors_0.34.0      Biobase_2.56.0        BiocGenerics_0.42.0  
    ## [10] rstudioapi_0.13      
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] fgsea_1.22.0           colorspace_2.0-3       ggtree_3.4.0          
    ##   [4] rjson_0.2.21           ellipsis_0.3.2         qvalue_2.28.0         
    ##   [7] XVector_0.36.0         aplot_0.1.6            farver_2.1.0          
    ##  [10] graphlayouts_0.8.0     ggrepel_0.9.1          bit64_4.0.5           
    ##  [13] fansi_1.0.3            scatterpie_0.1.7       splines_4.2.0         
    ##  [16] cachem_1.0.6           GOSemSim_2.22.0        knitr_1.39            
    ##  [19] polyclip_1.10-0        jsonlite_1.8.0         GO.db_3.15.0          
    ##  [22] png_0.1-7              ggforce_0.3.3          BiocManager_1.30.17   
    ##  [25] compiler_4.2.0         httr_1.4.3             assertthat_0.2.1      
    ##  [28] Matrix_1.4-1           fastmap_1.1.0          lazyeval_0.2.2        
    ##  [31] cli_3.3.0              tweenr_1.0.2           htmltools_0.5.2       
    ##  [34] tools_4.2.0            igraph_1.3.1           gtable_0.3.0          
    ##  [37] glue_1.6.2             GenomeInfoDbData_1.2.8 reshape2_1.4.4        
    ##  [40] DO.db_2.9              fastmatch_1.1-3        Rcpp_1.0.8.3          
    ##  [43] enrichplot_1.16.1      vctrs_0.4.1            Biostrings_2.64.0     
    ##  [46] ape_5.6-2              nlme_3.1-157           ggraph_2.0.5          
    ##  [49] xfun_0.31              stringr_1.4.0          lifecycle_1.0.1       
    ##  [52] XML_3.99-0.9           DOSE_3.22.0            zlibbioc_1.42.0       
    ##  [55] MASS_7.3-57            scales_1.2.0           tidygraph_1.2.1       
    ##  [58] parallel_4.2.0         RColorBrewer_1.1-3     yaml_2.3.5            
    ##  [61] memoise_2.0.1          gridExtra_2.3          ggplot2_3.3.6         
    ##  [64] downloader_0.4         ggfun_0.0.6            yulab.utils_0.0.4     
    ##  [67] stringi_1.7.6          RSQLite_2.2.13         tidytree_0.3.9        
    ##  [70] BiocParallel_1.30.0    GenomeInfoDb_1.32.2    rlang_1.0.2           
    ##  [73] pkgconfig_2.0.3        bitops_1.0-7           evaluate_0.15         
    ##  [76] lattice_0.20-45        purrr_0.3.4            treeio_1.20.0         
    ##  [79] patchwork_1.1.1        shadowtext_0.1.2       bit_4.0.4             
    ##  [82] tidyselect_1.1.2       plyr_1.8.7             magrittr_2.0.3        
    ##  [85] R6_2.5.1               generics_0.1.2         DBI_1.1.2             
    ##  [88] pillar_1.7.0           KEGGREST_1.36.2        RCurl_1.98-1.6        
    ##  [91] tibble_3.1.7           crayon_1.5.1           utf8_1.2.2            
    ##  [94] rmarkdown_2.14         viridis_0.6.2          grid_4.2.0            
    ##  [97] data.table_1.14.2      blob_1.2.3             digest_0.6.29         
    ## [100] tidyr_1.2.0            gridGraphics_0.5-1     munsell_0.5.0         
    ## [103] viridisLite_0.4.0      ggplotify_0.1.0

``` r
##Remove data objects which are not needed for further processing:
rm(list=setdiff(ls(), c("ewp.rectum.res", "ewp.ileum.res")))
```

### Last, we create a Jupyter notebook from this script

``` r
#Jupyter Notebook file
if(!"devtools" %in% installed.packages()) BiocManager::install("devtools")
devtools::install_github("mkearney/rmd2jupyter", force=TRUE)
```

    ## 
    ## * checking for file ‘/tmp/RtmpV4lWq4/remotes6aef20db5eeb/mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION’ ... OK
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
