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
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
```

## Importing dataset

The data will be read for the disease on two biopsy locations

``` r
## Select a disorder to analyse (options; CD or UC)
disorder <- "CD"
##Obtain data from step 3:
setwd('..')
work_DIR <- getwd()

if (disorder == "CD") {
  dataset_CD <- read.delim("3-identifier_mapping/output/IDMapping_CD.tsv")
  #filter out  unused columns, we select Entrez.ID, log2FC and pvalue, remove NA values: #background genes to be used in enrichment analysis
  dataset <- na.omit(subset( dataset_CD, select = c(3,2,4:7)))
  print("Selected disorder is Crohn's disease")
}else if(disorder == "UC"){ 
  dataset_UC <- read.delim("3-identifier_mapping/output/IDMapping_UC.tsv")
  #filter out  unused columns, we select Entrez.ID, log2FC and pvalue
  dataset <- na.omit(subset( dataset_UC, select = c(3,2,4:7)))
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
up.genes.ileum   <- dataset[dataset$log2FC_ileum >= 0.58 & dataset$pvalue_ileum < 0.05, 2] 
down.genes.ileum <- dataset[dataset$log2FC_ileum <= -0.58 & dataset$pvalue_ileum < 0.05, 2] 
deg.ileum <- unique(dataset[!is.na(dataset$ENTREZ.ID) & !is.na(dataset$pvalue_ileum) & dataset$pvalue_ileum < 0.05 & abs(dataset$log2FC_ileum) > 0.58,c(1:4)])
write.table(deg.ileum,file = paste0("output/DEGs_",disorder,"_ileum.tsv"),sep="\t", quote=FALSE, row.names = FALSE)

#for rectum location
up.genes.rectum   <- dataset[dataset$log2FC_rectum >= 0.58 & dataset$pvalue_rectum < 0.05, 2] 
down.genes.rectum <- dataset[dataset$log2FC_rectum <= -0.58 & dataset$pvalue_rectum < 0.05, 2] 
deg.rectum <- unique(dataset[!is.na(dataset$ENTREZ.ID) & !is.na(dataset$pvalue_rectum) & dataset$pvalue_rectum < 0.05 & abs(dataset$log2FC_rectum) > 0.58,c(1,2,5,6)])
write.table(deg.rectum, file=paste0("output/DEGs_",disorder,"_rectum.tsv"),sep="\t", quote=FALSE, row.names = FALSE)
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

#Count all significant pathways that have p.adjust value lower than 0.05 and qvalue<0.02
ileum.sign <- ewp.ileum.res[(ewp.ileum.res$p.adjust<0.05)&(ewp.ileum.res$qvalue<0.02),]

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

    ## [1] "The number of genes measured in all pathways is: 6026"

``` r
# number of DEGs in all pathways
paste0("The number of DEGs measured in all pathways is: ", length(deg.ileum$ENTREZ.ID[deg.ileum$ENTREZ.ID %in% unique(wp2gene$gene)]))
```

    ## [1] "The number of DEGs measured in all pathways is: 678"

``` r
#number of enriched pathways
paste0("The number of enriched pathways is: ", num.pathways.ileum <- dim(ewp.ileum.res)[1])
```

    ## [1] "The number of enriched pathways is: 471"

``` r
#number of significantly enriched pathways
paste0("The number of significantly enriched pathways is: ", num.pathways.ileum.sign <- dim(ileum.sign)[1])
```

    ## [1] "The number of significantly enriched pathways is: 43"

``` r
#exporting results to the file
write.table(ewp.ileum.res, file=paste0("output/enrichResults_ORA_",disorder,"_ileum.tsv"),
            sep = "\t" ,quote = FALSE, row.names = FALSE)

##################RECTUM location#######################
ewp.rectum <- clusterProfiler::enricher(
  deg.rectum$ENTREZ.ID, # a vector of gene IDs; options
  universe = as.character(dataset$ENTREZ.ID), #background genes to be used in enrichment analysis
  pAdjustMethod = "fdr",#you can change it as BH, bonferronni etc.
  pvalueCutoff = 1, #padjust cutoff
  qvalueCutoff = 1, #q value cutoff 
  TERM2GENE = wpid2gene,
  TERM2NAME = wpid2name)
ewp.rectum.res <- as.data.frame(ewp.rectum) 

#Count all significant pathways that have p.adjust value lower than 0.05 and qvalue<0.02
rectum.sign <- ewp.rectum.res[(ewp.rectum.res$p.adjust<0.05)&(ewp.rectum.res$qvalue<0.02),]

##Print location:
paste0("Pathways enrichment results for disorder: ", disorder , ", location: RECTUM")
```

    ## [1] "Pathways enrichment results for disorder: CD, location: RECTUM"

``` r
# number of genes measured in all pathways
paste0("The number of genes measured in all pathways is: ", length(ewp.rectum@universe))
```

    ## [1] "The number of genes measured in all pathways is: 6026"

``` r
# number of DEGs in all pathways
paste0("The number of DEGs measured in all pathways is: ", length(deg.rectum$ENTREZ.ID[deg.rectum$ENTREZ.ID %in% unique(wp2gene$gene)]))
```

    ## [1] "The number of DEGs measured in all pathways is: 642"

``` r
#number of enriched pathways
paste0("The number of enriched pathways is: ", num.pathways.rectum <- dim(ewp.rectum.res)[1])
```

    ## [1] "The number of enriched pathways is: 483"

``` r
#number of significantly enriched pathways
paste0("The number of significantly enriched pathways is: ", num.pathways.rectum.sign <- dim(rectum.sign)[1])
```

    ## [1] "The number of significantly enriched pathways is: 62"

``` r
#exporting results to the file
write.table(ewp.rectum.res, file=paste0("output/enrichResults_ORA_",disorder,"_rectum.tsv"),
            sep = "\t" ,quote = FALSE, row.names = FALSE)

#if the new data file exist, remove it (so it does not conflict with running the code against the local file) 
if(pathway_data == "new") file.remove(wp.hs.gmt)
```

##Print session info and remove large datasets:

``` r
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
    ##  [1] dplyr_1.1.0           clusterProfiler_4.6.1 rWikiPathways_1.18.0 
    ##  [4] org.Hs.eg.db_3.16.0   AnnotationDbi_1.60.0  IRanges_2.32.0       
    ##  [7] S4Vectors_0.36.2      Biobase_2.58.0        BiocGenerics_0.44.0  
    ## [10] rstudioapi_0.14      
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] nlme_3.1-162           bitops_1.0-7           ggtree_3.6.2          
    ##   [4] enrichplot_1.18.3      bit64_4.0.5            RColorBrewer_1.1-3    
    ##   [7] HDO.db_0.99.1          httr_1.4.5             GenomeInfoDb_1.34.9   
    ##  [10] tools_4.2.2            utf8_1.2.3             R6_2.5.1              
    ##  [13] lazyeval_0.2.2         DBI_1.1.3              colorspace_2.1-0      
    ##  [16] withr_2.5.0            gridExtra_2.3          tidyselect_1.2.0      
    ##  [19] bit_4.0.5              compiler_4.2.2         cli_3.6.0             
    ##  [22] scatterpie_0.1.8       shadowtext_0.1.2       scales_1.2.1          
    ##  [25] yulab.utils_0.0.6      stringr_1.5.0          digest_0.6.31         
    ##  [28] gson_0.0.9             rmarkdown_2.20         DOSE_3.24.2           
    ##  [31] XVector_0.38.0         pkgconfig_2.0.3        htmltools_0.5.4       
    ##  [34] fastmap_1.1.1          rlang_1.0.6            RSQLite_2.3.0         
    ##  [37] gridGraphics_0.5-1     generics_0.1.3         farver_2.1.1          
    ##  [40] jsonlite_1.8.4         BiocParallel_1.32.5    GOSemSim_2.24.0       
    ##  [43] RCurl_1.98-1.10        magrittr_2.0.3         ggplotify_0.1.0       
    ##  [46] GO.db_3.15.0           GenomeInfoDbData_1.2.9 patchwork_1.1.2       
    ##  [49] Matrix_1.5-3           Rcpp_1.0.10            munsell_0.5.0         
    ##  [52] fansi_1.0.4            ape_5.7                viridis_0.6.2         
    ##  [55] lifecycle_1.0.3        stringi_1.7.12         yaml_2.3.7            
    ##  [58] ggraph_2.1.0           MASS_7.3-58.2          zlibbioc_1.44.0       
    ##  [61] plyr_1.8.8             qvalue_2.30.0          grid_4.2.2            
    ##  [64] blob_1.2.3             parallel_4.2.2         ggrepel_0.9.3         
    ##  [67] crayon_1.5.2           lattice_0.20-45        graphlayouts_0.8.4    
    ##  [70] Biostrings_2.66.0      cowplot_1.1.1          splines_4.2.2         
    ##  [73] KEGGREST_1.38.0        knitr_1.42             pillar_1.8.1          
    ##  [76] fgsea_1.24.0           igraph_1.4.1           rjson_0.2.21          
    ##  [79] reshape2_1.4.4         codetools_0.2-19       fastmatch_1.1-3       
    ##  [82] XML_3.99-0.13          glue_1.6.2             evaluate_0.20         
    ##  [85] ggfun_0.0.9            downloader_0.4         data.table_1.14.8     
    ##  [88] BiocManager_1.30.20    treeio_1.22.0          png_0.1-8             
    ##  [91] vctrs_0.5.2            tweenr_2.0.2           gtable_0.3.1          
    ##  [94] purrr_1.0.1            polyclip_1.10-4        tidyr_1.3.0           
    ##  [97] cachem_1.0.7           ggplot2_3.4.1          xfun_0.37             
    ## [100] ggforce_0.4.1          tidygraph_1.2.3        tidytree_0.4.2        
    ## [103] viridisLite_0.4.1      tibble_3.1.8           aplot_0.1.9           
    ## [106] memoise_2.0.1

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
    ## ── R CMD build ─────────────────────────────────────────────────────────────────
    ##          checking for file 'C:\Users\duygu\AppData\Local\Temp\Rtmp2TkvOh\remotes3ffc42055f2b\mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION' ...  ✔  checking for file 'C:\Users\duygu\AppData\Local\Temp\Rtmp2TkvOh\remotes3ffc42055f2b\mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION'
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
rmd2jupyter("pathway_analysis.Rmd")
```
