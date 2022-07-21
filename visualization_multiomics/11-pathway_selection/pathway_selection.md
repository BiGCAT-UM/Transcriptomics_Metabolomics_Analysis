## Introduction

In this script, we select relevant pathways for both transcriptomics and
metabolomics data.

## Setup

``` r
# check if libraries are already installed > otherwise install it
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"dplyr" %in% installed.packages()) BiocManager::install("dplyr")
#Regular R packages:
if(!"ggplot2" %in% installed.packages()){install.packages("ggplot2")}
if(!"VennDiagram" %in% installed.packages()){install.packages("VennDiagram")}
if(!"RColorBrewer" %in% installed.packages()){install.packages("RColorBrewer")}
#load libraries
library(rstudioapi)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
# set your working environment to the location where your current source file is saved into.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd('../..')
work_DIR <- getwd()
```

##Obtain PW data for transcriptomics and metabolomics

``` r
#Set location to download data for transcriptomics pathway analysis:
filelocation_t <- paste0(work_DIR, "/transcriptomics_analysis/4-pathway_analysis/output/")
#Obtain data from step 4 (transcript PWs)
tPWs_CD_ileum <- read.delim(paste0(filelocation_t, 'enrichResults_ORA_CD_ileum.tsv'), sep = "\t", header = TRUE)
tPWs_CD_rectum <- read.delim(paste0(filelocation_t, 'enrichResults_ORA_CD_rectum.tsv'), sep = "\t",header = TRUE)
tPWs_UC_ileum <- read.delim(paste0(filelocation_t, 'enrichResults_ORA_UC_ileum.tsv'), sep = "\t", header = TRUE)
tPWs_UC_rectum <- read.delim(paste0(filelocation_t, 'enrichResults_ORA_UC_rectum.tsv'), sep = "\t",header = TRUE)

#Set location to download data for metabolomics pathway analysis:
filelocation_m <- paste0(work_DIR, "/metabolomics_analysis/9-metabolite_pathway_analysis/output/")
#Obtain data from step 9 (metabolite PWs)
mPWs_CD <- read.delim(paste0(filelocation_m, 'mbxPWdata_CD.csv'), sep = ",", na.strings=c("", "NA"))
mPWs_UC <- read.delim(paste0(filelocation_m, 'mbxPWdata_UC.csv'), sep = ",", na.strings=c("", "NA"))
#filter out unused columns
#mSet_CD <- mSet_CD [,c(1:4)]
#mSet_UC <- mSet_UC [,c(1:4)]

# Set Working Directory back to current folder
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
work_DIR <- getwd()
```

##Compare PWs from transcriptomics and metabolomics to find overlap

``` r
##Select a disorder:
disorder <- "CD" ##Options CD or UC

#First compare for significant PWs:
tPWs_CD_ileum_sign <- tPWs_CD_ileum[(tPWs_CD_ileum$p.adjust<0.05)&(tPWs_CD_ileum$qvalue<0.02),]
tPWs_CD_rectum_sign <- tPWs_CD_rectum[(tPWs_CD_rectum$p.adjust<0.05)&(tPWs_CD_rectum$qvalue<0.02),]
tPWs_UC_ileum_sign <- tPWs_UC_ileum[(tPWs_UC_ileum$p.adjust<0.05)&(tPWs_UC_ileum$qvalue<0.02),]
tPWs_UC_rectum_sign <- tPWs_UC_rectum[(tPWs_UC_rectum$p.adjust<0.05)&(tPWs_UC_rectum$qvalue<0.02),]

#Cutoff values for significant metabolomics PWs:
##p-value smaller than 0.05
mPWs_CD_sign <- mPWs_CD[(mPWs_CD$probabilities<0.05),]
mPWs_UC_sign <- mPWs_UC[(mPWs_UC$probabilities<0.05),]

#Cutoff values for other interesting metabolomics PWs:
## 3 or more metabolites in the PW, and 5 or more proteins.
mPWs_CD_interest <- mPWs_CD[(mPWs_CD$HMDBsInPWs>3)&(mPWs_CD$ProteinsInPWs>5),]
mPWs_UC_interest <- mPWs_UC[(mPWs_UC$HMDBsInPWs>3)&(mPWs_UC$ProteinsInPWs>5),]

#Cutoff values for significant && interesting metabolomics PWs:
##p-value smaller than 0.05, 3 or more metabolites in the PW, and 5 or more proteins.
mPWs_CD_sign_interest <- mPWs_CD[(mPWs_CD$probabilities<0.05)&(mPWs_CD$HMDBsInPWs>3)&(mPWs_CD$ProteinsInPWs>5),]
mPWs_UC_sign_interest <- mPWs_UC[(mPWs_UC$probabilities<0.05)&(mPWs_UC$HMDBsInPWs>3)&(mPWs_UC$ProteinsInPWs>5),]

set3 <- paste(mPWs_CD_sign_interest[,1] , sep="")
set4 <- paste(mPWs_UC_sign_interest[,1] , sep="")
##Compare both disorders with one another on metabolomics level:
mset_WP_IDs_overlap_sign_interest <- Reduce(intersect, list(set3, set4))

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
##Ignore log messages:
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
```

    ## NULL

``` r
##Check significant PWs only:
if(disorder == "CD"){# Generate 3 sets for comparison:
set1 <- paste(tPWs_CD_ileum_sign[,1] , sep="")
set2 <- paste(tPWs_CD_rectum_sign[,1], sep="")
set3 <- paste(mPWs_CD_sign[,1] , sep="")
set4 <- paste(mPWs_UC_sign[,1] , sep="")
namesDiagram <- c("CD ileum sign" , "CD rectum sign" , "CD metabolites sign")
filenameDiagram <- 'output/CD_comparison_sign.png'}else if(disorder == "UC"){# Generate 3 sets for comparison:
set1 <- paste(tPWs_UC_ileum_sign[,1] , sep="")
set2 <- paste(tPWs_UC_rectum_sign[,1], sep="")
set3 <- paste(mPWs_UC_sign[,1] , sep="")
set4 <- paste(mPWs_CD_sign[,1] , sep="")
namesDiagram <- c("UC ileum sign" , "UC rectum sign" , "UC metabolites sign")
filenameDiagram <- 'output/UC_comparison_sign.png'}else{print("Disorder not recognised.")}

# Chart
venn.diagram(
        x = list(set1, set2, set3),
        category.names = namesDiagram,
        filename = filenameDiagram,
        output=TRUE,
        
        # Output features
        imagetype="png" ,
        height = 480 , 
        width = 480 , 
        resolution = 300,
        compression = "lzw",
        
        # Circles
        lwd = 2,
        lty = 'blank',
        fill = myCol,
        
        # Numbers
        cex = .6,
        fontface = "bold",
        fontfamily = "sans",
        
        # Set names
        cat.cex = 0.6,
        cat.fontface = "bold",
        cat.default.pos = "outer",
        cat.pos = c(-27, 27, 135),
        cat.dist = c(0.055, 0.055, 0.085),
        cat.fontfamily = "sans",
        rotation = 1
)
```

    ## [1] 1

``` r
##Print overlapping PW IDs:
WP_IDs_overlap_sign <- Reduce(intersect, list(set1,set2,set3))
##Compare both disorders with one another on metabolomics level:
mset_WP_IDs_overlap_sign <- Reduce(intersect, list(set3, set4))

#Check interesting PWs only:
if(disorder == "CD"){# Generate 3 sets for comparison:
set1 <- paste(tPWs_CD_ileum_sign[,1] , sep="")
set2 <- paste(tPWs_CD_rectum_sign[,1], sep="")
set3 <- paste(mPWs_CD_interest[,1] , sep="")
set4 <- paste(mPWs_UC_interest[,1] , sep="")
namesDiagram <- c("CD ileum sign" , "CD rectum sign" , "CD metabolites interest")
filenameDiagram <- 'output/CD_comparison_interest.png'}else if(disorder == "UC"){# Generate 3 sets for comparison:
set1 <- paste(tPWs_UC_ileum_sign[,1] , sep="")
set2 <- paste(tPWs_UC_rectum_sign[,1], sep="")
set3 <- paste(mPWs_UC_interest[,1] , sep="")
set4 <- paste(mPWs_CD_interest[,1] , sep="")
namesDiagram <- c("UC ileum sign" , "UC rectum sign" , "UC metabolites interest")
filenameDiagram <- 'output/UC_comparison_interest.png'}else{print("Disorder not recognised.")}

# Chart
venn.diagram(
        x = list(set1, set2, set3),
        category.names = namesDiagram,
        filename = filenameDiagram,
        output=TRUE,
        
        # Output features
        imagetype="png" ,
        height = 480 , 
        width = 480 , 
        resolution = 300,
        compression = "lzw",
        
        # Circles
        lwd = 2,
        lty = 'blank',
        fill = myCol,
        
        # Numbers
        cex = .6,
        fontface = "bold",
        fontfamily = "sans",
        
        # Set names
        cat.cex = 0.6,
        cat.fontface = "bold",
        cat.default.pos = "outer",
        cat.pos = c(-27, 27, 135),
        cat.dist = c(0.055, 0.055, 0.085),
        cat.fontfamily = "sans",
        rotation = 1
)
```

    ## [1] 1

``` r
##Print overlapping PW IDs:
WP_IDs_overlap_interest <- Reduce(intersect, list(set1,set2,set3))
##Compare both disorders with one another on metabolomics level:
mset_WP_IDs_overlap_interest <- Reduce(intersect, list(set3, set4))

##Also compare all pathways:
if(disorder == "CD"){# Generate 3 sets for comparison:
set1 <- paste(tPWs_CD_ileum[,1] , sep="")
set2 <- paste(tPWs_CD_rectum[,1], sep="")
set3 <- paste(mPWs_CD[,1] , sep="")
set4 <- paste(mPWs_UC[,1] , sep="")
namesDiagram <- c("CD ileum" , "CD rectum " , "CD metabolites")
filenameDiagram <- 'output/CD_comparison_All.png'}else if(disorder == "UC"){# Generate 3 sets for comparison:
set1 <- paste(tPWs_UC_ileum[,1] , sep="")
set2 <- paste(tPWs_UC_rectum[,1], sep="")
set3 <- paste(mPWs_UC[,1] , sep="")
set4 <- paste(mPWs_CD[,1] , sep="")
namesDiagram <- c("UC ileum" , "UC rectum " , "UC metabolites")
filenameDiagram <- 'output/UC_comparison_All.png'}else{print("Disorder not recognised.")}

# Chart
venn.diagram(
        x = list(set1, set2, set3),
        category.names = namesDiagram,
        filename = filenameDiagram,
        output=TRUE,
        
        # Output features
        imagetype="png" ,
        height = 480 , 
        width = 480 , 
        resolution = 300,
        compression = "lzw",
        
        # Circles
        lwd = 2,
        lty = 'blank',
        fill = myCol,
        
        # Numbers
        cex = .6,
        fontface = "bold",
        fontfamily = "sans",
        
        # Set names
        cat.cex = 0.6,
        cat.fontface = "bold",
        cat.default.pos = "outer",
        cat.pos = c(-27, 27, 135),
        cat.dist = c(0.055, 0.055, 0.085),
        cat.fontfamily = "sans",
        rotation = 1
)
```

    ## [1] 1

``` r
##Print overlapping PW IDs:
WP_IDs_overlap_all <- Reduce(intersect, list(set1,set2,set3))
##Compare both disorders with one another on metabolomics level:
mset_WP_IDs_overlap_all <- Reduce(intersect, list(set3, set4))

##Paste comparison in as string for printing results to notebook.
string_WP_IDs_overlap_sign_interest_mset <- paste0(mset_WP_IDs_overlap_sign_interest, collapse = ', ')
string_WP_IDs_overlap_sign <- paste0(WP_IDs_overlap_sign, collapse = ', ')
string_WP_IDs_overlap_sign_mset <- paste0(mset_WP_IDs_overlap_sign, collapse = ', ')
string_WP_IDs_overlap_interest <- paste0(WP_IDs_overlap_interest, collapse = ', ')
string_WP_IDs_overlap_interest_mset <- paste0(mset_WP_IDs_overlap_interest, collapse = ', ')
string_WP_IDs_overlap_all <- paste0(WP_IDs_overlap_all, collapse = ', ')
string_WP_IDs_overlap_all_mset <- paste0(mset_WP_IDs_overlap_all, collapse = ', ')

paste0("All significant and interesting metabolic pathways are: ", string_WP_IDs_overlap_sign_interest_mset)
```

    ## [1] "All significant and interesting metabolic pathways are: WP15, WP4726"

``` r
paste0("All significant pathways that overlap for disorder ", disorder," are: ", string_WP_IDs_overlap_sign)
```

    ## [1] "All significant pathways that overlap for disorder CD are: WP1533, WP15"

``` r
paste0("All significant metabolic pathways are: ", string_WP_IDs_overlap_sign_mset)
```

    ## [1] "All significant metabolic pathways are: WP3925, WP15, WP4726, WP3580, WP4949, WP4313, WP4022, WP4504, WP4259, WP465, WP2533, WP4721, WP4719, WP4584, WP4722, WP167"

``` r
paste0("All interesting pathways that overlap for disorder ", disorder," are: ", string_WP_IDs_overlap_interest)
```

    ## [1] "All interesting pathways that overlap for disorder CD are: WP15"

``` r
paste0("All interesting metabolic pathways that overlap are: ", string_WP_IDs_overlap_interest_mset)
```

    ## [1] "All interesting metabolic pathways that overlap are: WP15, WP4723, WP2525, WP4726"

``` r
paste0("All pathways that overlap for disorder ", disorder," are: ", string_WP_IDs_overlap_all)
```

    ## [1] "All pathways that overlap for disorder CD are: WP2431, WP2848, WP5087, WP1533, WP15, WP5094, WP2197, WP98, WP4210, WP727, WP465, WP4341, WP706, WP4723, WP5122, WP4719, WP167, WP4483, WP5046, WP4947, WP3645, WP2435, WP5176, WP4595, WP4584, WP4545, WP3644, WP241, WP497, WP3925, WP3580, WP4224, WP5050, WP4583, WP2267, WP4721, WP2289, WP430, WP2526, WP4313, WP3940, WP4225, WP4290, WP4288, WP4726, WP4686, WP4906, WP2525, WP3869, WP3996, WP4216, WP1601, WP231, WP1982, WP4022"

``` r
paste0("All metabolic pathways that overlap are: ", string_WP_IDs_overlap_all_mset)
```

    ## [1] "All metabolic pathways that overlap are: WP3925, WP15, WP3940, WP4723, WP2525, WP4726, WP550, WP3580, WP4949, WP4313, WP497, WP241, WP706, WP4022, WP4504, WP4259, WP4288, WP4292, WP5122, WP106, WP4583, WP4156, WP4586, WP5190, WP1495, WP528, WP465, WP2533, WP4721, WP4719, WP4584, WP4722, WP167, WP3869, WP3996, WP1601, WP3933, WP5171, WP4220, WP5094, WP2197, WP2333, WP2848, WP98, WP4720, WP1600, WP2526, WP678, WP2431, WP2435, WP2513, WP4483, WP231, WP4216, WP4341"

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
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] RColorBrewer_1.1-3  VennDiagram_1.7.3   futile.logger_1.4.3
    ## [4] ggplot2_3.3.6       dplyr_1.0.9         rstudioapi_0.13    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] formatR_1.12         pillar_1.8.0         compiler_4.2.0      
    ##  [4] BiocManager_1.30.17  futile.options_1.0.1 tools_4.2.0         
    ##  [7] digest_0.6.29        evaluate_0.15        lifecycle_1.0.1     
    ## [10] tibble_3.1.7         gtable_0.3.0         pkgconfig_2.0.3     
    ## [13] rlang_1.0.4          cli_3.3.0            DBI_1.1.3           
    ## [16] yaml_2.3.5           xfun_0.31            fastmap_1.1.0       
    ## [19] withr_2.5.0          stringr_1.4.0        knitr_1.39          
    ## [22] generics_0.1.3       vctrs_0.4.1          tidyselect_1.1.2    
    ## [25] glue_1.6.2           R6_2.5.1             fansi_1.0.3         
    ## [28] rmarkdown_2.14       lambda.r_1.2.4       purrr_0.3.4         
    ## [31] magrittr_2.0.3       scales_1.2.0         ellipsis_0.3.2      
    ## [34] htmltools_0.5.3      assertthat_0.2.1     colorspace_2.0-3    
    ## [37] utf8_1.2.2           stringi_1.7.8        munsell_0.5.0

## Creating jupyter files

``` r
#Jupyter Notebook file
if(!"devtools" %in% installed.packages()) BiocManager::install("devtools")
devtools::install_github("mkearney/rmd2jupyter", force=TRUE)
```

    ## 
    ## * checking for file ‘/tmp/RtmpmoG4sl/remotes287256092117/mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION’ ... OK
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
rmd2jupyter("pathway_selection.Rmd")
```
