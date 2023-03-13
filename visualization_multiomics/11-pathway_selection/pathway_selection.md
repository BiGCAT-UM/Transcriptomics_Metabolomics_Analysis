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
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#setwd('../..')
setwd('..')
setwd('..')
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
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#setwd("visualization_multiomics")
#work_DIR <- getwd()
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
filenameDiagram <- paste0(work_DIR, "/visualization_multiomics/11-pathway_selection/output/CD_comparison_sign.png")
#filenameDiagram <- 'output/CD_comparison_sign.png'}else if(disorder == "UC"){# Generate 3 sets for comparison:
set1 <- paste(tPWs_UC_ileum_sign[,1] , sep="")
set2 <- paste(tPWs_UC_rectum_sign[,1], sep="")
set3 <- paste(mPWs_UC_sign[,1] , sep="")
set4 <- paste(mPWs_CD_sign[,1] , sep="")
namesDiagram <- c("UC ileum sign" , "UC rectum sign" , "UC metabolites sign")
filenameDiagram <- paste0(work_DIR, "/visualization_multiomics/11-pathway_selection/output/UC_comparison_sign.png")}else{print("Disorder not recognised.")}

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
filenameDiagram <- paste0(work_DIR, "/visualization_multiomics/11-pathway_selection/output/CD_comparison_interest.png")}else if(disorder == "UC"){# Generate 3 sets for comparison:
set1 <- paste(tPWs_UC_ileum_sign[,1] , sep="")
set2 <- paste(tPWs_UC_rectum_sign[,1], sep="")
set3 <- paste(mPWs_UC_interest[,1] , sep="")
set4 <- paste(mPWs_CD_interest[,1] , sep="")
namesDiagram <- c("UC ileum sign" , "UC rectum sign" , "UC metabolites interest")
filenameDiagram <- paste0(work_DIR, "/visualization_multiomics/11-pathway_selection/output/UC_comparison_interest.png")}else{print("Disorder not recognised.")}

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
filenameDiagram <- paste0(work_DIR, "/visualization_multiomics/11-pathway_selection/output/CD_comparison_all.png")}else if(disorder == "UC"){# Generate 3 sets for comparison:
set1 <- paste(tPWs_UC_ileum[,1] , sep="")
set2 <- paste(tPWs_UC_rectum[,1], sep="")
set3 <- paste(mPWs_UC[,1] , sep="")
set4 <- paste(mPWs_CD[,1] , sep="")
namesDiagram <- c("UC ileum" , "UC rectum " , "UC metabolites")
filenameDiagram <- paste0(work_DIR, "/visualization_multiomics/11-pathway_selection/output/UC_comparison_all.png")}else{print("Disorder not recognised.")}

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

    ## [1] "All significant pathways that overlap for disorder CD are: "

``` r
paste0("All significant metabolic pathways are: ", string_WP_IDs_overlap_sign_mset)
```

    ## [1] "All significant metabolic pathways are: WP3604, WP15, WP4726"

``` r
paste0("All interesting pathways that overlap for disorder ", disorder," are: ", string_WP_IDs_overlap_interest)
```

    ## [1] "All interesting pathways that overlap for disorder CD are: WP15"

``` r
paste0("All interesting metabolic pathways that overlap are: ", string_WP_IDs_overlap_interest_mset)
```

    ## [1] "All interesting metabolic pathways that overlap are: WP2525, WP4723, WP15, WP4726"

``` r
paste0("All pathways that overlap for disorder ", disorder," are: ", string_WP_IDs_overlap_all)
```

    ## [1] "All pathways that overlap for disorder CD are: WP15, WP4723, WP5176, WP3925, WP4726, WP2525"

``` r
paste0("All metabolic pathways that overlap are: ", string_WP_IDs_overlap_all_mset)
```

    ## [1] "All metabolic pathways that overlap are: WP3604, WP2525, WP4723, WP15, WP4726, WP550"

##Print session info:

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
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] RColorBrewer_1.1-3  VennDiagram_1.7.3   futile.logger_1.4.3
    ## [4] ggplot2_3.4.1       dplyr_1.1.0         rstudioapi_0.14    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] knitr_1.42           magrittr_2.0.3       munsell_0.5.0       
    ##  [4] tidyselect_1.2.0     colorspace_2.1-0     R6_2.5.1            
    ##  [7] rlang_1.0.6          fastmap_1.1.1        fansi_1.0.4         
    ## [10] tools_4.2.2          gtable_0.3.1         xfun_0.37           
    ## [13] utf8_1.2.3           cli_3.6.0            lambda.r_1.2.4      
    ## [16] withr_2.5.0          htmltools_0.5.4      yaml_2.3.7          
    ## [19] digest_0.6.31        tibble_3.1.8         lifecycle_1.0.3     
    ## [22] formatR_1.14         BiocManager_1.30.20  futile.options_1.0.1
    ## [25] vctrs_0.5.2          glue_1.6.2           evaluate_0.20       
    ## [28] rmarkdown_2.20       compiler_4.2.2       pillar_1.8.1        
    ## [31] generics_0.1.3       scales_1.2.1         pkgconfig_2.0.3

## Creating jupyter files

``` r
#Jupyter Notebook file
if(!"devtools" %in% installed.packages()) BiocManager::install("devtools")
devtools::install_github("mkearney/rmd2jupyter", force=TRUE)
```

    ## 
    ## ── R CMD build ─────────────────────────────────────────────────────────────────
    ##       ✔  checking for file 'C:\Users\duygu\AppData\Local\Temp\Rtmp0WUArv\remotes41ec2b037f70\mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION'
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
rmd2jupyter("pathway_selection.Rmd")
```
