## Introduction

In this workflow, heatmap visualization for enriched pathways will be
performed.

## R environment setup

``` r
# check if libraries are already installed > otherwise install it
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"RColorBrewer" %in% installed.packages()) BiocManager::install("RColorBrewer")  
if(!"dplyr" %in% installed.packages()) BiocManager::install("dplyr")
if(!"pheatmap" %in% installed.packages()) BiocManager::install("pheatmap")

#loading installed libraries
library(RColorBrewer) 
library(dplyr) 
library (pheatmap)
```

## Enriched pathway list will be imported.

``` r
##Obtain data from step 3:
setwd('..')

#we have four datasets in total
#read all pathway lists
CD.ileum <- read.delim("4-pathway_analysis/output/enrichResults_ORA_CD_ileum.tsv",sep = "\t", header = TRUE)
CD.rectum <- read.delim("4-pathway_analysis/output/enrichResults_ORA_CD_rectum.tsv", sep = "\t",header = TRUE)
UC.ileum <- read.delim("4-pathway_analysis/output/enrichResults_ORA_UC_ileum.tsv",sep = "\t", header = TRUE)
UC.rectum <- read.delim("4-pathway_analysis/output/enrichResults_ORA_UC_rectum.tsv", sep = "\t",header = TRUE)

# Set Working Directory back to current folder
setwd("5-create_heatmap")

#we need to get pathways that has p.adjust value lower than 0.05 and qvalue<0.02
#To prevent high false discovery rate (FDR) in multiple testing, q-values are also estimated for FDR control.
CD.ileum.f <- CD.ileum[(CD.ileum$p.adjust<0.05)&(CD.ileum$qvalue<0.02),]
CD.rectum.f <- CD.rectum[(CD.rectum$p.adjust<0.05)&(CD.rectum$qvalue<0.02),]
UC.ileum.f <- UC.ileum[(UC.ileum$p.adjust<0.05)&(UC.ileum$qvalue<0.02),]
UC.rectum.f <- UC.rectum[(UC.rectum$p.adjust<0.05)&(UC.rectum$qvalue<0.02),]
```

## Merge all pathways into a pathway data frame

``` r
#Filter out unused columns 
CD.ileum.f  <- CD.ileum.f [,c(2,6)]
CD.rectum.f <- CD.rectum.f [,c(2,6)]
UC.ileum.f  <- UC.ileum.f [,c(2,6)]
UC.rectum.f <- UC.rectum.f [,c(2,6)]

#first 20 max value of p.adjust pathways for each comparison: Note that if a dataset has less then 20 sign. changed PWs, less rows need to be selected (e.g. adapt c(1:20))
#merge CD pathways
all.pathways.1 <- merge(CD.ileum.f[c(1:20),], CD.rectum.f[c(1:20),],by.x="Description", by.y="Description",sort = TRUE, all.x = TRUE, all.y = TRUE)
#merge UC pathways
all.pathways.2 <- merge(UC.ileum.f, UC.rectum.f[c(1:20),],by.x="Description", by.y="Description",sort = TRUE, all.x = TRUE, all.y = TRUE)
#merge all of them
all.pathways <- merge(all.pathways.1 , all.pathways.2 ,by.x="Description", by.y="Description",sort = TRUE, all.x = TRUE, all.y = TRUE)
colnames(all.pathways) <- c("Description","CD.ileum.p.adjust","CD.rectum.p.adjust","UC.ileum.p.adjust","UC.rectum.p.adjust")
#remove unused variables
rm(all.pathways.1, all.pathways.2)
```

## Modify merged pathway list

``` r
#replace NA values with the values from the whole list
#### for CD ileum
#find pathways which does not occur in the filtered enriched pathway list of cd.ileum (p.adjust<0.05 & qvalue<0.02 )
#because we will not take into account not sig. enriched pathways, we will assign value of 1 for their p.adjust
notExist.CDileum <- setdiff(all.pathways$Description,CD.ileum.f$Description)
all.pathways[all.pathways$Description %in% notExist.CDileum,]$CD.ileum.p.adjust <- 1

#the rest NA values correspond to the sig.enriched pathways but not in the first 20 list.
#so we will replace NA values with the p.adjust values from the whole list
NA.indices <- which(is.na(all.pathways$CD.ileum.p.adjust), arr.ind = TRUE)
allIDs <- all.pathways[NA.indices,]$Description
df <- CD.ileum.f[CD.ileum.f$Description %in% allIDs,]
df <- df[order(df$Description),]
all.pathways[all.pathways$Description %in% df$Description,]$CD.ileum.p.adjust <- df$p.adjust

#### for CD rectum
#find pathways which does not occur in the filtered enriched pathway list of cd rectum (p.adjust<0.05 & qvalue<0.02 )
notExist.CDrectum<- setdiff(all.pathways$Description,CD.rectum.f$Description)
all.pathways[all.pathways$Description %in% notExist.CDrectum,]$CD.rectum.p.adjust <- 1

#replacing NA values with the values from the whole list
NA.indices <- which(is.na(all.pathways$CD.rectum.p.adjust), arr.ind = TRUE)
allIDs <- all.pathways[NA.indices,]$Description
df <- CD.rectum.f[CD.rectum.f$Description %in% allIDs,]
df <- df[order(df$Description),]
all.pathways[all.pathways$Description %in% df$Description,]$CD.rectum.p.adjust <- df$p.adjust

#### for UC ileum
#find pathways which does not occur in the filtered enriched pathway list of UC ileum (p.adjust<0.05 & qvalue<0.02 )
notExist.UCileum <- setdiff(all.pathways$Description,UC.ileum.f$Description)
all.pathways[all.pathways$Description %in% notExist.UCileum,]$UC.ileum.p.adjust <- 1

#### for UC rectum
#find pathways which does not occur in the filtered enriched pathway list of UC rectum (p.adjust<0.05 & qvalue<0.02 )
notExist.UCrectum<- setdiff(all.pathways$Description,UC.rectum.f$Description)
all.pathways[all.pathways$Description %in% notExist.UCrectum,]$UC.rectum.p.adjust <- 1

#replacing NA values with the values from the whole list
NA.indices <- which(is.na(all.pathways$UC.rectum.p.adjust), arr.ind = TRUE)
allIDs <- all.pathways[NA.indices,]$Description
df <- UC.rectum.f[UC.rectum.f$Description %in% allIDs,]
df <- df[order(df$Description),]
all.pathways[all.pathways$Description %in% df$Description,]$UC.rectum.p.adjust <- df$p.adjust
```

## Heatmap visualization

``` r
#take only required columns
row.names(all.pathways) <- all.pathways$Description
all.pathways  <- all.pathways[,2:5]
colnames(all.pathways) <- c("CD.ileum","CD.rectum","UC.ileum","UC.rectum")

#create output folder if not exist
if(!dir.exists("output")) dir.create("output")

## Select a size to visualize the heatmap with (options; large or small)
size_heatmap <- "large"

##Print labels large for paper, small for notebook:
fontsize_row_l = 30 
if (size_heatmap == "large") {
fontsize_col_l = 30 
fontsize_l = 30
width_l =2000 
height_l =2000 
name_heatmap_file <- "output/heatmap_log10_large.png"
}else if(size_heatmap == "small"){ 
fontsize_row_l = 10 
fontsize_col_l = 10 
width_l =1500 
height_l =1500 
fontsize_l = 10
name_heatmap_file <- "output/heatmap_log10_small.png"
}else{print("Size not Recognised")}

#normally darker value represent higher values light color represent smaller values
#when we use rev function higher ones are represented by light color
colMain <- colorRampPalette(rev(brewer.pal(9, "Blues")))(30)

my_heatmap <- pheatmap(as.matrix(log10(all.pathways)), scale = "none", color = colMain , 
                       legend = TRUE , legend_breaks = c(0, -5, -10, -15, min(log10(all.pathways))), 
                       main = "", 
                       legend_labels = c("adj. p-values \n", " -5", " -10", " -15", ""),
                       cellwidth = 80, treeheight_row = 200, fontsize = fontsize_l, fontsize_row= fontsize_row_l, 
                       fontsize_col = fontsize_col_l, cluster_rows = TRUE, cluster_cols = FALSE)
```

![](heatMap_files/figure-markdown_github/heatmap-1.png)

``` r
#save obtained heatmap
save_pheatmap_png <- function(x, filename, width = width_l, height = height_l) {
  png(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
##p-values are visalized on a log10 scale, to make them more discriminatory.
save_pheatmap_png(my_heatmap, name_heatmap_file)
```

    ## png 
    ##   2

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] pheatmap_1.0.12    dplyr_1.1.0        RColorBrewer_1.1-3
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] rstudioapi_0.14     knitr_1.42          magrittr_2.0.3     
    ##  [4] munsell_0.5.0       tidyselect_1.2.0    colorspace_2.1-0   
    ##  [7] R6_2.5.1            rlang_1.0.6         fastmap_1.1.1      
    ## [10] fansi_1.0.4         highr_0.10          tools_4.2.2        
    ## [13] grid_4.2.2          gtable_0.3.1        xfun_0.37          
    ## [16] utf8_1.2.3          cli_3.6.0           htmltools_0.5.4    
    ## [19] yaml_2.3.7          digest_0.6.31       tibble_3.1.8       
    ## [22] lifecycle_1.0.3     BiocManager_1.30.20 vctrs_0.5.2        
    ## [25] glue_1.6.2          evaluate_0.20       rmarkdown_2.20     
    ## [28] compiler_4.2.2      pillar_1.8.1        scales_1.2.1       
    ## [31] generics_0.1.3      pkgconfig_2.0.3

``` r
##Remove data objects which are not needed for further processing:
rm(list=setdiff(ls(), "all.pathways"))
```

## Last, we create a Jupyter notebook from this script

``` r
#Jupyter Notebook file
if(!"devtools" %in% installed.packages()) BiocManager::install("devtools")
devtools::install_github("mkearney/rmd2jupyter", force=TRUE)
```

    ## 
    ## ── R CMD build ─────────────────────────────────────────────────────────────────
    ##          checking for file 'C:\Users\duygu\AppData\Local\Temp\RtmpcLGaVU\remotes18a862fa58a6\mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION' ...  ✔  checking for file 'C:\Users\duygu\AppData\Local\Temp\RtmpcLGaVU\remotes18a862fa58a6\mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION'
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
rmd2jupyter("heatMap.Rmd")
```
