## Introduction

In this script, visualization of enriched pathway which include both
altered genes and metabolites is performed

## Setup

``` r
# check if libraries are already installed > otherwise install it
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"RCy3" %in% installed.packages()) BiocManager::install("RCy3")
if(!"rWikiPathways" %in% installed.packages()) BiocManager::install("rWikiPathways")
if(!"RColorBrewer" %in% installed.packages()) BiocManager::install("RColorBrewer")

#load libraries
library(rstudioapi)
library(RCy3)#connect cytoscape via R
library(rWikiPathways)#to get pathways from WikiPathways
library(RColorBrewer)#to manage colors with R

# set your working environment to the location where your current source file is saved into.
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
```

``` r
#make sure to launch cytoscape, if you get CyREST error you need to relaunch cytoscape
cytoscapePing()
#close all opened session before starting
closeSession(FALSE)# close all opened sessions data without saving
#pathway IDs
pathway.id <- "WP4723"# omega-3/omega-6 fatty acid synthesis pathway is one of the enriched pathways in our study
#so we will visualize the pathway in cytoscape

##### STEP-1 Open enriched pathway #########
#import pathways as network in cytoscape
RCy3::commandsRun(paste0('wikipathways import-as-pathway id=',pathway.id )) 
#RCy3::commandsRun(paste0('wikipathways import-as-network id=',pathway.id )) 
```

``` r
#upload gene data table from the file
geneData <- read.delim("data/transcriptomics")
#load data to the pathway in cytoscape
loadTableData(table = "node", data = geneData, data.key.column = "ENSEMBL.ID", table.key.column = "Ensembl")
```

    ## [1] "Success: Data loaded in defaultnode table"

``` r
#read mbx data
mbxData <- read.delim("data/metabolomics")
#load mbx data to the pathway in cytoscape
loadTableData(table = "node", data = mbxData, data.key.column = "CHEBI", table.key.column = "ChEBI")
```

    ## [1] "Success: Data loaded in defaultnode table"

``` r
##new visual style created ##
RCy3::copyVisualStyle("default","pathwayStyle")
RCy3::setVisualStyle("pathwayStyle")
```

    ##                 message 
    ## "Visual Style applied."

``` r
RCy3::lockNodeDimensions(TRUE, style.name="pathwayStyle")

#node shape mapping
RCy3::setNodeShapeMapping('Type',c('Protein','GeneProduct','Metabolite'),c('ELLIPSE','ELLIPSE','RECTANGLE'),style.name="pathwayStyle")
```

    ## NULL

``` r
#node height mapping
RCy3::setNodeHeightMapping('Type',c('Protein','GeneProduct','Metabolite'), c(20,20,25), mapping.type = "d", style.name = "pathwayStyle")
```

    ## NULL

``` r
#change node sizes
RCy3::setNodeSizeMapping('Type', c('Protein','GeneProduct','Metabolite'), c(30,30,40), mapping.type = "d", style.name = "pathwayStyle")
```

    ## NULL

``` r
#change node sizes
RCy3::setNodeWidthMapping('Type', c('Protein','GeneProduct','Metabolite'), c(55,55,150), mapping.type = "d", style.name = "pathwayStyle")
```

    ## NULL

``` r
#set node color based on log2FC
log2FCgenes = getTableColumns('node','log2FC_gene')  
min.logFC = min(log2FCgenes[,1],na.rm=TRUE)
max.logFC = max(log2FCgenes[,1],na.rm=TRUE)
data.values = c(min.logFC,0,max.logFC)
display.brewer.all(length(data.values), colorblindFriendly=TRUE, type="div") # div,qual,seq,all
```

![](visualization_files/figure-markdown_github/visualization-1.png)

``` r
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping('log2FC_gene', data.values, node.colors, style.name="pathwayStyle", default = "#D3D3D3")
```

    ## NULL

``` r
#set node border based on p-value
pvalues = getTableColumns('node','pvalue_gene')  
min.pval = min(pvalues[,1],na.rm=TRUE)
setNodeBorderWidthMapping('pvalue_gene', c(min.pval,0.05), c(5,1),style.name = "pathwayStyle")
```

    ## NULL

``` r
setNodeBorderColorMapping('pvalue_gene', c(min.pval,0.05), c('#00FF00','#00FF00'),style.name = "pathwayStyle")
```

    ## NULL

``` r
#color metabolite nodes byPass
all.nodes = getTableColumns(columns = c("name", "data.type"))  
met.nodes <- na.omit(all.nodes[(all.nodes$data.type == "metabolomics"),])
setNodeColorBypass(met.nodes[,1],"#FD39B8")

#save image 
png.file <- file.path(getwd(), "output/multi_omics_visualization.png")
exportImage(png.file,'PNG', zoom = 500)
```

    ##                                                                                                                                                                           file 
    ## "C:\\Users\\dedePC\\Desktop\\FNS_CLOUD\\GITHUB_codes\\Transcriptomics_Metabolomics_Analysis\\visualization_multiomics\\2-visualization\\output\\multi_omics_visualization.png"
