## Introduction

In this workflow, we will create protein-protein interaction (PPI)
networks for both biopsy locations ileum and rectum. Then these networks
will be extended with pathways from WikiPathways database to create
PPI-pathway networks. Finally, MCL (Markov Clustering) network
clustering algorithm will be applied to get clusters within the network.

## Setup

Installing and loading required libraries

## PPI network for selected biopsy location

``` r
#Obtain data from step 5:
setwd('..')
## Select a location to create network (options; ileum or rectum)
location <- "ileum"

#read dataset to be processed for ileum or rectum biopsy location
if (location == "ileum") {
  deg <- read.delim("5-extract-overlapped_genes/output/DEG.overlapped_ileum")
  print("Selected location is ileum")
}else if(location == "rectum"){ 
  deg <- read.delim("5-extract-overlapped_genes/output/DEG.overlapped_rectum")
  print("Selected location is rectum")
}
```

    ## [1] "Selected location is ileum"

``` r
#filter out genes that does not have ENTREZ ID 
deg <- deg %>% tidyr:: drop_na(ENTREZ)
#check that cytoscape is connected
cytoscapePing()
#close session before starting
closeSession(FALSE)
networkName = paste0("PPI_network_",location)
#create a PPI network using overlapped DEGs between CD and UC
#first take the deg input list as a query
x <- readr::format_csv(as.data.frame(deg$ENTREZ), col_names=F, escape = "double", eol =",")
#then use the below function to convert a command string into a CyREST query URL, executes a GET request, 
#and parses the result content into an R list object. Same as commandsGET
commandsRun(paste0('string protein query cutoff=0.7', ' newNetName=',networkName, ' query=',x,' limit=0'))
```

    ## [1] "Loaded network 'STRING network - PPI_network_ileum' with 189 nodes and 32 edges"

``` r
#get proteins (nodes) from the constructed network
proteins <- RCy3::getTableColumns(columns=c("query term", "display name"))
#get edges from the network
ppis     <- RCy3::getTableColumns(table="edge", columns=c("name"))
#split extracted edge information into source-target format
ppis     <- data.frame(do.call('rbind', strsplit(as.character(ppis$name),' (pp) ',fixed=TRUE)))
#merge obtained nodes and edges to get entrez IDs for each source genes 
ppis.2   <- merge(ppis, proteins, by.x="X1", by.y="display name", all.x=T)
#change column names
colnames(ppis.2) <- c("s", "t", "source")
#merge again to add entrez IDs of target genes 
ppis.3   <- merge(ppis.2, proteins, by.x="t", by.y="display name", all.x=T)
colnames(ppis.3)[4] <-"target"
#ppi3 stores interaction between all proteins so add new column represeting type of interaction
ppis.3$interaction <- "PPI"
#add col names to protein
colnames(proteins) <- c("id","label")
proteins$type <- "protein"

###############get all pathways from WIKIPATHWAYS #################
#below code should be performed first to handle the ssl certificate error uploading pathways 
options(RCurlOptions = list(cainfo = paste0( tempdir() , "/cacert.pem" ), ssl.verifypeer = FALSE))
wp.hs.gmt <- rWikiPathways::downloadPathwayArchive(organism="Homo sapiens", format = "gmt")
#wp.hs.gmt <- "wikipathways-20220210-gmt-Homo_sapiens.gmt"
#Now that we have the latest GMT file for human pathways, 
#all wp and gene information stored in wp2gene object
wp2gene   <- rWikiPathways::readPathwayGMT(wp.hs.gmt)
#filter out  pathways that does not consist of any differentially expressed genes 
wp2gene.filtered <- wp2gene [wp2gene$gene %in% deg$ENTREZ,]

#change column names 
colnames(wp2gene.filtered)[3] <- c("source")
colnames(wp2gene.filtered)[5] <- c("target")
#add new column for representing interaction type
wp2gene.filtered$interaction <- "Pathway-Gene"

#store only wp information 
pwy.filtered <- unique( wp2gene [wp2gene$gene %in% deg$ENTREZ,c(1,3)])
colnames(pwy.filtered) <- c("label", "id")
pwy.filtered$type <- "pathway"
colnames(pwy.filtered) <- c("label","id", "type")

#get genes 
genes <- unique(deg[,c(1,2)])
genes$type <- "gene"
colnames(genes) <- c("id","label","type")
genes$id <- as.character(genes$id)
#genes and pathways are separate nodes and they need to be merged
nodes.ppi <- dplyr::bind_rows(genes,pwy.filtered)
rownames(nodes.ppi) <- NULL
edges.ppi <- unique(dplyr::bind_rows(ppis.3[,c(3,4,5)], wp2gene.filtered[,c(3,5,6)]))
rownames(edges.ppi) <- NULL
```

## Create PPI-pathway network for selected biopsy location

``` r
#create a network name 
networkName = paste0("PPI_Pathway_Network_",location)

###########Create PPI-pathway network###
RCy3::createNetworkFromDataFrames(nodes= nodes.ppi, edges = edges.ppi, title=networkName, collection=location)
```

    ## networkSUID 
    ##      172820

``` r
RCy3::loadTableData(nodes.ppi, data.key.column = "label", table="node", table.key.column = "label")
```

    ## [1] "Success: Data loaded in defaultnode table"

``` r
RCy3::loadTableData(deg, data.key.column = "ENTREZ", table.key.column = "id")
```

    ## [1] "Success: Data loaded in defaultnode table"

``` r
###########Visual style#################
RCy3::copyVisualStyle("default","ppi")#Create a new visual style (ppi) by copying a specified style (default)
RCy3::setNodeLabelMapping("label", style.name="ppi")
```

    ## NULL

``` r
RCy3::lockNodeDimensions(TRUE, style.name="ppi")#Set a boolean value to have node width and height fixed to a single size value.

#threshold is set based of differential expressed gene criteria
data.values<-c(-0.58,0,0.58) 
#red-blue color schema chosen
node.colors <- c(brewer.pal(length(data.values), "RdBu"))
#nodes are split to show both log2fc values for both diseases 
RCy3::setNodeCustomHeatMapChart(c("log2FC_CD","log2FC_UC"), slot = 2, style.name = "ppi", colors = c("#CC3300","#FFFFFF","#6699FF","#CCCCCC"))
RCy3::setVisualStyle("ppi")
```

    ##                 message 
    ## "Visual Style applied."

``` r
# Saving output
if(!dir.exists("output")) dir.create("output")
outputName = paste0 ("output/PPI_Pathway_Network_",location,".png")
png.file <- file.path(getwd(), outputName)
exportImage(png.file,'PNG', zoom = 500)
```

    ##                                                                                                                                                                              file 
    ## "C:\\Users\\dedePC\\Desktop\\FNS_CLOUD\\GITHUB_codes\\Transcriptomics_Metabolomics_Analysis\\transcriptomics_analysis\\6-network_analysis\\output\\PPI_Pathway_Network_ileum.png"

## Clustering obtained networks

``` r
#we will continue with the same session used for ppi-pathway networks
#to check cytoscape is connected
cytoscapePing()

networkName = paste0 ("PPI_Pathway_Network_",location)
#to get network name of the rectum 
networkSuid = getNetworkSuid(networkName)
setCurrentNetwork(network=networkSuid)
#create cluster command
clustermaker <- paste("cluster mcl createGroups=TRUE showUI=TRUE network=SUID:",networkSuid, sep="")
#run the command in cytoscape
res <- commandsGET(clustermaker)
#total number of clusters 
cat("Total number of clusters for",location, as.numeric(gsub("Clusters: ", "", res[1])))
```

    ## Total number of clusters for ileum 68

``` r
#change pathway node visualization
pathways <- RCy3::selectNodes(nodes="pathway", by.col = "type")
RCy3::setNodeColorBypass(node.names = pathways$nodes, "#D3D3D3")
RCy3::setNodeBorderWidthBypass(node.names = pathways$nodes, 10)

#export image
outputName = paste0 ("output/PPI_Pathway_Network_",location,"_clustered",".png")
png.file <- file.path(getwd(), outputName)
exportImage(png.file,'PNG', zoom = 500)
```

    ##                                                                                                                                                                                        file 
    ## "C:\\Users\\dedePC\\Desktop\\FNS_CLOUD\\GITHUB_codes\\Transcriptomics_Metabolomics_Analysis\\transcriptomics_analysis\\6-network_analysis\\output\\PPI_Pathway_Network_ileum_clustered.png"

``` r
#save session
#cys.file <- file.path(getwd(), "output/PPI_Pathway_Network_",location,"_clustered",".cys")
#saveSession(cys.file) 
```

## Last, we create a Jupyter notebook from this script

``` r
#Jupyter Notebook file
if(!"devtools" %in% installed.packages()) BiocManager::install("devtools")
devtools::install_github("mkearney/rmd2jupyter", force=TRUE)
```

    ## 
    ## * checking for file 'C:\Users\dedePC\AppData\Local\Temp\RtmpgLo4e5\remotes4e1419657980\mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION' ... OK
    ## * preparing 'rmd2jupyter':
    ## * checking DESCRIPTION meta-information ... OK
    ## * checking for LF line-endings in source and make files and shell scripts
    ## * checking for empty or unneeded directories
    ## Omitted 'LazyData' from DESCRIPTION
    ## * building 'rmd2jupyter_0.1.0.tar.gz'
    ## 

``` r
library(devtools)
library(rmd2jupyter)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
rmd2jupyter("Network_analysis.Rmd")
```
