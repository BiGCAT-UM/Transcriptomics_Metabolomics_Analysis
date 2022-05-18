## Introduction

In this workflow, we link the metabolites of interest to pathway data
from WikiPathways, based on their HMDB identifier.

``` r
# Obtain Working Directory for step 8 to find processed data
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()
```

    ## [1] "/home/deniseslenter/Documents/GitHub/Transcriptomics_Metabolomics_Analysis/metabolomics_analysis/9-metabolite_pathway_analysis"

``` r
setwd('..')
work_DIR <- getwd()

#Obtain data from step 8
mSet_CD <- read.csv("8-significantly_changed_metabolites_analysis/output/mbxData_CD.csv", na.strings=c("", "NA"))
mSet_UC <- read.csv("8-significantly_changed_metabolites_analysis/output/mbxData_UC.csv", na.strings=c("", "NA"))

# Set Working Directory back to current folder
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
work_DIR <- getwd()

## Select a disorder to analyse (options; CD or UC)
disorder <- "CD"

if (disorder == "CD") {
  mSet = mSet_CD 
  print("Selected disorder is Crohn's disease")}else if(disorder == "UC"){ 
    mSet = mSet_UC
    print("Selected disorder is Ulcerative Colitis")}else{print("Disorder not Recognised")}
```

    ## [1] "Selected disorder is Crohn's disease"

Find pathways based on relevant labels column

``` r
if(!"SPARQL" %in% installed.packages()){
  install.packages("SPARQL")
}
library(SPARQL)
##Connect to Endpoint WikiPathways
endpointwp <- "https://sparql.wikipathways.org/sparql"
## 1. Query metadata:
queryMetadata <-
"SELECT DISTINCT ?dataset (str(?titleLit) as ?title) ?date ?license 
WHERE {
   ?dataset a void:Dataset ;
   dcterms:title ?titleLit ;
   dcterms:license ?license ;
   pav:createdOn ?date .
 }"
resultsMetadata <- SPARQL(endpointwp,queryMetadata,curl_args=list(useragent=R.version.string))
showresultsMetadata <- resultsMetadata$results
remove(queryMetadata, resultsMetadata)


## Create a list of HMDB IDs according to filtering criteria from step 8.
list_Relevant_HMDB_IDs <- list(mSet$relevant_ids)
vector_HMDB <- unlist(list_Relevant_HMDB_IDs) #convert list to array, for traversing the data to a SPARQL query later on
vector_HMDB <- vector_HMDB[!is.na(vector_HMDB)]
##Add the HMDb prefix IRI in front of all IDs.
query_HMDBs <- paste("ch:", vector_HMDB, sep="")
##Merge the individual entries in the vector into one string, separated by a space
string_HMDB <- paste(c(query_HMDBs), collapse=' ' )

##TODO: add column with nr. of Metabolites in PW (to calculate PW impact)

#For now, filter out Reactome PWs due to visualization issues in Cytoscape.
item1 = "PREFIX ch: <https://identifiers.org/hmdb/>
PREFIX cur: <http://vocabularies.wikipathways.org/wp#Curation:>
select distinct ?pathwayRes (str(?wpid) as ?pathway) (str(?title) as ?pathwayTitle) (count(distinct ?hmdbMetabolite) AS ?HMDBsInPWs) (count(distinct ?metaboliteDatanode) AS ?TotalMetabolitesinPW) (GROUP_CONCAT(DISTINCT fn:substring(?hgnc,37);separator=' ') AS ?Proteins) (count(distinct ?hgnc) AS ?ProteinsInPWs) (GROUP_CONCAT(DISTINCT fn:substring(?hmdbMetabolite,30);separator=' ') AS ?includedHMDBs)
where {
VALUES ?hmdbMetabolite {"
item2 = "}
 
 ?metaboliteDatanode    a wp:Metabolite ;
                        dcterms:isPartOf ?pathwayRes .
 
 ?datanode  a wp:Metabolite ;          
            wp:bdbHmdb  ?hmdbMetabolite ;
            dcterms:isPartOf ?pathwayRes .
 ?pathwayRes a wp:Pathway ;
             wp:organismName 'Homo sapiens' ; 
            dcterms:identifier ?wpid ;
            dc:title ?title .
            
 ?datanode2 wp:bdbHgncSymbol ?hgnc ;
            dcterms:isPartOf ?pathwayRes .
            
  #?pathwayRes wp:ontologyTag cur:Reactome_Approved . 
  ?pathwayRes wp:ontologyTag cur:AnalysisCollection .           
}
ORDER BY DESC(?HMDBsInPWs)"
query_CombinePWs <- paste(item1,string_HMDB,item2)
remove(item1, item2)

results_CombinePWs <- SPARQL(endpointwp,query_CombinePWs,curl_args=list(useragent=R.version.string))
showresults_CombinePWs <- results_CombinePWs$results
remove(query_CombinePWs,results_CombinePWs)
#Print table with first 5 relevant pathways (if less than 5 are found, print only those)
if(nrow(showresults_CombinePWs) < 5){
print(showresults_CombinePWs[1:nrow(showresults_CombinePWs),c(2:5)])
}else{print(showresults_CombinePWs[1:5,c(2:5)])}
```

    ##   pathway                               pathwayTitle HMDBsInPWs
    ## 1  WP3925                      Amino acid metabolism          7
    ## 2  WP4723     Omega-3 / omega-6 fatty acid synthesis          6
    ## 3    WP15             Selenium micronutrient network          6
    ## 4  WP3940 One-carbon metabolism and related pathways          6
    ## 5   WP661                        Glucose homeostasis          6
    ##   TotalMetabolitesinPW
    ## 1                  108
    ## 2                   38
    ## 3                  110
    ## 4                   41
    ## 5                   21

``` r
remove(cleaned_string_HMDB, list_Relevant_HMDB_IDs, query_HMDBs)
```

Calculate the ORA score for each pathway, using the Fishers exact test.

``` r
##Based on: https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
##TODO: finish this section!

#Create a dataframe to store the required numbers in.
Contingency_table <- data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("WP.ID", "x", "m", "n", "k"))))
counter = 1
for (i in 1:nrow(showresults_CombinePWs)) {
   Contingency_table[counter,1] <- (showresults_CombinePWs[i,2]) #WP.ID
   Contingency_table[counter,2] <- (showresults_CombinePWs[i,4]) ##x <- (number4) #Total differentially changed metabolites, also in a PW. (HMDBsInPWs)
   Contingency_table[counter,3] <- (showresults_CombinePWs[i,5]) ##m <- (number) #Total Metabolites in PW (TotalMetabolitesinPW)
   Contingency_table[counter,4] <- (length(unique(mSet[,1])) - showresults_CombinePWs[i,4]) ##n <- (number2) #Total Metabolites measured not in PW (DISTINCT all_HMDB - HMDBsInPWs)
   Contingency_table[counter,5] <- length(unique(vector_HMDB)) ##k <- (number3) #Total differentially changed metabolites. (DISTINCT vector_HMDB)

   counter <- counter + 1
}

# Calculate hypergeometric density p-value for all pathways.
i <- 1:nrow(Contingency_table)
probabilities <- dhyper(Contingency_table[i,2], Contingency_table[i,3], Contingency_table[i,4], Contingency_table[i,5], log = FALSE)

pathwayAnalysis_results <- cbind(showresults_CombinePWs[, c(2:4)], probabilities, showresults_CombinePWs[, c(6,7)])
colnames(pathwayAnalysis_results)[5] <- "HGNCs"
colnames(pathwayAnalysis_results)[6] <- "ProteinsInPWs"

##Sort PW results based on 1. highest amount of #HMDBs in PW, 2. lowest p-values,  3. highest amouny of proteins in PW (which might be relevant for transcriptomics analysis later)
pathwayAnalysis_results_sorted <- pathwayAnalysis_results[  with(pathwayAnalysis_results, order(-HMDBsInPWs, probabilities, -ProteinsInPWs)),]

print(pathwayAnalysis_results_sorted[1:5,])
```

    ##   pathway                               pathwayTitle HMDBsInPWs probabilities
    ## 1  WP3925                      Amino acid metabolism          7  1.109534e-04
    ## 3    WP15             Selenium micronutrient network          6  2.236586e-05
    ## 4  WP3940 One-carbon metabolism and related pathways          6  1.079858e-01
    ## 2  WP4723     Omega-3 / omega-6 fatty acid synthesis          6  1.322064e-01
    ## 5   WP661                        Glucose homeostasis          6  1.428674e-01
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                HGNCs
    ## 1                  ACAA1 ACADM ACLY ACO2 ACSS1 ADH1C ADH4 ADH5 ADH7 ALDH18A1 ALDH1A1 ALDH7A1 AOC3 ARG1 ARG2 ASNS ASS1 AUH BCAT1 BHMT CAD CBS CPS1 CS CTH DBH DDC DLD DLST EHHADH EPRS1 FAH FARSB FH FTCD G6PC2 GCLM GLS GLUD1 GLUL GOT1 GOT2 GPT2 GSR GSS HADH HAL HDC HIBADH HIBCH HMGCL HMGCS2 HNMT IARS1 IDH1 LARS2 LDHA MAOA MARS2 MCCC1 MDH1 MDH2 MMUT MPST OAT ODC1 OGDH OTC P4HA2 PC PCK1 PDHA1 PDHX PDK4 PKM PNMT PPM1L PYCR1 RARS1 SDHA SDS SMS SRM SUCLG1 TAT TDO2 TH TPH1 TPO VARS1 WARS1
    ## 3 ABCA1 ALB ALOX15B ALOX5 ALOX5AP APOA1 APOB CAT CBS CCL2 CRP CTH DIO1 DIO2 DIO3 F2 F7 FGA FGB FGG FLAD1 GGT1 GPX1 GPX2 GPX3 GPX4 GPX6 GSR HBA1 HBB ICAM1 IFNG IL1B IL6 INS INSR KMO KYNU LDLR MPO MSRB1 MTHFR MTR NFKB1 NFKB2 PLAT PLG PNPO PRDX1 PRDX2 PRDX3 PRDX4 PRDX5 PTGS1 PTGS2 RELA RFK SAA1 SAA2 SAA3P SAA4 SCARB1 SELENOF SELENOH SELENOI SELENOK SELENOM SELENON SELENOO SELENOP SELENOS SELENOT SELENOV SELENOW SEPHS2 SERPINA3 SERPINE1 SOD1 SOD2 SOD3 TNF TXN TXNRD1 TXNRD2 TXNRD3 XDH
    ## 4                                                                                                                                                                                                           AGXT2 AHCYL1 BAAT BCAT1 BCAT2 BHMT BHMT2 CBSL CDO1 CEPT1 CHDH CHKA CHKB CHPT1 CSAD CTH DHFR2 DMGDH DNM1 DNMT3A ETNK1 ETNK2 GAD1 GAD2 GCLC GCLM GNMT GPX1 GPX2 GPX3 GPX4 GPX5 GPX6 GPX7 GSR GSS MAT1A MAT2A MTHFR MTR PCYT1A PCYT1B PCYT2 PEMT PLD1 SARDH SHMT1 SHMT2 SOD1 SOD2 SOD3 TYMS
    ## 2                                                                                                                                                                                                                                                                                                                                                                                                        ACOT2 ACOX1 ACOX3 ACSL1 ACSL3 ACSL4 ELOVL2 ELOVL5 FADS1 FADS2 PLA2G4A PLA2G4B PLA2G5 PLA2G6
    ## 5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                INS
    ##   ProteinsInPWs
    ## 1            91
    ## 3            86
    ## 4            52
    ## 2            14
    ## 5             1

Export the pathway data:

``` r
##Save the data file
nameDataFile <- paste0("output/mbxPWdata_", disorder ,".csv")
write.table(pathwayAnalysis_results_sorted, nameDataFile, sep =",", row.names = FALSE)
```

Print significantly changed metabolites which were not in a pathway, by
ID and name:

``` r
#TODO
```

### Last, we create a Jupyter notebook and markdown file from this script

``` r
#Jupyter Notebook file
if(!"devtools" %in% installed.packages()) BiocManager::install("devtools")
devtools::install_github("mkearney/rmd2jupyter", force=TRUE)
```

    ## 
    ## * checking for file ‘/tmp/RtmpIgYty7/remotes70f7468ca1a2/mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION’ ... OK
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
rmd2jupyter("metabolomics_pathway_analysis.Rmd")

#markdown_file <- "metabolomics_pathway_analysis.md"
#if (file.exists(markdown_file)) {
#   unlink(markdown_file, recursive=TRUE)#first delete the existing one
# }
#If this next line trows an error, build the md file with knittr manual selection (file, Knit document, or ctrl+shift+k -keyboard shortcut).
#rmarkdown::render("metabolomics_pathway_analysis.Rmd", "md_document")

##Clean up data
remove(counter, i, probabilities, string_HMDB, vector_HMDB, Contingency_table, pathwayAnalysis_results, showresults_CombinePWs, showresultsMetadata)
```
