## Introduction

In this workflow, we link the metabolites of interest to pathway data
from WikiPathways, based on their HMDB identifier.

``` r
# Obtain Working Directory for step 8 to find processed data
setwd('..')

#Obtain data from step 8
mSet_CD <- read.csv("8-significantly_changed_metabolites_analysis/output/mbxData_CD.csv", na.strings=c("", "NA"))
mSet_UC <- read.csv("8-significantly_changed_metabolites_analysis/output/mbxData_UC.csv", na.strings=c("", "NA"))

## Select a disorder to analyse (options; CD or UC)
disorder <- "CD"

if (disorder == "CD") {
  mSet = mSet_CD 
  print("Selected disorder is Crohn's disease")}else if(disorder == "UC"){ 
    mSet = mSet_UC
    print("Selected disorder is Ulcerative Colitis")}else{print("Disorder not Recognised")}
```

    ## [1] "Selected disorder is Crohn's disease"

## Find pathways based on relevant IDs column

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
 #below code should be performed first to handle the ssl certificate error
options(RCurlOptions = list(cainfo = paste0( tempdir() , "/cacert.pem" ), ssl.verifypeer = FALSE))
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
select distinct ?pathwayRes (str(?wpid) as ?pathway) (str(?title) as ?pathwayTitle) (count(distinct ?hmdbMetabolite) AS ?HMDBsInPWs) 
(GROUP_CONCAT(DISTINCT fn:substring(?hmdbMetabolite,30);separator=' ') AS ?includedHMDBs)
where {
VALUES ?hmdbMetabolite {"
item2 = "}
 
 ?datanode  a wp:Metabolite ;          
            wp:bdbHmdb  ?hmdbMetabolite ;
            dcterms:isPartOf ?pathwayRes .
            
 ?pathwayRes a wp:Pathway ;
             wp:organismName 'Homo sapiens' ; 
            dcterms:identifier ?wpid ;
            dc:title ?title .
            
  #?pathwayRes wp:ontologyTag cur:Reactome_Approved . 
  ?pathwayRes wp:ontologyTag cur:AnalysisCollection .           
}
ORDER BY DESC(?HMDBsInPWs)"
query_CombinePWs <- paste(item1,string_HMDB,item2)
remove(item1, item2)

results_CombinePWs <- SPARQL(endpointwp,query_CombinePWs,curl_args=list(useragent=R.version.string))
showresults_CombinePWs <- results_CombinePWs$results
remove(query_CombinePWs,results_CombinePWs)

##Top pathway cuttoff threshold (can be defined by users)
pathway_cutoff_metabolites = 10

#Keep and print table within threshold (if less than threshold are found, print only those)
if(nrow(showresults_CombinePWs) < pathway_cutoff_metabolites){
print(showresults_CombinePWs[1:nrow(showresults_CombinePWs),c(2:4)])
}else{
  #delete rows below threshold
  showresults_CombinePWs <- showresults_CombinePWs[-c((pathway_cutoff_metabolites+1):nrow(showresults_CombinePWs)),]
  print(showresults_CombinePWs[1:5,c(2:4)])}
```

    ##   pathway                                                  pathwayTitle
    ## 1  WP3604                                  Biochemical pathways: part I
    ## 2  WP2525 Trans-sulfuration, one-carbon metabolism and related pathways
    ## 3  WP4723                        Omega-3 / omega-6 fatty acid synthesis
    ## 4  WP3925                                         Amino acid metabolism
    ## 5   WP661                                           Glucose homeostasis
    ##   HMDBsInPWs
    ## 1         17
    ## 2          8
    ## 3          7
    ## 4          7
    ## 5          6

``` r
remove(cleaned_string_HMDB, list_Relevant_HMDB_IDs, query_HMDBs)
```

Retrieve additional relevant pathway data, for example protein and gene
data

``` r
string_WP_IDs_list <- paste0("'", showresults_CombinePWs$pathway ,"'")
string_WP_IDs <- paste(c(string_WP_IDs_list), collapse=' ' )

item1 = "
PREFIX ch: <https://identifiers.org/hmdb/>
PREFIX cur: <http://vocabularies.wikipathways.org/wp#Curation:>
select distinct ?pathwayRes (count(distinct ?metaboliteDatanode) AS ?TotalMetabolitesinPW) (GROUP_CONCAT(DISTINCT fn:substring(?hgnc,37);separator=' ') AS ?Proteins) (count(distinct ?hgnc) AS ?ProteinsInPWs)
where {
VALUES ?wpid {
"
item2 = "
}
 
 ?metaboliteDatanode    a wp:Metabolite ;
                       dcterms:isPartOf ?pathwayRes .
 
 ?pathwayRes a wp:Pathway ;
             wp:organismName 'Homo sapiens' ; 
             dcterms:identifier ?wpid ;
             dc:title ?title ;
             wp:ontologyTag cur:AnalysisCollection .   
  OPTIONAL{         
 ?datanode2 wp:bdbHgncSymbol ?hgnc ;
            dcterms:isPartOf ?pathwayRes .
  }
        
}
"
query_CombinePWs_gene <- paste(item1,string_WP_IDs,item2)
remove(item1, item2)

results_CombinePWs_gene <- SPARQL(endpointwp,query_CombinePWs_gene,curl_args=list(useragent=R.version.string))
showresults_CombinePWs_gene <- results_CombinePWs_gene$results
remove(query_CombinePWs_gene,results_CombinePWs_gene)

##Merge the two dataframes together:
showresults_CombinePW_data <- merge(x = showresults_CombinePWs, y = showresults_CombinePWs_gene, by = "pathwayRes", all.x = TRUE)
##Reorder data to fit previous format:
showresults_CombinePW_data <- showresults_CombinePW_data[, c(1:4, 6:8, 5)]
```

## Calculate the ORA score for each pathway, using the Fishers exact test.

``` r
##Based on: https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/

#Create a dataframe to store the required numbers in.
Contingency_table <- data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("WP.ID", "x", "m", "n", "k"))))
counter = 1
for (i in 1:nrow(showresults_CombinePW_data)) {
   Contingency_table[counter,1] <- (showresults_CombinePW_data[i,2]) #WP.ID
   Contingency_table[counter,2] <- (showresults_CombinePW_data[i,4]) ##x <- (number4) #Total differentially changed metabolites, also in a PW. (HMDBsInPWs)
   Contingency_table[counter,3] <- (showresults_CombinePW_data[i,5]) ##m <- (number) #Total Metabolites in PW (TotalMetabolitesinPW)
   Contingency_table[counter,4] <- (length(unique(mSet[,1])) - showresults_CombinePW_data[i,4]) ##n <- (number2) #Total Metabolites measured not in PW (DISTINCT all_HMDB - HMDBsInPWs)
   Contingency_table[counter,5] <- length(unique(vector_HMDB)) ##k <- (number3) #Total differentially changed metabolites. (DISTINCT vector_HMDB)

   counter <- counter + 1
}

# Calculate hypergeometric density p-value for all pathways.
i <- 1:nrow(Contingency_table)
probabilities <- dhyper(Contingency_table[i,2], Contingency_table[i,3], Contingency_table[i,4], Contingency_table[i,5], log = FALSE)

pathwayAnalysis_results <- cbind(showresults_CombinePW_data[, c(2:4)], probabilities, showresults_CombinePW_data[, c(6,7)])
colnames(pathwayAnalysis_results)[5] <- "HGNCs"
colnames(pathwayAnalysis_results)[6] <- "ProteinsInPWs"

##Sort PW results based on 1. highest amount of #HMDBs in PW, 2. lowest p-values,  3. highest amouny of proteins in PW (which might be relevant for transcriptomics analysis later)
pathwayAnalysis_results_sorted <- pathwayAnalysis_results[  with(pathwayAnalysis_results, order(-HMDBsInPWs, probabilities, -ProteinsInPWs)),]

print(pathwayAnalysis_results_sorted[1:5,])
```

    ##   pathway                                                  pathwayTitle
    ## 4  WP3604                                  Biochemical pathways: part I
    ## 2  WP2525 Trans-sulfuration, one-carbon metabolism and related pathways
    ## 5  WP3925                                         Amino acid metabolism
    ## 6  WP4723                        Omega-3 / omega-6 fatty acid synthesis
    ## 1    WP15                                Selenium micronutrient network
    ##   HMDBsInPWs probabilities
    ## 4         17  3.544334e-14
    ## 2          8  7.972883e-02
    ## 5          7  1.109534e-04
    ## 6          7  1.551589e-01
    ## 1          6  2.236586e-05
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                HGNCs
    ## 4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
    ## 2                                                                                                               AGXT2 AHCY AHCYL1 AHCYL2 AMT BAAT BCAT1 BCAT2 BHMT BHMT2 CBS CBSL CDO1 CEPT1 CHDH CHKA CHKB CHPT1 CSAD CTH DHFR DHFR2 DMGDH DNM1 DNMT3A DNMT3B DNMT3L ETNK1 ETNK2 GAD1 GAD2 GCLC GCLM GNMT GPX1 GPX2 GPX3 GPX4 GPX5 GPX6 GPX7 GSR GSS MAT1A MAT2A MAT2B MTHFD1 MTHFD1L MTHFD2 MTHFD2L MTHFR MTR PCYT1A PCYT1B PCYT2 PEMT PHGDH PLD1 PSAT1 PSPH SARDH SHMT1 SHMT2 SOD1 SOD2 SOD3 TYMS
    ## 5                  ACAA1 ACADM ACLY ACO2 ACSS1 ADH1C ADH4 ADH5 ADH7 ALDH18A1 ALDH1A1 ALDH7A1 AOC3 ARG1 ARG2 ASNS ASS1 AUH BCAT1 BHMT CAD CBS CPS1 CS CTH DBH DDC DLD DLST EHHADH EPRS1 FAH FARSB FH FTCD G6PC2 GCLM GLS GLUD1 GLUL GOT1 GOT2 GPT2 GSR GSS HADH HAL HDC HIBADH HIBCH HMGCL HMGCS2 HNMT IARS1 IDH1 LARS2 LDHA MAOA MARS2 MCCC1 MDH1 MDH2 MMUT MPST OAT ODC1 OGDH OTC P4HA2 PC PCK1 PDHA1 PDHX PDK4 PKM PNMT PPM1L PYCR1 RARS1 SDHA SDS SMS SRM SUCLG1 TAT TDO2 TH TPH1 TPO VARS1 WARS1
    ## 6                                                                                                                                                                                                                                                                                                                                                                                                  ACOT1 ACOT2 ACOX1 ACOX3 ACSL1 ACSL3 ACSL4 ELOVL2 ELOVL5 FADS1 FADS2 PLA2G4A PLA2G4B PLA2G5 PLA2G6
    ## 1 ABCA1 ALB ALOX15B ALOX5 ALOX5AP APOA1 APOB CAT CBS CCL2 CRP CTH DIO1 DIO2 DIO3 F2 F7 FGA FGB FGG FLAD1 GGT1 GPX1 GPX2 GPX3 GPX4 GPX6 GSR HBA1 HBB ICAM1 IFNG IL1B IL6 INS INSR KMO KYNU LDLR MPO MSRB1 MTHFR MTR NFKB1 NFKB2 PLAT PLG PNPO PRDX1 PRDX2 PRDX3 PRDX4 PRDX5 PTGS1 PTGS2 RELA RFK SAA1 SAA2 SAA3P SAA4 SCARB1 SELENOF SELENOH SELENOI SELENOK SELENOM SELENON SELENOO SELENOP SELENOS SELENOT SELENOV SELENOW SEPHS2 SERPINA3 SERPINE1 SOD1 SOD2 SOD3 TNF TXN TXNRD1 TXNRD2 TXNRD3 XDH
    ##   ProteinsInPWs
    ## 4             0
    ## 2            67
    ## 5            91
    ## 6            15
    ## 1            86

## Export the pathway data:

``` r
##Save the data file
nameDataFile <- paste0("output/mbxPWdata_", disorder ,".csv")
write.table(pathwayAnalysis_results_sorted, nameDataFile, sep =",", row.names = FALSE)
```

## Print significantly changed metabolites which were not in a pathway, by ID and name:

``` r
##Find Missing Biomarkers (not part of any Human pathway model)
item1 = "PREFIX ch: <https://identifiers.org/hmdb/>
SELECT DISTINCT ?HMDBMetabolite WHERE {
  VALUES ?HMDBMetabolite {"
item2 = "}
  ?pathwayRes  a wp:Pathway ;
                wp:organismName 'Homo sapiens' .
  
  ?metabolite   a wp:Metabolite ;
                dcterms:identifier ?id ;
                dcterms:isPartOf ?pathwayRes .
  ?metabolite wp:bdbHmdb ?HMDBMetabolite.
}"
queryMissingBiomarkers <- paste(item1,string_HMDB,item2)
remove(item1,item2)
resultsMissingBiomarkers <- SPARQL(endpointwp,queryMissingBiomarkers,curl_args=list(useragent=R.version.string))
listMissingBiomarkers <- c(resultsMissingBiomarkers$results) #safe results as list for comparison.
remove(queryMissingBiomarkers,resultsMissingBiomarkers)
HMDBs_inPWs <- gsub("[<https://identifiers.org/hmdb/>]", "", listMissingBiomarkers) #HMDB IDs IRI cleanup
intersectingHMDB <- setdiff(vector_HMDB, HMDBs_inPWs)

string_intersectingHMDB <- paste(c(intersectingHMDB), collapse=', ' )

#Find names for missing Biomarkers based on HMDB ID (to help with data understanding and curation)
missingNames <- list()
for (j in 1:length(intersectingHMDB)){
  for (i in 1:nrow(mSet)){
    if(!is.na(mSet[i,5]) & mSet[i,5] == intersectingHMDB[j]){
       missingNames[j] <- mSet[i,6]
      }
    else{next}
  }
}
remove(i,j)
#Save list on one string for reporting purposes
string_missingNames <- do.call(paste, c(as.list(missingNames), sep = ", "))
#Print relevant information:
if(length(intersectingHMDB) == 0 ){print("All relevant biomarkers are in a pathway!")} else{
  print(paste0("For the disorder ", disorder, ", ", length(intersectingHMDB), " biomarkers are not in a pathway; with the following HMDB IDs: " , string_intersectingHMDB, "; with the following Database names: ", string_missingNames))}
```

    ## [1] "For the disorder CD, 50 biomarkers are not in a pathway; with the following HMDB IDs: HMDB0000479, HMDB0000610, HMDB0000779, HMDB0000792, HMDB0000848, HMDB0000885, HMDB0000932, HMDB0001906, HMDB0002250, HMDB0002815, HMDB0004161, HMDB0005015, HMDB0005065, HMDB0005066, HMDB0005462, HMDB0005476, HMDB0006726, HMDB0006731, HMDB0006733, HMDB0007199, HMDB0007871, HMDB0007970, HMDB0008006, HMDB0008038, HMDB0010169, HMDB0010370, HMDB0010379, HMDB0010383, HMDB0010384, HMDB0010391, HMDB0010393, HMDB0010395, HMDB0010404, HMDB0010407, HMDB0011208, HMDB0011212, HMDB0011241, HMDB0011243, HMDB0011252, HMDB0011310, HMDB0011503, HMDB0011507, HMDB0011520, HMDB0012097, HMDB0012104, HMDB0013122, HMDB0013287, HMDB0013325, HMDB0015070, HMDB0030180; with the following Database names: 3-methylhistidine, C18:2 CE, phenylalanine, sebacate, C18 carnitine, C16:0 CE, tauro-alpha-muricholate/tauro-beta-muricholate, aminoisobutyric acid/GABA, C12 carnitine, C18:1 LPC, urobilin, gabapentin, C18:1 carnitine, C14 carnitine, C56:7 TAG, C58:10 TAG, C20:4 CE, C20:5 CE, C22:6 CE, C38:5 DAG, C32:0 PC, C34:0 PC, C34:3 PC, C36:1 PC, C16:0 SM, C18:3 CE, C16:0 LPC, C16:1 LPC, C18:0 LPC, C20:1 LPC, C20:3 LPC, C20:4 LPC, C22:6 LPC, C16:1 LPC plasmalogen, C34:1 PC plasmalogen, C34:4 PC plasmalogen, C36:1 PC plasmalogen, C36:2 PC plasmalogen, C38:4 PC plasmalogen, C36:4 PC plasmalogen, C16:0 LPE, C18:2 LPE, C22:0 LPE, C14:0 SM, C22:1 SM, C18:1 LPC plasmalogen, N6,N6-dimethyllysine, C10:2 carnitine, oxymetazoline, crustecdysone"

## Print session info:

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
    ## [1] SPARQL_1.16     RCurl_1.98-1.10 XML_3.99-0.13  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_4.2.2  fastmap_1.1.1   cli_3.6.0       tools_4.2.2    
    ##  [5] htmltools_0.5.4 rstudioapi_0.14 yaml_2.3.7      rmarkdown_2.20 
    ##  [9] knitr_1.42      xfun_0.37       digest_0.6.31   bitops_1.0-7   
    ## [13] rlang_1.0.6     evaluate_0.20

## Last, we create a Jupyter notebook file from this script:

``` r
#Jupyter Notebook file
if(!"devtools" %in% installed.packages()) BiocManager::install("devtools")
devtools::install_github("mkearney/rmd2jupyter", force=TRUE)
```

    ## 
    ## ── R CMD build ─────────────────────────────────────────────────────────────────
    ##          checking for file 'C:\Users\duygu\AppData\Local\Temp\RtmpQTnLcr\remotes49dc34c469f0\mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION' ...  ✔  checking for file 'C:\Users\duygu\AppData\Local\Temp\RtmpQTnLcr\remotes49dc34c469f0\mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION'
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
rmd2jupyter("metabolomics_pathway_analysis.Rmd")

##Clean up data
remove(counter, i, probabilities, string_HMDB, vector_HMDB, Contingency_table, pathwayAnalysis_results, showresults_CombinePWs, showresultsMetadata)
```
