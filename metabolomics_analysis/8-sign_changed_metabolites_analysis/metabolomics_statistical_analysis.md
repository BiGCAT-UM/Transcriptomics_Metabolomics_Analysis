## Introduction

In this workflow, we will apply statistical analysis on metabolomics
data and link the metabolites of interest to pathway data from
WikiPathways.

## First we locate the metabolomics data from step-7 (preprocessing).

``` r
setwd('..')

#Obtain data from step 7
mSet_CD <- read.csv("7-metabolite_data_preprocessing/output/mbxDataCD_nonIBD.csv", na.strings=c("", "NA"))
mSet_UC <- read.csv("7-metabolite_data_preprocessing/output/mbxDataUC_nonIBD.csv", na.strings=c("", "NA"))

## Select a disorder to analyse (options; CD or UC)
disorder <- "CD"

if (disorder == "CD") {
  mSet = mSet_CD 
  print("Selected disorder is Crohn's disease")}else if(disorder == "UC"){ 
    mSet = mSet_UC
  print("Selected disorder is Ulcerative Colitis")}else{print("Disorder not Recognised")
}
```

    ## [1] "Selected disorder is Crohn's disease"

## Second, we perform data extraction from the file, and process the data

``` r
#Merge column headers: disorder_patientID
names(mSet) <- paste(mSet [1, ], names(mSet), sep = "_")
mSet <- mSet [-1,]

##DATA PROCESSING:

#Remove metabolites with < 50% data (located in columns 8-553.
columns <- ncol(mSet)
rowsData <- nrow(mSet)
removeLines <- rowSums(is.na(mSet[,3-columns])) #HERE WE DELETE A COLUMN WHICH SHOULD NOT BE LIKE THAT -it should be like-> mSet[,3:columns]
fifty_percent <- floor((columns)/2)

mSet_MissingDataCounted <- cbind(mSet, removeLines)
mSet_NoMissingData <- subset(mSet_MissingDataCounted, removeLines <= fifty_percent)
#Remove last column for further processing.
mSet_NoMissingData <- subset(mSet_NoMissingData, select=-c(removeLines))

#Convert intensity data to numeric values                         
mSet_NoMissingData[, c(3:columns)] <- apply(mSet_NoMissingData[, c(3:columns)],2, function(x) as.numeric(as.character(x)))

##normalization (see https://doi.org/10.1177%2F1469066720918446 and https://www.statology.org/transform-data-in-r/)
##Users can select different transformation styles here:
transformation <- "log_2" #log_2 transformation on intensity data selected for IBD dataset, options are: cube_root, square_root, log_2, log_10

if(transformation == "cube_root"){
    mSet_transformed <- cbind(mSet_NoMissingData[,c(1,2)], mSet_NoMissingData[,3:columns]^(1/3))
}else if(transformation == "square_root"){
    mSet_transformed <- cbind(mSet_NoMissingData[,c(1,2)], mSet_NoMissingData[,3:columns]^(1/2))
}else if(transformation == "log_2"){
    mSet_transformed <- cbind(mSet_NoMissingData[,c(1,2)], log2(mSet_NoMissingData[,3:columns]))
}else if(transformation == "log_10"){
    mSet_transformed <- cbind(mSet_NoMissingData[,c(1,2)], log10(mSet_NoMissingData[,3:columns]))
}else{print("Warning: name for transformation not recognized")}

colnames(mSet_transformed)[1]="HMDB.ID"
colnames(mSet_transformed)[2]="Compound.Name"

## Visualize the data after the transformation (for one sample to get an idea of suitability of transformation:
#create histogram for original distribution for first column with data
hist(mSet_NoMissingData[,3], col='steelblue', main='Original')
```

![](metabolomics_statistical_analysis_files/figure-markdown_github/data_import-1.png)

``` r
#create histogram for log-transformed distribution 
hist(mSet_transformed[,3], col='coral2', main=transformation)
```

![](metabolomics_statistical_analysis_files/figure-markdown_github/data_import-2.png)

``` r
## Testing if the transformation creates a normally distributed dataset (alpha >= 0.05)
##Calculate all Shapiro values for raw and transformed data:
mSet_NoMissingData_Shapiro <- lapply(mSet_NoMissingData[,3:columns], shapiro.test)
mSet_transformed_Shapiro <- lapply(mSet_transformed[,3:columns], shapiro.test)

#Obtain the p-values for raw and transformed data
mSet_NoMissingData_Shapiro_pvalues <- do.call(rbind, mSet_NoMissingData_Shapiro)
mSet_transformed_Shapiro_pvalues <- do.call(rbind, mSet_transformed_Shapiro)

## Count how often the p-value is above 0.05, to obtain an estimate of achieved normality due to transformation
mSet_NoMissingData_Shapiro_pvalues_sum <- sum(mSet_NoMissingData_Shapiro_pvalues[,2] >= 0.05, na.rm=TRUE)
mSet_transformed_Shapiro_pvalues_sum <- sum(mSet_transformed_Shapiro_pvalues[,2] >= 0.05, na.rm=TRUE)

eighty_percent <- floor(((columns)/10)*8)

#Print relevant information:
if(mSet_transformed_Shapiro_pvalues_sum[1] > eighty_percent ){paste0("Data after ", transformation ," transformation seems to follow a normal distribution for more then 80% of your data")} else{
  print("Advised to select a different data transformation procedure")}
```

    ## [1] "Data after log_2 transformation seems to follow a normal distribution for more then 80% of your data"

``` r
remove(mSet_HMDB_NEW, mSet_MissingDataCounted, mSet_NoMissingData, mSet_NoMissingData_Shapiro, mSet_NoMissingData_Shapiro_pvalues, mSet_transformed_Shapiro, mSet_transformed_Shapiro_pvalues, eighty_percent, fifty_percent, mSet_NoMissingData_Shapiro_pvalues_sum, mSet_transformed_Shapiro_pvalues_sum, removeLines, columns, rowsData)
```

## Calculate logFC and p-value

``` r
#Create backup of data
mSet_transformed.b <- mSet_transformed
##Order columns based on disease abbreviation (located in column names); adding the control group nonIBD as last columns:
if (disorder == "CD") {
mSet_FINAL <- mSet_transformed[ , order(names(mSet_transformed))]}else{mSet_FINAL <- mSet_transformed[ , order(names(mSet_transformed), decreasing=TRUE)]}

##Move Name column back to column 2
columnNumber <- which(colnames(mSet_FINAL)=="Compound.Name")
mSet_FINAL <- mSet_FINAL[,c(columnNumber,1:ncol(mSet_FINAL)-1)]## WHILE APPLYING THIS STEP ONE COLOUMN WAS DELETED ->  FOR CD "nonIBD_MSM9VZMM" is deleted
#we need to change it as mSet_FINAL <- mSet_FINAL[,c(columnNumber,1:ncol(mSet_FINAL))
mSet_FINAL <- mSet_FINAL [,-(columnNumber+1)]

##Move HMDB column back to the start
columnNumber <- which(colnames(mSet_FINAL)=="HMDB.ID")
mSet_FINAL <- mSet_FINAL[,c(columnNumber,1:ncol(mSet_FINAL)-1)]
mSet_FINAL <- mSet_FINAL [,-(columnNumber+1)]

#Find relevant columns per group.
library(stringr)
columns_disorders <- sum(str_count(colnames(mSet_FINAL), disorder))
end_Disorders <- (columns_disorders+2) ##add two, to compensate for HMDB.ID and Compound.Names column.

##calculate logFC for 2 groups (CD vs control, UC vs control), ignoring missing values (NAs) when calculating the mean.  
disease = apply(mSet_FINAL[,3:end_Disorders], 1, mean, na.rm=TRUE)
control_IBD = apply(mSet_FINAL[,(end_Disorders+1):ncol(mSet_FINAL)], 1, mean, na.rm=TRUE)

#because the metabolomics data is already log2 transformed, we need to take the difference between the means (iso dividing the means over one another), since log2 Fold Change or log2 Ratio == log2(condition / control). Note: if the transformation step applied is cube_root or square_root, one needs to divide control over disease for this step!
if(transformation == "log_2" | transformation == "log10"){
    foldchange_disorder <-  disease - control_IBD
}else{
  foldchange_disorder <- (disease /control_IBD )}

##ADD HMDB column at start, add fold change columns.
mSet_AnalysisReady <- cbind(mSet_FINAL$HMDB.ID, mSet_FINAL$Compound.Name, foldchange_disorder)
##Rename first column and second column
colnames(mSet_AnalysisReady)[1] <- "HMDB_ID"
colnames(mSet_AnalysisReady)[2] <- "Compound_Name"

##Calculate p-value for two groups based on t-test (comparing control to disease).
##general function to store p-values for multiple rows:
ttest_mSet <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}
p_values_disorder <- apply(mSet_FINAL, 1, ttest_mSet, grp1 = c(3:end_Disorders), grp2 = c((end_Disorders+1):ncol(mSet_FINAL)))

##Add p_values column to analysis dataset:
mSet_AnalysisReady <- cbind(mSet_AnalysisReady, p_values_disorder)

#Convert logFC and p-values columns to numeric values            
mSet_AnalysisReady <- as.data.frame(mSet_AnalysisReady)
mSet_AnalysisReady[ , c(3,4)] <- apply(mSet_AnalysisReady[ , c(3,4)], 2, function(x) as.numeric(as.character(x)))

remove(mSet_transformed, columns_disorders, columnNumber, end_Disorders, foldchange_disorder, p_values_disorder, control_IBD, disease, ttest_mSet, mSet_FINAL)
```

## Volcano plot

The next step will visualize relevant metabolites in a Volcano plot,
with thresholds defined for the log2FCs and p-values. Since we are
dealing with non-targeted LC_MS data, we are applying stringent criteria
here (p-value below 0.05, and \|log2FC\| \>= 1.5)

``` r
##Inspired by: https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html
if(!"ggplot2" %in% installed.packages()){install.packages("ggplot2")}
library('ggplot2')

##Define the thresholds for log2 (Fold Change) and p-values
#For cut-off value uncertainties, see https://doi.org/10.1039/C6AN01342B .
log2FC_min <- -1
log2FC_max <- 1
p_value_threshold <- 0.05

##Create column with HMDB_IDs, only if the data is relevant
mSet_AnalysisReady$relevant_labels <- mSet_AnalysisReady$HMDB_ID
mSet_AnalysisReady$relevant_labels[!((mSet_AnalysisReady$foldchange_disorder <= log2FC_min | mSet_AnalysisReady$foldchange_disorder >= log2FC_max) &  mSet_AnalysisReady$p_values_disorder <= p_value_threshold)] <- NA

##Duplication issues:
if(!"dplyr" %in% installed.packages()){install.packages("dplyr")}
library(dplyr)

##Check for duplicate HMDB.IDs, however with different compound name (example: HMDB0011130, called LPE, LPE-A (both relevant), and LPE-B (not relevant) & HMDB0000252, called sphingosine isomer 1, 2 and 3.)
mSet_AnalysisReady_Duplicates <- mSet_AnalysisReady %>% group_by(HMDB_ID) %>% distinct(Compound_Name, .keep_all = TRUE)
#Remove all rows with relevant data, which will be added later by their average in case of duplicate HMDB.IDs:
mSet_AnalysisReady_Duplicates <- subset(mSet_AnalysisReady_Duplicates, is.na(mSet_AnalysisReady_Duplicates$relevant_labels))
mSet_AnalysisReady_Duplicates <- subset(mSet_AnalysisReady_Duplicates[,1:4])

##Check for duplicates HMDB.IDs, calculate the average for these if both significant (example: HMDB0010384).
mSet_AnalysisReady_FC <- mSet_AnalysisReady %>% group_by(HMDB_ID, Compound_Name) %>% filter(!is.na(relevant_labels)) %>% summarize(foldchange_disorder=mean(foldchange_disorder)) 
mSet_AnalysisReady_p <- mSet_AnalysisReady %>% group_by(HMDB_ID, Compound_Name)  %>% filter(!is.na(relevant_labels)) %>% summarize(p_values_disorder=mean(p_values_disorder))

#Merge FC and p-value data together based on HMDB.ID
mSet_AnalysisReady_FCandp <-
left_join(mSet_AnalysisReady_FC, mSet_AnalysisReady_p, by=c("HMDB_ID", "Compound_Name")) %>%
  rowwise()  
  colnames(mSet_AnalysisReady_FCandp)[3] <- "foldchange_disorder"
  colnames(mSet_AnalysisReady_FCandp)[4] <- "p_values_disorder"

##Combine datasets again, before continuing to next steps:
mSet_AnalysisFinal <- rbind(mSet_AnalysisReady_Duplicates, mSet_AnalysisReady_FCandp)

##Create column with HMDB_IDs again, only if the data is relevant, for visualization in Volcano Plot
mSet_AnalysisFinal$relevant_ids <- mSet_AnalysisFinal$HMDB_ID
mSet_AnalysisFinal$relevant_ids[!((mSet_AnalysisFinal$foldchange_disorder <= log2FC_min | mSet_AnalysisFinal$foldchange_disorder >= log2FC_max) &  mSet_AnalysisFinal$p_values_disorder <= p_value_threshold)] <- NA

##Create another column with Compound names, only if the data is relevant, for visualization in Volcano Plot
mSet_AnalysisFinal$relevant_labels <- mSet_AnalysisFinal$Compound_Name
mSet_AnalysisFinal$relevant_labels[!((mSet_AnalysisFinal$foldchange_disorder <= log2FC_min | mSet_AnalysisFinal$foldchange_disorder >= log2FC_max) &  mSet_AnalysisFinal$p_values_disorder <= p_value_threshold)] <- NA

##Select visualization of HMDB-IDs(relevant_ids), or chemical compound names(relevant_labels)
selectViz <- "relevant_labels"

# Library needed to space out the labels
if(!"ggrepel" %in% installed.packages()){install.packages("ggrepel")}
library(ggrepel)

if(selectViz == "relevant_labels"){
##volcanoPlot_Disorder 
volcanoPlot_disorder <- ggplot(data=mSet_AnalysisFinal, aes(x=foldchange_disorder, y=-log10(p_values_disorder), label=relevant_labels)) + geom_point() + theme_minimal() + geom_text_repel()
}else if(selectViz == "relevant_ids"){
##volcanoPlot_Disorder 
volcanoPlot_disorder <- ggplot(data=mSet_AnalysisFinal, aes(x=foldchange_disorder, y=-log10(p_values_disorder), label=relevant_ids)) + geom_point() + theme_minimal() + geom_text_repel()
}else{print("Column name not recognized for label visualization.")}

## Add vertical lines for FoldChange and P-value thresholds:
volcanoPlot_disorder <- volcanoPlot_disorder + geom_vline(xintercept=c(log2FC_min, log2FC_max), col="blue") +
    geom_hline(yintercept=-log10(p_value_threshold), col="red") + theme(plot.background = element_rect(fill = "white"))

if (disorder == "UC") {disorderName <- "Ulcerative Colitis (UC)"}else{disorderName <- "Crohn's disease (CD)"}

titleVolcano <- paste0("Volcano plot of ", transformation, " transformed data for ", disorderName )
verticalAxisTitle <- paste0(transformation, " Fold Change, ", disorderName, " versus control ")

## Add title and update axis labels:
volcanoPlot_disorder <- volcanoPlot_disorder + ggtitle(titleVolcano) + labs(y = "-log10(p-value)", x = verticalAxisTitle)

# Show the Volcano plot in the notebook output:
volcanoPlot_disorder
```

![](metabolomics_statistical_analysis_files/figure-markdown_github/volcano_plot-1.png)
\## Export the data and Volcano Plot:

``` r
##Save the data file
nameDataFile <- paste0("output/mbxData_", disorder ,".csv")
write.table(mSet_AnalysisFinal, nameDataFile, sep =",", row.names = FALSE)

if(!"svglite" %in% installed.packages()){install.packages("svglite")}
library('svglite')

##Save the Volcano plot:
imageType <- "png" ##Options are: svg, png, eps, ps, tex, pdf, jpeg, tiff, png, bmp, svg or wmf
nameVolcano <- paste0("output/", disorder, "_", selectViz, "_VolcanoPlot_absLogFC_", log2FC_max, "_pValue_", p_value_threshold, ".", imageType)

ggsave(nameVolcano)
```

## Last, we create a Jupyter notebook and markdown file from this script

``` r
#Jupyter Notebook file
if(!"devtools" %in% installed.packages()) BiocManager::install("devtools")
devtools::install_github("mkearney/rmd2jupyter", force=TRUE)
```

    ## 
    ## ── R CMD build ─────────────────────────────────────────────────────────────────
    ##       ✔  checking for file 'C:\Users\duygu\AppData\Local\Temp\Rtmpyu29C9\remotes620025544512\mkearney-rmd2jupyter-d2bd2aa/DESCRIPTION'
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
rmd2jupyter("metabolomics_statistical_analysis.Rmd")

##Clean up R-studio environment
remove(list_Relevant_HMDB_IDs, mSet_AnalysisFinal, mSet_AnalysisReady, mSet_AnalysisReady_Duplicates, mSet_AnalysisReady_FC, mSet_AnalysisReady_FCandp, mSet_AnalysisReady_p, mSet_transformed.b, volcanoPlot_disorder, disorderName, imageType, log2FC_max, log2FC_min, nameDataFile, nameVolcano, p_value_threshold, selectViz, titleVolcano, transformation, verticalAxisTitle)
```
