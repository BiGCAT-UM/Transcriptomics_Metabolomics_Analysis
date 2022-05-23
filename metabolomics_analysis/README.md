The Metabolomics workflow consist out of three steps in total:

7. [Data preprocessing](https://github.com/BiGCAT-UM/Transcriptomics_Metabolomics_Analysis/tree/master/metabolomics_analysis/7-metabolite_data_preprocessing)
- Filters features not annotated with an HMDB identifier
- Filter features with the name “redundant ion”
- Harmonizes HMDB IDs to HMDB-4.0 structure
- Removes special characters (e.g. an asterisks, *) in IDs.
- Filters samples on the first time point where metabolomics data was collected as baseline.

8. [Significantly Changed Metabolite Analysis](https://github.com/BiGCAT-UM/Transcriptomics_Metabolomics_Analysis/tree/master/metabolomics_analysis/8-significantly_changed_metabolites_analysis)

- Removes missing values for annotated metabolites with >50% missing value 
- Normalization by cube root, square root, log2 or log10 (effect tested by a Shapiro-Wilk normality test)
- Calculates Fold Change (FC) and p-value (based on a t-test).
- Sets threshold values for significantly changed metabolites.
- Handles duplicate IDs for significantly changed metabolites, by calculating the average of the FC and p-values was calculated
- Visualizes metabolite data in a Volcano Plot

9. [Metabolite Pathway Analysis (based on HMDB IDs)](https://github.com/BiGCAT-UM/Transcriptomics_Metabolomics_Analysis/tree/master/metabolomics_analysis/9-metabolite_pathway_analysis)

- Retrieves significantly changed metabolites as HMDB IDs
- Compared sign. changed metabolites against the content of the pathway in the WikiPathways SPARQL endpoint (filtering out pathways from the Reactome database) 
- Performs Enrichment analysis of the metabolomics data by Over Representation Analysis (ORA) using a Fisher's exact test (p-value) by means of a hypergeometric density calculation.
- Sorts results based on the highest number of matching individual significantly changed metabolites, lowest p-value, and highest number of proteins for each pathway.

The data from step 9 can be combined with the data from the transcriptomics analysis (step 2), after which multi-omics visualization is possible (step 10 and 11).
