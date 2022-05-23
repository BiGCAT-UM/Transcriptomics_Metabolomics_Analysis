The Metabolomics workflow consist out of three steps in total:

7. [Data preprocessing](https://github.com/BiGCAT-UM/Transcriptomics_Metabolomics_Analysis/tree/master/metabolomics_analysis/7-metabolite_data_preprocessing)
- Filters features not annotated with an HMDB identifier
- Filter features with the name “redundant ion”
- Harmonizes HMDB IDs to HMDB-4.0 structure
- Removes special characters (e.g. an asterisks, *) in IDs.
- Filters samples on the first time point where metabolomics data was collected as baseline.

9. [Significantly Changed Metabolite Analysis](https://github.com/BiGCAT-UM/Transcriptomics_Metabolomics_Analysis/tree/master/metabolomics_analysis/8-significantly_changed_metabolites_analysis)


11. [Metabolite Pathway Analysis (based on HMDB IDs)](https://github.com/BiGCAT-UM/Transcriptomics_Metabolomics_Analysis/tree/master/metabolomics_analysis/9-metabolite_pathway_analysis)

The data from step 9 can be combined with the data from the transcriptomics analysis (step 2), after which multi-omics visualization is possible (step 10 and 11).
