## A semi-automated workflow for functional analysis of transcriptomic and metabolomic data applied to understand inflammatory bowel diseases

This workflow, developed in R markdown files as well as and Jupyter notebook files, includes differential gene expression analysis, statistical analysis of metabolomics data, as well as pathway enrichment analysis for both transcriptomics and metabolomics data followed by integration of this data through network analysis to identify disease-related processes and visualization of multi-omics data. A publicly available (https://ibdmdb.org/) gut-transcriptomic and stool-metabolome dataset of the gut microbial ecosystem in inflammatory bowel diseases was used to test the proposed workflow. 

Transcriptomics analysis:  
[1-Data preprocessing](https://github.com/BiGCAT-UM/Transcriptomics_Metabolomics_Analysis/tree/master/transcriptomics_analysis/1-data_preprocessing)<br /> 
[2-Differential gene expression analysis](/transcriptomics_analysis/2-differential_gene_expression_analysis/)<br />
[3-Gene Pathway Analysis (ORA)](/transcriptomics_analysis/3-pathway_analysis/)<br />
[4-Heatmap creation](/transcriptomics_analysis/4-create_heatmap/)<br />
[5-Overlapped genes extraction](/transcriptomics_analysis/5-extract-overlapped_genes/)<br />
[6-Network analysis](/transcriptomics_analysis/6-network_analysis)<br />
Metabolomics analysis:  
[7-Data preprocessing](metabolomics_analysis/7-metabolite_data_preprocessing/)<br />
[8-Significantly changed metabolites analysis](metabolomics_analysis/8-significantly_changed_metabolites_analysis/)<br />
[9-Metabolite Pathway Analysis (ORA)](metabolomics_analysis/9-metabolite_pathway_analysis/)<br />
Multi-omics visualization<br />
[10-Identifier mapping](visualization_multiomics/10-identifier_mapping/)<br />
[11-Visualizaiton of multi-omics](visualization_multiomics/11-visualization/)<br />


The setup of this project has been tested with:
- OS Windows 10, R-studio 2021.09.02, R 4.1.3.
- OS Linux (Debian), R-studio 2022.02.2, R 4.2.

The workflow is an example of how to bring together different software tools and methods to analyze transcriptomics and metabolomics data to reveal the underlying mechanism behind IBD disease as a use case study.

![workflow](https://user-images.githubusercontent.com/65600609/169603139-1732ce3c-fd87-4c8e-9925-9ce7134653d8.png)

### Acknowledgment

This research was undertaken by Maastricht University (UM, Netherlands), a beneficiary in FNS-Cloud, which has received funding from the European Union’s Horizon 2020 Research and Innovation programme (H2020-EU.3.2.2.3. – A sustainable and competitive agri-food industry) under Grant Agreement No. 863059. 
