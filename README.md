# BCSIR-LUAD-Study

# Epigenetic Regulation of lncRNAs in Lung Adenocarcinoma (LUAD)

This repository contains the R scripts used for the integrative analysis of DNA methylation and gene expression data to identify epigenetically deregulated lncRNAs in Lung Adenocarcinoma (LUAD), as described in our study. 

## Analysis Workflow

The analysis is organized into the following steps:

1.  **Differential Expression Analysis**: `01_Differential_Expression_Analysis.R`
    - Identifies differentially expressed genes (DEGs) using `DESeq2`.
2.  **Visualization**: `02_Visualization_Volcano.R`
    - Generates volcano plots for visualizing differential expression results.
3.  **Differential Methylation Analysis**: `03_Differential_Methylation_Analysis.R`
    - Identifies differentially methylated CpG sites (DMCs) using `ChAMP`. Methylation preprocessing and DMP calling were performed using the ChAMP pipeline on IDAT files. Although the ChAMP workflow uses the IlluminaHumanMethylation450kanno.ilmn12.hg19 package for array-level compatibility, final CpG coordinate mapping and promoter assignment in the integrative analysis were performed separately using the hg38-based HM450_hg38_manifest_gencode_v36 annotation resource.
4.  **Integration & Correlation**: `04_Integration_and_Correlation.R`
    - Performs Spearman correlation analysis between promoter methylation and lncRNA expression.
5.  **Survival Analysis**: `05_Survival_Analysis_KaplanMeier.R`
    - Evaluates the prognostic value of candidate lncRNAs using Kaplan-Meier curves.
6.  **Prognostic Modeling**: `06_Prognostic_Modeling_LASSO_ROC.R`
    - Constructs a prognostic signature using LASSO Cox regression and validates it using ROC analysis.
7.  **Network Analysis**: `07_Network_Construction_Correlation.R`
    - Analyzes correlations within the lncRNA-miRNA-mRNA regulatory network.
  
#Code Disclaimer:
These scripts are a programmatic representation of the methodology outlined in the manuscript. Because they were compiled post-analysis, they are provided primarily for methodological transparency and inspection of the analytical framework, and may require structural adjustments depending on individual R environments and data formats.
