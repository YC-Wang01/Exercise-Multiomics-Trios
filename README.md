# FinlandSports V2.0: ATM Multi-omics Analysis Pipeline

## Project Overview
This repository contains the complete analytical pipeline for the "FinlandSports (ATM)" multi-omics research project. Our study aims to explore the molecular mechanisms of exercise-induced systemic responses and inter-tissue crosstalk (Adipose-Muscle-Serum) through the integration of transcriptomics, proteomics, metabolomics, and DNA methylation data.

## Directory Structure
The repository is systematically organized into data, scripts, and results:

```text
C:.
│  .gitignore
│  .Rhistory
│  Exercise-Multiomics-Trios.Rproj
│  project_structure.txt
│  README.md
│  
├─00_Raw_Data (Not tracked via Git due to size limits)
├─01_Clean_Data (Not tracked via Git due to size limits)
│      
├─02_Scripts
│  ├─S0_Data_Preparation
│  │      00_Clinical_Purification.R
│  │      01_Omics_QC_and_Normalization.R
│  │      
│  ├─S1_Main_Figures
│  │      01_Fig1D_CV_Analysis.R
│  │      02_Fig1D_PCA_Distance.R
│  │      03_Fig1H_Clinical_Phenotypes.R
│  │      04_Fig2AB_Differential_Analysis.R
│  │      05_Fig2CD_S2_GSEA_Network.R
│  │      06_Fig2EF_Mechanisms.R
│  │      
│  └─S2_Supplementary_Figures
│          01_FigS1_PCA_Limma.R
│          02_FigS2_Cyclebar.R
│          03_FigS2_Clinical_GSEA_Crosstalk.R
│          
└─03_Results
    ├─Fig_1 (Baseline Landscape, PCA & Phenotypes)
    ├─Fig_2 (Exercise Response, GSEA, Network & Mechanisms)
    ├─Fig_S1 (PCA Trajectories & Significant Proportions)
    ├─Fig_S2 (Clinical GSEA Crosstalk & Cyclebars)
    └─Summary_Stats (Limma Results & Loadings)