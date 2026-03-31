# FinlandSports: Multi-Omics Landscape of Exercise Intervention (V2.0)
[![DOI](https://zenodo.org/badge/1195134340.svg)](https://doi.org/10.5281/zenodo.19343366)
[![R Version](https://img.shields.io/badge/R-v4.4.0-blue.svg)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Reproducibility](https://img.shields.io/badge/Reproducibility-100%25-brightgreen.svg)]()

## 📖 Project Overview
This repository contains the complete, standardized computational pipeline (V2.0) for the **FinlandSports** multi-omics project. This study investigates the systemic biological responses to exercise intervention across multiple tissues (Adipose, Muscle) and biofluids (Serum) using transcriptomic, proteomic, and metabolomic profiling, alongside deep clinical phenotyping.

The V2.0 pipeline was entirely refactored to ensure strict reproducibility, completely decoupling raw data processing from downstream visualization, and employing dynamic in-memory delta-value engines.

## ✨ Key Pipeline Features
* **Absolute Reproducibility**: All scripts operate on dynamic relative paths originating from the unified `01_Clean_Data` directory.
* **Clinical Data Engine**: Subject mapping relies entirely on a single ground-truth metadata file (`Clinical_Master_Strict.csv`), eliminating sample misalignment.
* **Robust Statistical Modeling**: Incorporates comprehensive Two-way ANOVA and `limma`-based linear models, strictly adjusted for Age, Gender, and subject-specific baseline covariates (e.g., blood cell composition).
* **Publication-Ready Visualization**: Automated generation of scaled, high-resolution vector graphics (PDFs) with customized academic palettes (Nature NPG).

## 📂 Repository Structure
```text
.
├── 00_Raw_Data/                  # Original raw matrices and clinical questionnaires
├── 01_Clean_Data/                # Standardized, quality-controlled CSVs and Sample Sheets
├── 02_Scripts/                   # Modularized R scripts for the entire analysis
│   ├── S0_Data_Preparation/      # Data purification and baseline normalization
│   ├── S1_Main_Figures/          # Scripts generating Main Figs 1 through 6
│   └── S2_Supplementary_Figures/ # Scripts for Supplementary Figs
├── 03_Results/                   # Auto-generated outputs (PDFs, Tables)
│   ├── Fig_1/ ... Fig_6/         # Figure-specific visual components and stats
│   ├── Fig_S1/ ... Fig_S4/       # Supplementary plots
│   └── Summary_Stats/            # Global statistical matrices
├── Exercise-Multiomics-Trios.Rproj
└── README.md
