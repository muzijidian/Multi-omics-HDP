# Integration of metagenomics and metabolomics for preeclampsia and gestational hypertension: longitudinal multi-omics analyses

This repository contains the R and Python source code for the computational analysis, multi-omics integration, predictive modeling, and visualization presented in our study on hypertensive disorders of pregnancy (HDP).

## Project Overview
This study integrates longitudinal maternal gut metagenomics, mycobiome, and metabolomics data to explore their interplay in preeclampsia (PE) and gestational hypertension (GH). Our goal is to uncover disease-specific multi-omics signatures and develop predictive models for early pregnancy risk assessment.

## Repository Structure
The scripts are organized according to the analytical workflow and corresponding main figures in the manuscript.

### 1. Microbiome Diversity and Covariate Analysis
* **`Fig2ab.R`**: Analyzes and visualizes longitudinal alpha diversity (richness and Shannon index) of the maternal gut bacteriome and mycobiome.
* **`Fig2cde.R`**: Performs beta diversity analyses (Bray-Curtis distance, PCoA, and PERMANOVA) to evaluate longitudinal compositional shifts and temporal stability.
* **`Fig2f.R`**: Calculates and visualizes the variance in the gut microbiome explained by maternal baseline characteristics, lifestyle, clinical history, and dietary covariates.

### 2. Differential Signatures and Multi-Omics Integration
* **`Fig3a.R`**: Generates comprehensive heatmaps of differentially abundant microbial taxa (bacteria/fungi) across different gestational stages.
* **`Fig3f.R`**: Constructs forest plots for subgroup analyses and interaction testing of key microbial features.
* **`Fig3g.py`**: Visualizes correlations between specific microbial species and clinical indicators (e.g., blood pressure).
* **`Fig4a.py`**: Generates a Circos plot detailing the longitudinal associations of differential metabolites with PE and GH.
* **`Fig4b.R`**: Analyzes the longitudinal dynamics of differential metabolites, visualizing intersections and class distributions.
* **`Fig4c.R`**: Evaluates the cross-cohort consistency of identified metabolic signatures.
* **`Fig4de.R`**: Quantifies and visualizes the variance of metabolites explained by the gut microbiome (bacteria/fungi) and host metadata.

### 3. Predictive Modeling and Validation
* **`Figure5_prediction.R`**: Constructs random forest models to predict PE and GH outcomes, implementing repeated cross-validation, upsampling, and external validation routines.
* **`Figure5.R`**: Visualizes the performance of predictive models.