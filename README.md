# DP-RST Workload: Supplementary Materials and Reproducibility

This repository contains the complete code and supplementary materials necessary to reproduce the results presented in the manuscript: **"Bayesian Nonparametric Clustering of Spatial Transcriptomics Data with Complex Tissue Geometry"**.

It includes scripts for data simulation, preprocessing pipelines, model implementation, and downstream analysis for both standard and high-resolution spatial transcriptomics platforms.

## Repository Structure

The repository is organized into the following primary components, reflecting the analyses conducted in the study:

### 1. Simulation Studies
This directory contains scripts to generate synthetic data and benchmark the DP-RST model against competing methods across four distinct topological scenarios:
* **Swiss-Roll:** A non-convex spiral domain designed to test manifold learning capabilities.
* **U-shape:** A nested non-convex structure.
* **Curl:** A highly convoluted domain with skewed-normal expression profiles.
* **Ellipse:** Disjoint spatial domains mimicking separated tissue sections.

### 2. Real Data Applications
This section includes the preprocessing and analysis workflows for the biological datasets analyzed in the manuscript.

#### Main Analysis: Mouse Intestine
* **Preprocessing:** Scripts for quality control, normalization, and PCA on the Swiss-roll colon tissue sample.
* **Clustering:** Implementation of DP-RST to identify the five histological layers.
* **Downstream Analysis:** Code for Differential Gene Expression (DGE) analysis to identify marker genes.

#### Benchmark Datasets
Scripts to reproduce the performance benchmarks on seven additional datasets spanning various technologies:
* **10x Visium:** Human DLPFC, Breast Cancer, and Prostate Cancer.
* **Imaging-based:** STARmap (Mouse Visual Cortex) and osmFISH (Mouse Somatosensory Cortex).

### 3. Scalable Inference (High-Resolution)
This directory contains the implementation of the **Consensus Monte Carlo (ConMC)** framework designed for large-scale datasets. It includes:
* **Data Partitioning:** Scripts to split massive datasets into spatial shards (e.g., 28 shards for Lung, 100 shards for Colon).
* **Parallel Inference:** Code to run parallel MCMC chains on individual shards.
* **Aggregation:** Implementation of the Optimal Transport alignment strategy to merge local results into a global consensus.
* **Target Datasets:** Workflows for the 10x Visium HD Colon Cancer dataset (~230,000 observations) and 10x Xenium Lung dataset (~50,000 cells).

### 4. Sensitivity Analysis
Code to reproduce the sensitivity analyses presented in the Supplementary Materials, assessing model stability with respect to:
* The Dirichlet Process concentration parameter ($\alpha$).
* The initial number of spatial clusters ($K$).

## Dependencies

The analysis relies on the `DP.RST` package. Please refer to the main package repository for installation instructions. Additional dependencies for data handling and visualization include `Seurat`, `SingleCellExperiment`, and standard spatial analysis libraries.
