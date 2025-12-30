**DP‑RST Workload: Supplementary Materials and Reproducibility**

This repository accompanies the manuscript **“Bayesian Non‑parametric Clustering of Spatial Transcriptomics Data with Complex Tissue Geometry.”** It contains the code, synthetic and real data, and analysis scripts necessary to reproduce the results and figures reported in the paper. The directory layout mirrors the structure of the supplementary materials and separates simulation experiments from real data analyses, scalability experiments and sensitivity studies.

## Top‑level layout

| Directory / file               | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| ------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `.gitattributes`, `.gitignore` | Standard repository configuration files.                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| **`README.md`**                | High‑level description of the repository (this document).                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| **`Simulations/`**             | Scripts and data for generating synthetic spatial transcriptomics datasets and benchmarking DP‑RST against competing methods across four topologies: *Swiss‑Roll*, *U‑shape*, *Curl* and *Ellipse*. Each subfolder contains an R script to generate the geometry and expression profiles, PDF summaries of spatial and density features, optional CSV/zip files of domain contours, and a `*_methods` sub‑directory containing method‑specific code for DP‑RST, BAST, BayesSpace, DR.SC and other competitors. |
| **`Datasets Results/`**        | Real data analysis workflows. Subdirectories correspond to individual datasets—**Brain**, **Breast**, **Gut** and **Lung_Xenium**—and contain `Methods/` (R scripts for preprocessing, clustering with DP‑RST/BAST and competing methods, and downstream analyses) and `Results/` (intermediate and final outputs such as cluster assignments, differential gene expression tables and diagnostics). Additional scripts compute Normalised Mutual Information (NMI) scores across methods.                     |
| **`ConsensusMCM/`**            | Implementation of the *Consensus Monte Carlo* (ConMC) framework for scaling DP‑RST to very large datasets. The directory contains compressed input data (`Colon.zip`, `Lung.zip`), R scripts for preparing spatial shards and running DP‑RST/BAST on each shard (`Colon_scripts/`, `Lung_scripts/`), and a `ReferenceBasedAlignment/` folder with Python and R code for aligning shard‑level results via optimal transport and aggregating them into a global consensus.                                       |
| **`Figures/`**                 | R scripts used to recreate the figures in the paper. Files such as `3D_PCA_plot – Figure 1.R`, `Change&Hyper moves plot – Figure 2.R` and `Swiss_Roll_Results – Figures 1 and 4.R` generate the main manuscript figures. A subfolder `PCs_figures/` contains pre‑generated PNG files of principal component expression/density plots used in the paper. Additional scripts and `.RData` objects support figure creation.                                                                                       |
| **`Sensitivity/`**             | Code and data for sensitivity analyses. Scripts such as `SR_sim_SensAlpha.R` and `SR_sim_SensK.R` vary the Dirichlet Process concentration parameter and the number of initial clusters and save results as `.RData` and PDF plots.                                                                                                                                                                                                                                                                            |
| **`Extra_fun/`**               | A collection of helper functions written in R. These functions compute and plot minimum‑spanning‑tree (MST) cluster assignments, merge MST components and create PDF summaries.                                                                                                                                                                                                                                                                                                                                |
| **`Legacy BAST/`**             | Legacy implementation of the Bayesian Adaptive Spatial Tree (BAST) clustering algorithm. The scripts `BASTFun_2.R`, `ComplexDomainFun.R` and `FEMFun.R` are retained for reproducibility.                                                                                                                                                                                                                                                                                                                      |
| **`profiling.R`**              | Example script showing how to profile the DP‑RST implementation using the `profvis` package on a Swiss‑Roll simulation. It loads simulated data, constructs an initial MST and cluster assignments, sets hyper‑parameters, runs a short MCMC chain and saves an interactive profiling report.                                                                                                                                                                                                                  |

## Simulation studies

The *Simulations* directory allows researchers to reproduce the synthetic experiments reported in the paper. For each topology (**Curl**, **Ellipse**, **Swiss‑Roll** and **U‑shape**) there are:

* **Data generation scripts:** these create spatial coordinates, boundaries and synthetic expression profiles with user‑specified parameters. For example, `Curl_DataGeneration.R` produces a convoluted domain with skewed‑normal gene expression and saves the coordinates, true cluster labels and boundaries.
* **Feature summaries:** PDF files such as `Curl_SpatialFeatures.pdf` and `Curl_DensityFeatures.pdf` summarise the spatial layout and gene‑expression densities for each synthetic dataset.
* **Method benchmarks:** the `*_methods` sub‑directories provide R scripts to run DP‑RST and competing clustering methods (BAST, BayesSpace, DR.SC, etc.), with versions tailored to 3‑PC and 10‑PC feature sets. Some variants implement batch‑wise inference to handle larger simulated datasets.

A reference image (`Simulations_reference.png`) compares the four topologies at a glance.

## Real data analyses

The *Datasets Results* directory contains reproducible pipelines for the biological data analysed in the study. Each dataset subfolder is organised as follows:

* **Space‑Ranger data:** compressed files (e.g., `Space_Ranger_Data_Brain.zip`) provide the raw 10× Visium outputs for each sample and can be used to reproduce the preprocessing steps.
* **Methods:** scripts under `Methods/` perform quality control, normalisation, PCA, DP‑RST clustering, BAST initialisation, BayesSpace/DR.SC runs and differential gene expression analyses for the given dataset. Files are often stratified by the number of principal components.
* **Results:** subdirectories under `Results/` store intermediate and final outputs, including cluster assignments for each method, NMI scores, gene‑expression summaries and figures. The top level of *Datasets Results* also includes scripts that summarise method performance across datasets.

## Scalability experiments

The *ConsensusMCM* directory implements the consensus Monte Carlo strategy used to scale DP‑RST to high‑resolution spatial transcriptomics data. It contains:

* **Compressed input data:** `Colon.zip` and `Lung.zip` contain pre‑processed spatial coordinates and expression matrices for the 10× Visium HD colon cancer and 10× Xenium lung datasets.
* **Shard‑wise scripts:** `Colon_scripts/` and `Lung_scripts/` comprise R functions to partition the large datasets into spatial shards, run DP‑RST or BAST independently on each shard and post‑process the results. Scripts are further broken into preparation (`BAST_prep_*.R`), shard‑level inference (`run_*_shard_fun.R`) and list‑level aggregation (`run_list_*.R`).
* **Reference‑based alignment:** the `ReferenceBasedAlignment/` folder provides a Python implementation of the optimal transport alignment used to merge shard‑level clusters, an R wrapper for processing the aligned results and compressed archives of the final aligned cluster assignments.

## Figure generation

Scripts in the *Figures* folder recreate every figure in the manuscript. They include:

* **Primary figures:** R scripts for 3‑D PCA plots, hyperparameter move visualisations and Swiss‑Roll simulation summaries.
* **PCA figure assets:** a `PCs_figures/` subfolder containing PNG images of principal component densities and expression patterns across simulation scenarios. These images can be directly embedded in the supplementary materials.
* **Supporting data:** `.RData` files storing simulation data used in the figure scripts.

## Sensitivity analyses

The *Sensitivity* directory holds scripts and results used to explore the robustness of DP‑RST to the concentration parameter and the number of initial clusters. Key components include:

* `SR_sim_SensAlpha.R` and `SR_sim_SensK.R` for varying the Dirichlet Process concentration parameter (`α`) and the initial number of clusters (`K`).
* `.RData` files containing simulated data and MCMC outputs for multiple repetitions under different hyperparameter settings.
* PDF reports summarising sensitivity results (`Sensetivity_Alpha.pdf`, `Sensetivity_SpatialClusters.pdf`).

## Utility and legacy functions

* **Extra_fun/** – additional helper functions for manipulating MST cluster assignments, shredding clusters, merging MST components and plotting cluster structures in PDF format.
* **Legacy BAST/** – an archived implementation of earlier BAST routines (`BASTFun_2.R`, `ComplexDomainFun.R`, `FEMFun.R`) preserved for reference and reproducibility.
* **profiling.R** – demonstration of profiling DP‑RST with `profvis` on a Swiss‑Roll simulation.

## Dependencies

The analyses in this repository require the [DP.RST](https://github.com/rvaleriya/DP.RST) package and standard R libraries for single‑cell and spatial transcriptomics analysis. Specific dependencies include `Seurat`, `SingleCellExperiment`, `BayesSpace`, `DR.SC`, `profvis`, `igraph`, `ggplot2` and others. See the individual scripts for detailed package requirements. Before running the scripts, ensure that the DP.RST package and its dependencies are installed and available in your R environment.
