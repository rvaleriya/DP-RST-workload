# Load necessary libraries
library(SingleCellExperiment)
library(BayesSpace)
library(ggplot2)
library(mclust)
library(paletteer)
library(rprojroot)

# Set working directory to the root of the R project
project_root <- find_root(is_rstudio_project)
setwd(project_root)
getwd()

# Set seed
set.seed(9362)

# Load custom functions
source("./Codes/Custom Functions/ComplexDomainFun.R")
source("./Codes/Custom Functions/plotting_fun.R")

### Load Data ###
load("./Simulations/U-shape_6k/UshapeSim_coords_boundary1.RData")
load("./Simulations/U-shape_6k/U_sim_PT_30reps.RData")

# Get seeds from the BAST run for initialization
seeds <- c()
for (i in 1:length(U_sim_PT_reps)){
  seeds[i] <- U_sim_PT_reps[[i]][["seed"]]
}
rm(U_sim_PT_reps) # Remove BASTinit data

# Extract first three PCAs
Y_sample = sim_data[, 4:6]
loc = sim_data[, 1:2]

n = nrow(Y_sample) # number of observations
p = ncol(Y_sample) # number of variables in Y

# Standardize coordinates, Y values
coords <- apply(loc, 2, scale)
Y_std <- apply(Y_sample, 2, scale)

rownames(Y_std) <- paste0("Cell", 1:n)
colnames(Y_std) <- paste0("PC", 1:p)
rownames(coords) <- rownames(Y_std)

# Create the SingleCellExperiment object
sce <- SingleCellExperiment(
  reducedDims = list(PCA = Y_std),  # Add PCA data as an assay
  colData = coords        # Add spatial coordinates as column data
)

# Rename the coordinates columns to match BayesSpace expectations
colnames(colData(sce)) <- c("col", "row")

# Create a list to store results from different runs
U_sim_BayesSpace_reps <- list()

for (s in 1:length(seeds)){
  set.seed(seeds[s])
  
  BS_results <- spatialCluster(sce, use.dimred = "PCA", d = 3, q = 6, platform = "ST")
  
  U_sim_BayesSpace_reps[[s]] <- BS_results
}

# Save the results of the repetitive runs
save(U_sim_BayesSpace_reps, file = "./Simulations/U-shape_6k/U_sim_BayesSpace_30reps.RData")
