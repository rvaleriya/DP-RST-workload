# Load necessary libraries
library(ggplot2)
library(mclust)
library(paletteer)
library(dplyr)
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

# Create a list to store results from different runs
U_sim_kmeans_reps <- list()

for (s in 1:length(seeds)){
  
  set.seed(seeds[s])
  
  kmeans_fit <- kmeans(Y_std, 6)$cluster
  U_sim_kmeans_reps[[s]] <- kmeans_fit
}

# Save the results of the repetitive runs
save(U_sim_kmeans_reps, file = "./Simulations/U-shape_6k/U_sim_kmeans_30reps.RData")

