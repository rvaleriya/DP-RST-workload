# Load necessary libraries
library(SingleCellExperiment)
library("SC.MEB")
library(ggplot2)
library(mclust)
library(paletteer)
library(rprojroot)

# Set seed
set.seed(9362)

# Set working directory to the root of the R project
project_root <- find_root(is_rstudio_project)
setwd(project_root)
getwd()

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

# Rename the coordinates columns
colnames(colData(sce)) <- c("col", "row")

### Set the parameters ###
beta_grid = seq(0, 5, 0.2) # Vector specifying the smoothness of Random Markov Field
K_set= 2:10 # Numbers of mixture components
parallel=TRUE #Logical value specifing the run the model in parallel
num_core = 3
PX = TRUE # Logical value for paramter expansion in EM algorithm
maxIter_ICM = 10 # Maximum iteration of ICM algorithm
maxIter = 50 # Maximum iteration of EM algorithm

#### Calculating the neighborhood from coordinates
Adj_sp  <- find_neighbors2(sce, platform = "ST")

# Create a list to store results from different runs
U_sim_SC.MEB_reps <- list()

for (s in 1:length(seeds)){
  
  set.seed(seeds[s])
  
  #### Run the SC-MEB in parallel
  SC.MEB_fit = SC.MEB(y = Y_std, Adj_sp, beta_grid = beta_grid, K_set = K_set, 
                      parallel = parallel, num_core = num_core, PX = PX, 
                      maxIter_ICM = maxIter_ICM, maxIter = maxIter)
  
  U_sim_SC.MEB_reps[[s]] <- SC.MEB_fit
}

# Save the results of the repetitive runs
save(U_sim_SC.MEB_reps, file = "./Simulations/U-shape_6k/U_sim_SC.MEB_30reps.RData")
