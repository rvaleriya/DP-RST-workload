# Load necessary libraries

library(SingleCellExperiment)
library("SC.MEB")

#-------------------------------------------------------------------------------

# Set working directory to the root of the R project
setwd("~/Desktop/DP-RST-workload")

# Set seed
set.seed(9362)

### Load Data ###
Curl_sim_data <- readr::read_csv("./Simulations/Curl/Curl_sim_data.csv")

#-------------------------------------------------------------------------------

seeds <- c(141, 549, 75, 492, 676, 179, 587, 592, 601, 916,
           518, 339, 921, 423, 330, 388, 273, 286, 61, 807,
           283, 127, 165, 952, 311, 597, 473, 594, 605, 101)

#-------------------------------------------------------------------------------

# Extract first three PCAs
Y_sample = Curl_sim_data[, 5:7] # 5:7 for PCs 1-3 and 5:14 for PCs 1-10
head(Y_sample)

loc = Curl_sim_data[, 1:2]
head(loc)

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

#-------------------------------------------------------------------------------

### Set the parameters ###
beta_grid = seq(0, 5, 0.2) # Vector specifying the smoothness of Random Markov Field
K_set = 10 # Numbers of mixture components
parallel=TRUE #Logical value specifing the run the model in parallel
num_core = 3
PX = TRUE # Logical value for paramter expansion in EM algorithm
maxIter_ICM = 10 # Maximum iteration of ICM algorithm
maxIter = 50 # Maximum iteration of EM algorithm

#### Calculating the neighborhood from coordinates
Adj_sp  <- find_neighbors2(sce, platform = "ST")

# Create a list to store results from different runs
Curl_sim_SC.MEB_p3_reps <- list()

for (s in 1:length(seeds)){
  
  set.seed(seeds[s])
  
  #### Run the SC-MEB in parallel
  SC.MEB_fit = SC.MEB(y = Y_std, Adj_sp, beta_grid = beta_grid, K_set = K_set, 
                      parallel = parallel, num_core = num_core, PX = PX, 
                      maxIter_ICM = maxIter_ICM, maxIter = maxIter)
  
  Curl_sim_SC.MEB_p3_reps[[s]] <- SC.MEB_fit
}

# Save the results of the repetitive runs
save(Curl_sim_SC.MEB_p3_reps, file = "./Simulations/Curl/Curl_results/Curl_sim_SC.MEB_p3_30reps.RData")
