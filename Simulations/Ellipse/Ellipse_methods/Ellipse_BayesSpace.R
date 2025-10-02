# Load necessary libraries
library(SingleCellExperiment)
library(BayesSpace)

#-------------------------------------------------------------------------------

# Set working directory to the root of the R project
setwd("~/Desktop/DP-RST-workload")

# Set seed
set.seed(9362)

### Load Data ###
Ellipse_sim_data <- readr::read_csv("./New_Simulations/Ellipse/Ellipse_sim_data.csv")

#-------------------------------------------------------------------------------

seeds <- c(141, 549, 75, 492, 676, 179, 587, 592, 601, 916,
           518, 339, 921, 423, 330, 388, 273, 286, 61, 807,
           283, 127, 165, 952, 311, 597, 473, 594, 605, 101)

#-------------------------------------------------------------------------------

# Extract first three PCAs
Y_sample = Ellipse_sim_data[, 5:7] # 5:7 for PCs 1-3 and 5:14 for PCs 1-10
head(Y_sample)

loc = Ellipse_sim_data[, 1:2]
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

# Rename the coordinates columns to match BayesSpace expectations
colnames(colData(sce)) <- c("col", "row")

# Create a list to store results from different runs
Ellipse_sim_BayesSpace_p3_reps <- list()

for (s in 1:length(seeds)){
  set.seed(seeds[s])
  
  BS_results <- spatialCluster(sce, use.dimred = "PCA", d = p, q = 5, platform = "ST")
  
  Ellipse_sim_BayesSpace_p3_reps[[s]] <- BS_results
}

# Save the results of the repetitive runs
save(Ellipse_sim_BayesSpace_p3_reps, file = "./New_Simulations/Ellipse/Ellipse_results/Ellipse_sim_BayesSpace_p3_30reps.RData")
