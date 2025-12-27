# Load necessary libraries
library(SingleCellExperiment)
library(BayesSpace)
library(readr)

#-------------------------------------------------------------------------------

# Set working directory
setwd("~/Desktop/DP-RST-workload")

# Load Data
Curl_sim_data <- read_csv("./Simulations/Curl/Curl_sim_data.csv")

# Define Seeds
seeds <- c(141, 549, 75, 492, 676, 179, 587, 592, 601, 916,
           518, 339, 921, 423, 330, 388, 273, 286, 61, 807,
           283, 127, 165, 952, 311, 597, 473, 594, 605, 101)

#-------------------------------------------------------------------------------
# CONFIGURATION LOOP
# Define the column ranges for the two scenarios
# p3 uses columns 5:7, p10 uses columns 5:14
configs <- list(
  p3  = 5:7,
  p10 = 5:14
)

# Loop through both configurations
for (config_name in names(configs)) {
  
  cat(paste0("\nProcessing configuration: ", config_name, "...\n"))
  
  # 1. Prepare Data for this specific run
  col_range <- configs[[config_name]]
  Y_sample <- Curl_sim_data[, col_range]
  
  loc <- Curl_sim_data[, 1:2]
  
  n <- nrow(Y_sample)
  p <- ncol(Y_sample) # Automatically detects 3 or 10
  
  # Standardize
  coords <- apply(loc, 2, scale)
  Y_std <- apply(Y_sample, 2, scale)
  
  rownames(Y_std) <- paste0("Cell", 1:n)
  colnames(Y_std) <- paste0("PC", 1:p)
  rownames(coords) <- rownames(Y_std)
  
  # Create SCE Object
  sce <- SingleCellExperiment(
    reducedDims = list(PCA = Y_std),
    colData = coords
  )
  colnames(colData(sce)) <- c("col", "row")
  
  # 2. Run the 30 Reps
  # Initialize storage list
  results_list <- list()
  
  for (s in 1:length(seeds)) {
    set.seed(seeds[s])
    
    # Run BayesSpace
    # use.dimred = "PCA", d = p (3 or 10), q = 10 clusters
    BS_results <- spatialCluster(sce, use.dimred = "PCA", d = p, q = 10, platform = "ST")
    
    results_list[[s]] <- BS_results
    
    if(s %% 5 == 0) cat(paste("  Finished seed", s, "of 30\n"))
  }
  
  # 3. Dynamic Saving
  # Construct the specific object name (e.g., Curl_sim_BayesSpace_p3_reps)
  object_name <- paste0("Curl_sim_BayesSpace_", config_name, "_reps")
  
  # Assign the results list to that specific name
  assign(object_name, results_list)
  
  # Construct filename
  file_name <- paste0("./Simulations/Curl/Curl_results/", object_name, "_30reps.RData")
  
  # Save using the specific variable name
  save(list = object_name, file = file_name)
  
  cat(paste("Saved:", file_name, "\n"))
}

cat("\nAll runs completed successfully.\n")