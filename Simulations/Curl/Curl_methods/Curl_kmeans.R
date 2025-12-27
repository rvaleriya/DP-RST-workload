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

# Create a list to store results from different runs
Curl_sim_kmeans_3p_reps <- list()

for (s in 1:length(seeds)){
  
  set.seed(seeds[s])
  
  kmeans_fit <- kmeans(Y_std, 10)$cluster
  Curl_sim_kmeans_3p_reps[[s]] <- kmeans_fit
}

# Save the results of the repetitive runs
save(Curl_sim_kmeans_3p_reps, file = "./Simulations/Curl/Curl_results/Curl_sim_kmeans_3p_30reps.RData")

#-------------------------------------------------------------------------------
# Extract first three PCAs
Y_sample = Curl_sim_data[, 5:14] # 5:7 for PCs 1-3 and 5:14 for PCs 1-10
head(Y_sample)

loc = Curl_sim_data[, 1:2]
head(loc)

n = nrow(Y_sample) # number of observations
p = ncol(Y_sample) # number of variables in Y

# Standardize coordinates, Y values
coords <- apply(loc, 2, scale)
Y_std <- apply(Y_sample, 2, scale)

# Create a list to store results from different runs
Curl_sim_kmeans_10p_reps <- list()

for (s in 1:length(seeds)){
  
  set.seed(seeds[s])
  
  kmeans_fit <- kmeans(Y_std, 10)$cluster
  Curl_sim_kmeans_10p_reps[[s]] <- kmeans_fit
}

# Save the results of the repetitive runs
save(Curl_sim_kmeans_10p_reps, file = "./Simulations/Curl/Curl_results/Curl_sim_kmeans_10p_30reps.RData")
