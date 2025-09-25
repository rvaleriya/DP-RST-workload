writeLines("R script started", stderr())

start.time_file <- Sys.time()

options(repos = c(CRAN = "https://cran.rstudio.com"))

Sys.setlocale("LC_ALL", "en_US.UTF-8")

# Set a writable library path
my_lib_path <- ""
# if (!dir.exists(my_lib_path)) {
#   dir.create(my_lib_path)
# }
.libPaths(my_lib_path)

library(rlang)
library(MASS)
library(igraph)
library(fields)
library(mclust)
library(dplyr)
library(fdaPDE)

print("Libraries are loaded")

setwd("/BASTION_HPRC")

source("BASTFun_2.R")
source("ComplexDomainFun.R")

# set seed
set.seed(9362)

# # For debugging purposes
# options(error=recover)

load("/Brain/DFPLC_151510_data_for_iIMPACT.RData")
load("/Brain/brain_boundary.RData")

brain_BS_results <- readRDS("/Brain/NewRuns_2025/BayesSpace_Brain_clustering_results.rds")
colnames(brain_BS_results)[colnames(brain_BS_results) == "imagerow"] <- "x"
colnames(brain_BS_results)[colnames(brain_BS_results) == "imagecol"] <- "y"
head(brain_BS_results)

# Extract first three PCAs
Y_sample = Y[, 1:3]
bnd = boundary

n = nrow(Y_sample) # number of observations
p = ncol(Y_sample) # number of variables in Y

# Standardize coordinates, Y values
coords <- apply(loc, 2, scale)
Y_std <- apply(Y_sample, 2, scale)

#Standartize the boundary given the params of coordinates
bnd_scaled = list()
bnd_scaled$x <- (bnd$x - mean(loc[,1]))/sd(loc[,1])
bnd_scaled$y <- (bnd$y - mean(loc[,2]))/sd(loc[,2])

# Add small random variable to the scaled boundary to avoid duplicates
bnd_scaled_r = list()
bnd_scaled_r$x <- bnd_scaled$x + c(0, runif(length(bnd_scaled$x)-2, min = 1e-04, max = 5e-04), 0)
bnd_scaled_r$y <- bnd_scaled$y + c(0, runif(length(bnd_scaled$y)-2, min = 1e-04, max = 5e-04),0)

# Create dataframe for plotting
df_subset <- data.frame(coords, Y_std)


##### SET UP INITIAL VALUES #####

coords_init <- loc %>%
  left_join(brain_BS_results, by = c("x", "y"))

coords_init[, c("x", "y")] <- coords

# Check the result
head(coords_init)


compute_MST_spatial <- function(data,
                                coords_cols = c("x", "y"),
                                cluster_col = "cluster",
                                bnd = NULL,
                                threshold = 5000,
                                penalty = NULL) {

  # Convert coordinates to a matrix
  coords <- as.matrix(data[, coords_cols])

  # If no boundary is supplied, compute the convex hull and form a boundary
  # if (is.null(bnd)) {
  #   hull_idx <- chull(coords)
  #   # Ensure boundary is closed by repeating the first point at the end
  #   bnd <- list(x = coords[c(hull_idx, hull_idx[1]), 1],
  #               y = coords[c(hull_idx, hull_idx[1]), 2])
  # }

  # Generate a mesh (requires your gen2dMesh function)
  # mesh <- gen2dMesh(coords, bnd)
  load("/Brain/brain_mesh.RData")

  # Create the graph using constrainedDentri (assumes function available)
  graph0 <- constrainedDentri(n = nrow(coords), mesh = mesh, threshold = threshold)
  E(graph0)$eid = c(1:ecount(graph0))  # edge id
  V(graph0)$vid = c(1:vcount(graph0))  # vertex id

  # Assign the supplied clustering as a vertex attribute
  V(graph0)$cluster <- data[[cluster_col]]

  # Get the original edge weights (computed in constrainedDentri)
  orig_weights <- E(graph0)$weight

  # If penalty is not provided, set it to max weight + 1
  if (is.null(penalty)) {
    penalty <- max(orig_weights, na.rm = TRUE) + 1
  }

  # Extract edge endpoints and their cluster memberships
  edge_ends <- as.data.frame(get.edgelist(graph0))
  colnames(edge_ends) <- c("V1", "V2")
  edge_ends$V1 <- as.numeric(edge_ends$V1)
  edge_ends$V2 <- as.numeric(edge_ends$V2)

  clusters_v1 <- V(graph0)$cluster[edge_ends$V1]
  clusters_v2 <- V(graph0)$cluster[edge_ends$V2]

  # Identify edges connecting different clusters
  inter_cluster <- clusters_v1 != clusters_v2

  # Adjust edge weights: add penalty for inter-cluster edges
  new_weights <- orig_weights + penalty * inter_cluster
  E(graph0)$new_weight <- new_weights

  # Compute the MST using the penalized weights
  mst_graph <- mst(graph0, weights = E(graph0)$new_weight)

  # Identify spatial components:
  # Remove MST edges that connect points with different initial clusters.
  head_idx <- as.numeric(head_of(mst_graph, E(mst_graph)))
  tail_idx <- as.numeric(tail_of(mst_graph, E(mst_graph)))
  diff_clusters <- V(mst_graph)$cluster[head_idx] != V(mst_graph)$cluster[tail_idx]

  mst_same_label <- delete_edges(mst_graph, which(diff_clusters))

  # Compute connected components (each component is a "spatial cluster")
  comp <- igraph::components(mst_same_label)

  # Append the spatial cluster membership to the data frame
  data$spatial_cluster <- comp$membership

  # Return a list containing the MST and the updated data frame
  return(list(mst = mst_graph, data = data))
}


init_partition <- compute_MST_spatial(coords_init,
                                      coords_cols = c("x", "y"),
                                      cluster_col = "cluster_q7_pc3",
                                      bnd = bnd_scaled)

# groups_plot(init_partition$data[, c("x", "y")], init_partition$data$spatial_cluster, bnd_scaled)

# Get mesh and triangulation
# mesh = gen2dMesh(coords, bnd_scaled_r)
load("/Brain/brain_mesh.RData")
graph0 = constrainedDentri(n, mesh)
E(graph0)$eid = c(1:ecount(graph0))  # edge id
V(graph0)$vid = c(1:vcount(graph0))  # vertex id
graph0 = delete_edge_attr(graph0, 'weight')

mstgraph0 = delete_edge_attr(init_partition$mst, 'weight')

# Define the "temperatures"
temp = seq(0.1:1, by = 0.1)
# Now M is not the number of weak learners but the number of temperatures
M = length(temp)

mu = list() # list for initial values of mu
sigmasq_y = list()
mstgraph_lst = list()  # initial spanning trees

clusters_bast = init_partition$data$spatial_cluster
k_max = length(unique(init_partition$data$spatial_cluster))
print(k_max)

cluster = matrix(rep(clusters_bast), nrow = n, ncol = M, byrow = F)
cluster_means_matrix <- matrix(0, nrow = k_max, ncol = p, byrow = F)

# Get the means of each cluster
for (i in 1:k_max) {
  cluster_means_matrix[i,] <- colMeans(Y_std[c(which(clusters_bast == i)), , drop = FALSE])
}

# For each temperature, we have same initial cluster assignment
for(m in 1:M) {
  mu[[m]] = cluster_means_matrix
  sigmasq_y[[m]] = cov(Y_std) #Covariance matrix
  mstgraph_lst[[m]] = mstgraph0
}


# find the scale matrix for the Inverse-Wishart prior distribution of the data cbraince matrix
nu = p

init_val = list()
init_val[['trees']] = mstgraph_lst
init_val[['mu']] = mu
init_val[['cluster']] = cluster #matrix(close_truth, nrow = length(close_truth), ncol = M, byrow = FALSE) #matrix(rep(kmeans(Y_std, 4)$cluster), nrow = n, ncol = M, byrow = F) #cluster #as.matrix(true_clusters) # as.matrix(kmeans(Y_std, 4)$cluster) # cluster
init_val[['sigmasq_y']] = sigmasq_y

hyperpar = list()
hyperpar[['sigmasq_mu']] = (0.5/(2*sqrt(1)))^2 #(0.5/(2*sqrt(M)))^2 # leave it to be as in univariate case (just in the code write as a matrix)
hyperpar[['lambda_s']] = diag(1, p) # lambda_s
hyperpar[['nu']] = nu
hyperpar[['lambda_k']] = 20
hyperpar[['M']] = M
hyperpar[['k_max']] = k_max


##### TRAIN DPM BAST WITH SEED 7394 #####

# MCMC parameters
# number of posterior samples = (MCMC - BURNIN) / THIN
MCMC = 30000 # MCMC iterations
BURNIN = 25000 # burnin period length
THIN = 5       # thinning intervals

dynamic_part <- "Brain_BAST_PT_Version1_newBS"

print("All the parameters are set")

start.time_mcmc <- Sys.time()
Brain_BAST_PT_Version1_newBS = fitBAST(Y_std, graph0, init_val, hyperpar, temp,
                                MCMC, BURNIN, THIN,
                                PT = TRUE, seed = 7394,
                                backup_d = sprintf("/Brain/NewRuns_2025/Outputs_BAST/%s_backup.RData", dynamic_part))
end.time_mcmc <- Sys.time()

mcmc_time <- end.time_mcmc - start.time_mcmc
mcmc_time

# save(Brain_BAST_PT_Version1, file = "/Brain/NewRuns_2025/Outputs_BAST/Brain_BAST_PT_Version1_OutputOnly.RData", compress = TRUE)
saveRDS(Brain_BAST_PT_Version1_newBS, file = "/Brain/NewRuns_2025/Outputs_BAST/Brain_BAST_PT_Version1_newBS_OutputOnly.rds")

end.time_file <- Sys.time()

file_time <- end.time_file - start.time_file
file_time
