writeLines("R script started", stderr())

start.time_file <- Sys.time()

# Set a writable library path
my_lib_path <- "/scratch/user/varogovchenko/Rlibs"
.libPaths(my_lib_path)

library(rlang)
library(MASS)
library(igraph)
library(fields)
library(mclust)
library(dplyr)

print("Libraries are loaded")

setwd("/scratch/user/varogovchenko/BASTION_HPRC")
out_dir <- "/scratch/user/varogovchenko/BASTION_HPRC/GSE169749_RAW-2/NewRuns_2025/Sensitivity_K_PC10/Outputs_BAST"

num_clusters <- 30 # c(5, 10, 15, 25, 30)
dynamic_part <- sprintf("Gut_BAST_PT_k%02d", num_clusters)
backup_path  <- file.path(out_dir, paste0(dynamic_part, "_backup.RData"))
out_path     <- file.path(out_dir, paste0(dynamic_part, "_OutputOnly.RData"))

source("BASTFun_2.R")
source("ComplexDomainFun.R")

# set seed
set.seed(9362)

# # For debugging purposes
# options(error=recover)

gut_df_wt_muscle = readRDS("GSE169749_RAW-2/NewRuns_2025/gut_df_wt_muscle.rds")
load("/scratch/user/varogovchenko/BASTION_HPRC/GSE169749_RAW-2/swiss_roll_wt_muscle_boundary.RData")

# Extract first three PCAs
Y_sample = gut_df_wt_muscle[, 1:10]
print(head(Y_sample))
loc = gut_df_wt_muscle[, 16:17]
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

# Get mesh and triangulation
# mesh = gen2dMesh(coords, bnd_scaled_r)
load("/scratch/user/varogovchenko/BASTION_HPRC/GSE169749_RAW-2/swiss_roll_wt_muscle_mesh.RData")
graph0 = constrainedDentri(n, mesh)
E(graph0)$eid = c(1:ecount(graph0))  # edge id
V(graph0)$vid = c(1:vcount(graph0))  # vertex id
mstgraph = mst(graph0)  # initial spanning tree
graph0 = delete_edge_attr(graph0, 'weight')
mstgraph0 = delete_edge_attr(mstgraph, 'weight')

# Define the "temperatures"
temp = seq(0.1:1, by = 0.1)
# Now M is not the number of weak learners but the number of temperatures
M = length(temp)

# Compute the geodesic distance matrix
geodesic_dist_matrix <- distances(graph0, v = V(graph0), to = V(graph0),
                                  weights = E(graph0)$weight)

geodesic_dist_matrix[is.infinite(geodesic_dist_matrix)] <- max(geodesic_dist_matrix[is.finite(geodesic_dist_matrix)], na.rm = TRUE) * 2

# Perform hierarchical clustering on the distance matrix
hc <- fastcluster::hclust(as.dist(geodesic_dist_matrix), method = "ward.D2")

# Cut the tree to form clusters
print("Number of spatial clusters")
print(num_clusters)

clusters <- cutree(hc, k = num_clusters)

# Compute centroids for each cluster
centroids <- aggregate(Y_std, by = list(clusters), FUN = mean)

# Store the centroids and cluster assignments 
centroids <- centroids[, -1]  # Remove cluster labels

# Initialize from one cluster only
cluster = matrix(rep(clusters), nrow = n, ncol = M) # initial cluster memberships
cluster_ids <- unique(clusters)
k_max = length(cluster_ids) # maximum number of clusters

mu = list() # list for initial values of mu
sigmasq_y = list()
mstgraph_lst = list()  # initial spanning trees

# For each temperature, we have same initial cluster assignment
for(m in 1:M) {
  mu[[m]] = centroids
  sigmasq_y[[m]] = cov(Y_std)
  mstgraph_lst[[m]] = mstgraph0
}

# find the scale matrix for the Inverse-Wishart prior distribution of the data cbraince matrix
nu = p

init_val = list()
init_val[['trees']] = mstgraph_lst
init_val[['mu']] = mu
init_val[['cluster']] = cluster 
init_val[['sigmasq_y']] = sigmasq_y

hyperpar = list()
hyperpar[['sigmasq_mu']] = (0.5/(2*sqrt(1)))^2 
hyperpar[['lambda_s']] = diag(1, p) # lambda_s
hyperpar[['nu']] = nu
hyperpar[['lambda_k']] = num_clusters
hyperpar[['M']] = M
hyperpar[['k_max']] = k_max


##### TRAIN DPM BAST WITH SEED 7394 #####

# MCMC parameters
# number of posterior samples = (MCMC - BURNIN) / THIN
MCMC = 50000 # MCMC iterations
BURNIN = 40000 # burnin period length
THIN = 5      # thinning intervals

print("All the parameters are set")

start.time_mcmc <- Sys.time()
Gut_BAST_PT = fitBAST(Y_std, graph0, init_val, hyperpar, temp,
                      MCMC, BURNIN, THIN,
                      PT = TRUE, seed = 7394,
                      backup_d = backup_path)
end.time_mcmc <- Sys.time()

mcmc_time <- end.time_mcmc - start.time_mcmc
mcmc_time

save(Gut_BAST_PT, file = out_path)

end.time_file <- Sys.time()

file_time <- end.time_file - start.time_file
file_time