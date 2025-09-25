writeLines("R script started", stderr())

start.time_file <- Sys.time()

# Set a writable library path
my_lib_path <- ""
# if (!dir.exists(my_lib_path)) {
#   dir.create(my_lib_path)
# }
.libPaths(my_lib_path)

library(DP.RST)
library(MASS)
library(igraph)
library(fields)
# library(fossil)
library(mclust)
library(fields)

print("Libraries are loaded")

#-------------------------------------------------------------------------------
setwd("")

source("BASTFun_2.R")
source("ComplexDomainFun.R")

Ellipse_sim_data <- readr::read_csv("./Simulations/Ellipse/Ellipse_sim_data.csv")
head(Ellipse_sim_data)

# set seed
set.seed(9362)

# # For debugging purposes
# options(error=recover)

#-------------------------------------------------------------------------------

# Extract first three PCAs
Y_sample = Ellipse_sim_data[, 5:7]
loc = Ellipse_sim_data[, 1:2]

n = nrow(Y_sample) # number of observations
p = ncol(Y_sample) # number of variables in Y

# Standardize coordinates, Y values
coords <- apply(loc, 2, scale)
Y_std <- apply(Y_sample, 2, scale)

#-------------------------------------------------------------------------------

mesh <- fdaPDE::create.mesh.2D(coords)

# Apply trimming function to cut long edges
trim_graph <- function(n, mesh, threshold = 100) {
  coords <- mesh$nodes[1:n, ]
  edge_list <- mesh$edges  # No filtering by bnd_edges, since you have no boundary edges
  
  # Compute distances
  distance <- sqrt(rowSums((coords[edge_list[, 1], ] - coords[edge_list[, 2], ])^2))
  
  # Drop long edges
  rid_drop <- distance > threshold
  edge_list <- edge_list[!rid_drop, ]
  distance <- distance[!rid_drop]
  
  # Build graph
  graph0 <- graph_from_edgelist(edge_list, directed = FALSE)
  E(graph0)$weight <- distance
  return(graph0)
}

graph0 <- trim_graph(n = nrow(coords), mesh = mesh, threshold = 0.6)

# plot(graph0, layout = as.matrix(coords), vertex.size = 1, vertex.label = NA, edge.width = 0.5)


##### SET UP INITIAL VALUES #####
# Get mesh and triangulation
# mesh = gen2dMesh(coords, bnd_scaled)
# graph0 = constrainedDentri(n, mesh)
E(graph0)$eid = c(1:ecount(graph0))  # edge id
V(graph0)$vid = c(1:vcount(graph0))  # vertex id
mstgraph = mst(graph0)  # initial spanning tree
graph0 = delete_edge_attr(graph0, 'weight')
mstgraph0 = delete_edge_attr(mstgraph, 'weight')

# # plot spatial graph
# spatial_graph <- plotGraph(coords, graph0) +
#   # geom_boundary(bnd_scaled) +
#   scale_y_reverse() +
#   labs(x = 'Coordinate X', y = 'Coordinate Y') +
#   ggtitle('Spatial Graph')
# spatial_graph

#-------------------------------------------------------------------------------

# Define the "temperatures"
temp = c(1) #seq(0.1:1, by = 0.1)
# Now M is not the number of weak learners but the number of temperatures
M = length(temp)

#-------------------------------------------------------------------------------

# Compute the geodesic distance matrix for the selected batch of nodes
geodesic_dist_matrix <- distances(graph0, v = V(graph0), to = V(graph0),
                                  weights = E(graph0)$weight)

# Replace Inf with a large value
geodesic_dist_matrix[is.infinite(geodesic_dist_matrix)] <- max(geodesic_dist_matrix[is.finite(geodesic_dist_matrix)], na.rm = TRUE) * 2

# Perform hierarchical clustering on the combined distance matrix
hc <- fastcluster::hclust(as.dist(geodesic_dist_matrix), method = "ward.D2")

clusters <- cutree(hc, k = 15)

# init_clust_plot <- ggplot() +
#   # geom_boundary(bnd) +
#   scale_y_reverse() +
#   geom_point(aes(x = X, y = Y, col = factor(clusters)), data = Ellipse_sim_data) +
#   labs(x = 'Coordinate X', y = 'Coordinate Y', colour = "True label") +
#   ggtitle('Initial Spatial Cluseters') +
#   theme(legend.position = "none")
# init_clust_plot

#-------------------------------------------------------------------------------

# Initialize from one cluster only
cluster = matrix(clusters, nrow = n, ncol = M) # initial cluster memberships
cluster_ids <- unique(cluster)
k_max = length(cluster_ids) # maximum number of clusters

mu = list() # list for initial values of mu
sigmasq_y = list()
mstgraph_lst = list()  # initial spanning trees

cluster_means_matrix <- matrix(0, nrow = length(cluster_ids), ncol = p)
# Get the means of each cluster
for (i in 1:length(cluster_ids)) {
  cluster_means_matrix[i, ] <- colMeans(Y_std[cluster[,1] == cluster_ids[i], , drop = FALSE])
}

# For each temperature, we have same initial cluster assignment
for(m in 1:M) {
  mu[[m]] = cluster_means_matrix
  sigmasq_y[[m]] = cov(Y_std) #Covariance matrix
  mstgraph_lst[[m]] = mstgraph0
}


init_val = list()
init_val[['trees']] = mstgraph_lst
init_val[['mu']] = mu
init_val[['cluster']] = cluster #matrix(close_truth, nrow = length(close_truth), ncol = M, byrow = FALSE) #matrix(rep(kmeans(Y_std, 4)$cluster), nrow = n, ncol = M, byrow = F) #cluster #as.matrix(true_clusters) # as.matrix(kmeans(Y_std, 4)$cluster) # cluster
init_val[['sigmasq_y']] = sigmasq_y

hyperpar = list()
hyperpar[['sigmasq_mu']] = (0.5/(2*sqrt(1)))^2 #(0.5/(2*sqrt(M)))^2 # leave it to be as in univariate case (just in the code write as a matrix)
hyperpar[['lambda_s']] = diag(1, p) # lambda_s
hyperpar[['nu']] = p
hyperpar[['lambda_k']] = 10
hyperpar[['M']] = M
hyperpar[['k_max']] = k_max


##### TRAIN DPM BAST WITH SEED 7394 #####

# MCMC parameters
# number of posterior samples = (MCMC - BURNIN) / THIN
MCMC = 15000 # MCMC iterations
BURNIN = 2000 # burnin period length
THIN = 5       # thinning intervals

dynamic_part <- "Ellipse_sim_rep"

seeds <- sample.int(1000, size = 30, replace = FALSE)
seeds

Ellipse_sim_BAST_p3_reps <- list()

for (s in 1:length(seeds)){

  Ellipse_sim_BAST = fitBAST(Y_std, graph0, init_val, hyperpar, temp,
                      MCMC, BURNIN, THIN,
                      PT = FALSE, seed = seeds[s],
                      backup_d = sprintf("./Simulations/Ellipse/%s_backup.RData", dynamic_part))

  Ellipse_sim_BAST$seed <- seeds[s]

  Ellipse_sim_BAST_p3_reps[[s]] <- Ellipse_sim_BAST

  cat('Loop', s, 'done\n')
}

# Save the results of the repetitive runs
save(Ellipse_sim_BAST_p3_reps, file = "Simulations/Ellipse/Ellipse_sim_BAST_p3_30reps.RData")

end.time_file <- Sys.time()

file_time <- end.time_file - start.time_file
file_time
