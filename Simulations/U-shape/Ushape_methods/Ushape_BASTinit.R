# Load necessary libraries
library(MASS)
library(igraph)
library(fields)
library(ggplot2)
library(colorspace)
library(ggpubr)
library(fossil)
library(gridExtra)
library(patchwork)
library(mclust)
library(fields)
library(rprojroot)

# Set working directory to the root of the R project
project_root <- find_root(is_rstudio_project)
setwd(project_root)
getwd()

# Load custom functions
source("./Codes/Custom Functions/ComplexDomainFun.R")
source("./Codes/Custom Functions/plotting_fun.R")
source("./Codes/Custom Functions/BASTFun_2.R")

# set seed
set.seed(9362)

# For debugging purposes
options(error=recover)

# Load Data
load("./Simulations/U-shape_6k/UshapeSim_coords_boundary1.RData")

# Extract first three PCAs
Y_sample = sim_data[, 4:6]
loc = sim_data[, 1:2]
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

# Create dataframe for plotting
df_subset <- data.frame(coords, Y_std)

Y1_plot <- ggplot() +
  geom_boundary(bnd_scaled) +
  geom_point(aes(x = x, y = y, col = Y1), data = df_subset) +
  labs(x = 'Coordinate X', y = 'Coordinate Y', colour = "True label") +
  ggtitle('First PCA') +
  theme(legend.position = "none")
Y1_plot

clust_plot <- ggplot() +
  geom_boundary(bnd) +
  geom_point(aes(x = x, y = y, col = factor(cluster)), data = sim_data) +
  labs(x = 'Coordinate X', y = 'Coordinate Y', colour = "True label") +
  ggtitle('Spatial Cluseters') +
  theme(legend.position = "none")
clust_plot


##### SET UP INITIAL VALUES #####
# Get mesh and triangulation
mesh = gen2dMesh(coords, bnd_scaled)
graph0 = constrainedDentri(n, mesh)
E(graph0)$eid = c(1:ecount(graph0))  # edge id
V(graph0)$vid = c(1:vcount(graph0))  # vertex id
mstgraph = mst(graph0)  # initial spanning tree
graph0 = delete_edge_attr(graph0, 'weight')
mstgraph0 = delete_edge_attr(mstgraph, 'weight')

# plot spatial graph
spatial_graph <- plotGraph(coords, graph0) +
  geom_boundary(bnd_scaled) +
  labs(x = 'Coordinate X', y = 'Coordinate Y') +
  ggtitle('Spatial Graph')
spatial_graph

# Define the "temperatures"
temp = c(1) # no parallel tempering
M = length(temp)

# Initialize from one cluster only
cluster = matrix(rep(1), nrow = n, ncol = M) # initial cluster memberships
cluster_ids <- unique(cluster[,1])
k_max = 20 # maximum number of clusters

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
init_val[['cluster']] = cluster 
init_val[['sigmasq_y']] = sigmasq_y

hyperpar = list()
hyperpar[['sigmasq_mu']] = (0.5/(2*sqrt(1)))^2 
hyperpar[['lambda_s']] = diag(1, p) 
hyperpar[['nu']] = n-p+1
hyperpar[['lambda_k']] = 20
hyperpar[['M']] = M
hyperpar[['k_max']] = k_max


# MCMC parameters
# number of posterior samples = (MCMC - BURNIN) / THIN
MCMC = 7000 # MCMC iterations
BURNIN = 5000 # burnin period length
THIN = 5       # thinning intervals

dynamic_part <- "U_sim_rep" # For backup files

# Generate seeds to run algorithm multiple times
seeds <- sample.int(1000, size = 30, replace = FALSE)
seeds

# Create a list to save results from all the runs
U_sim_PT_reps <- list()

# Run the BAST algorithm to obtain initialization for DP-RST
for (s in 1:length(seeds)){
  
  U_sim_PT = fitBAST(Y_std, graph0, init_val, hyperpar, temp, 
                      MCMC, BURNIN, THIN, 
                      PT = FALSE, seed = seeds[s],
                      backup_d = sprintf("./Simulations/U-shape_6k/%s_backup.RData", dynamic_part))
  U_sim_PT$seed <- seeds[s]
  
  U_sim_PT_reps[[s]] <- U_sim_PT
}

# Save the results of the repetitive runs
save(U_sim_PT_reps, file = "./Simulations/U-shape_6k/U_sim_PT_30reps.RData")
