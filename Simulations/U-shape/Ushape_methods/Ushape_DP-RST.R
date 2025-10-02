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
library(paletteer)
library(rprojroot)

# Set working directory to the root of the R project
project_root <- find_root(is_rstudio_project)
setwd(project_root)
getwd()

# Load custom functions
source("./Codes/Custom Functions/ComplexDomainFun.R")
source("./Codes/Custom Functions/plotting_fun.R")
source("./Codes/Custom Functions/DP.RST_fun.R")

# set seed
set.seed(9362)

# For debugging purposes
options(error=recover)

### Load Data ###
load("./Simulations/U-shape_6k/UshapeSim_coords_boundary1.RData")
load("./Simulations/U-shape_6k/U_sim_PT_30reps.RData")

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
  ggtitle(paste("First Y value, n = ", n, ", p = ", p, sep = ""))
Y1_plot


my_palette <- paletteer_dynamic("cartography::multi.pal", 20)

clust_plot <- ggplot() +
  geom_boundary(bnd) +
  geom_point(aes(x = x, y = y, col = factor(cluster)), data = sim_data) +
  scale_color_manual(values = my_palette) +
  labs(x = 'Coordinate X', y = 'Coordinate Y', colour = "True label") +
  ggtitle(paste("True Cluseters, k = ", length(unique(sim_data$cluster)), sep = "")) +
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

# Define the "temperatures"
temp = seq(0.1:1, by = 0.1)
M = length(temp)

### Obtain the seeds from the initialization run
seeds <- c()
for (i in 1:length(U_sim_PT_reps)){
  seeds[i] <- U_sim_PT_reps[[i]][["seed"]]
}

# Create a list to store results from different runs
U_sim_DPM_reps <- list()

## Run DPM for each initialization from PT
for (s in 1:length(seeds)){
  
  ### Choose the iteration from BAST posterior output ###
  results <- U_sim_PT[[s]]
  
  groups_assign = results$cluster_out; groups_number = results$k_out
  # List to save co-clustering matrices
  W_cum <- matrix(0, nrow = n, ncol = n)
  
  # Get the mode of the team
  tbl <- table(groups_number)
  max_freq <- max(tbl)
  modes <- as.numeric(names(tbl[tbl == max_freq]))
  
  # Get the indices where teams_number equals the mode
  mode_indices <- which(groups_number == modes)
  
  # Run the loop to create co-clustering matrix for over the selected indices
  for (i in 1:length(mode_indices)) {
    # Binary matrix for groups and teams
    X <- table(sequence(length(groups_assign[mode_indices[i], ])), groups_assign[mode_indices[i], ])
    # Compute W such that w_ij = 1 if observation i and observation j share same team
    W <- X %*% t(X)
    W_cum <- W_cum + W
  }
  
  mean_matrix <- W_cum/length(mode_indices)
  
  # Calculate the vector with Frobenius norms
  norm = c()
  
  for (i in 1:length(mode_indices)) {
    # Binary matrix for groups and teams
    X <- table(sequence(length(groups_assign[mode_indices[i], ])), groups_assign[mode_indices[i], ])
    # Get the membership of observations in each team
    mat <- X %*% t(X)
    diff_matrix <- mat - mean_matrix
    norm[i] <- norm(diff_matrix, type = "F")
  }
  
  # Fins the minimum Frobenius norm index
  min_index <- mode_indices[which.min(norm)]
  groups_assign_out <- groups_assign[min_index, ]
  rm(W, W_cum, mean_matrix)
  
  
  # Initialize values for DP-RST
  
  mu = list() # list for initial values of mu
  sigmasq_y = list()
  mstgraph_lst = list()  # initial spanning trees
  
  clusters_bast <- groups_assign_out
  clusters_bast_number <- length(unique(clusters_bast))
  k_max = clusters_bast_number
  
  cluster = matrix(rep(clusters_bast), nrow = n, ncol = M, byrow = F) 
  cluster_means_matrix <- matrix(0, nrow = k_max, ncol = p)
  # Get the means of each cluster
  for (i in 1:k_max) {
    cluster_means_matrix[i,] <- colMeans(Y_std[c(which(clusters_bast == i)), , drop = FALSE])
  }
  
  # Initialize all the groups from one team
  teams = matrix(rep(1, k_max), ncol = M, nrow = k_max, byrow = F)
  mu_teams = list()
  # For each temperature, we have same initial cluster assignment
  for(m in 1:M) {
    mu[[m]] = cluster_means_matrix #rep(0, p)
    mu_teams[[m]] = colMeans(Y_std)
    sigmasq_y[[m]] = cov(Y_std) #COvariance matrix
    mstgraph_lst[[m]] = mstgraph0
  }
  
  init_val = list()
  init_val[['trees']] = mstgraph_lst
  init_val[['mu']] = mu 
  init_val[['cluster']] = cluster 
  init_val[['teams']] = teams
  init_val[['mu_teams']] = mu_teams
  init_val[['sigmasq_y']] = sigmasq_y
  
  hyperpar = list()
  hyperpar[['sigmasq_mu']] = (0.5/(2*sqrt(1)))^2 
  hyperpar[['lambda_s']] = diag(1, p) 
  hyperpar[['nu']] = n-p+1
  hyperpar[['M']] = M
  hyperpar[['k_max']] = k_max
  hyperpar[['j_max']] = k_max
  
  # MCMC parameters
  # number of posterior samples = (MCMC - BURNIN) / THIN
  MCMC = 5000 # MCMC iterations
  BURNIN = 3000 # burnin period length
  THIN = 5     # thinning intervals
  
  hyperpar[['alpha']] = 0.5 # Concentration parameter for the DPM model
  
  dynamic_part <- "U_sim_rep" # To save backup files
  
  U_sim_DPM = DP.RST(Y_std, graph0, init_val, hyperpar, temp, 
                                  MCMC, BURNIN, THIN, 
                                  PT = TRUE, seed = seeds[s],
                                  backup_d = sprintf("./Simulations/U-shape_6k/%s_backup.RData", dynamic_part))
  U_sim_DPM$seed <- seeds[s]
  U_sim_DPM_reps[[s]] <- U_sim_DPM
  print("Loop is done")
} #end of loop over seeds

## Save the results of the repetitive runs
save(U_sim_DPM_reps, file = "./Simulations/U-shape_6k/U_sim_DPM_30reps.RData")
