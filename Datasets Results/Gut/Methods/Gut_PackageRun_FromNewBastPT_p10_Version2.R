start.time_file <- Sys.time()

# options(repos = c(CRAN = "https://cran.rstudio.com"))

# Set a writable library path
my_lib_path <- "/Rlibs"
# if (!dir.exists(my_lib_path)) {
#   dir.create(my_lib_path)
# }
.libPaths(my_lib_path)

library(DP.RST)
library(igraph)
library(dplyr)
print("Libraries are loaded")

# set seed
set.seed(729)

gut_df_wt_muscle = readRDS("/GSE169749_RAW-2/NewRuns_2025/gut_df_wt_muscle.rds")
load("/GSE169749_RAW-2/swiss_roll_wt_muscle_boundary.RData")
load("/GSE169749_RAW-2/NewRuns_2025/Outputs_BAST/Gut_BAST_PT_OutputOnly.RData")

# Extract first three PCAs
Y_sample = gut_df_wt_muscle[, 1:10]
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

# Create dataframe for plotting
df_subset <- data.frame(coords, Y_std)

#### SET UP INITIAL VALUES #####
results = Gut_BAST_PT

##### CHOOSE THE ITERATION #####
groups_assign = results$cluster_out; groups_number = results$k_out
# List to save co-clustering matrices
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


# Get mesh and triangulation
# mesh = gen2dMesh(coords, bnd_scaled)
load("/GSE169749_RAW-2/swiss_roll_wt_muscle_mesh.RData")
graph0 = constrainedDentri(n, mesh)
E(graph0)$eid = c(1:ecount(graph0))  # edge id
V(graph0)$vid = c(1:vcount(graph0))  # vertex id
graph0 = delete_edge_attr(graph0, 'weight')
# mstgraph0 = delete_edge_attr(init_partition$mst, 'weight')


# Define the "temperatures"
temp = seq(0.1:1, by = 0.1)
M = length(temp)


# Initialize parameters for DP-RST

mu = list() # list for initial values of mu
sigmasq_y = list()
mstgraph_lst = list()  # initial spanning trees

clusters_bast <- groups_assign_out
clusters_bast_number <- length(unique(clusters_bast))
k_max = clusters_bast_number
print(k_max)

cluster = matrix(rep(clusters_bast), nrow = n, ncol = M, byrow = F)
cluster_means_matrix <- matrix(0, nrow = k_max, ncol = p, byrow = F)
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
  mstgraph_lst[[m]] = results[["tree_out"]][[min_index]]
}

init_val = list()
init_val[['mstgraph_lst']] = mstgraph_lst
init_val[['mu']] = mu
init_val[['cluster']] = cluster
init_val[['teams']] = teams
init_val[['mu_teams']] = mu_teams
init_val[['sigmasq_y']] = sigmasq_y

hyperpar = list()
hyperpar[['sigmasq_mu']] = 0.3
hyperpar[['lambda_s']] = diag(1, p)
hyperpar[['nu']] = p
hyperpar[['M']] = M
hyperpar[['temp']] = temp
hyperpar[['k_max']] = k_max
hyperpar[['j_max']] = k_max

# MCMC parameters
# number of posterior samples = (MCMC - BURNIN) / THIN
MCMC = 25000 # MCMC iterations
BURNIN = 20000 # burnin period length
THIN = 1     # thinning intervals

hyperpar[['alpha']] = 0.5 # Concentration parameter for the DPM model

print("All the parameters are set")

start.time_mcmc <- Sys.time()
Gut_DP.RST_FromNewBastPT_p10_Version2 = DP.RST(Y_std, graph0, init_val, hyperpar,
                                                MCMC, BURNIN, THIN,
                                                PT = TRUE, PT_diff = 0.1, seed = 5284)
end.time_mcmc <- Sys.time()

mcmc_time <- end.time_mcmc - start.time_mcmc
mcmc_time

save(Gut_DP.RST_FromNewBastPT_p10_Version2, file = "/GSE169749_RAW-2/NewRuns_2025/Outputs_DP-RST/Gut_DP.RST_FromNewBastPT_p10_Version2_OutputOnly.RData")

end.time_file <- Sys.time()

file_time <- end.time_file - start.time_file
file_time
