library(DP.RST)
library(ggplot2)
library(igraph)
library(gridExtra)
library(paletteer)

# set seed
set.seed(9362)

setwd("~/Desktop/DP-RST-workload")

#-------------------------------------------------------------------------------

### Load Data ###
load("Sensitivity/Swiss_Roll_sim_3k_data_1.RData")
load("Sensitivity/SR_sim_PT_30reps.RData")

# Extract first three PCAs
Y_sample = sim_data[, 5:7]
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

#-------------------------------------------------------------------------------

##### SET UP INITIAL VALUES #####
# Get mesh and triangulation
mesh = gen2dMesh(coords, bnd_scaled)
graph0 = constrainedDentri(n, mesh)
E(graph0)$eid = c(1:ecount(graph0))  # edge id
V(graph0)$vid = c(1:vcount(graph0))  # vertex id
mstgraph = mst(graph0)  # initial spanning tree
graph0 = delete_edge_attr(graph0, 'weight')
mstgraph0 = delete_edge_attr(mstgraph, 'weight')

#-------------------------------------------------------------------------------

# Define the "temperatures"
temp = seq(0.1:1, by = 0.1)
# Now M is not the number of weak learners but the number of temperatures
M = length(temp)

#-------------------------------------------------------------------------------

seed <- SR_sim_PT_reps[[1]][["seed"]]

#### Obrain the initialization from BAST ###
results <- SR_sim_PT_reps[[1]]

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


#-------------------------------------------------------------------------------

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
init_val[['mstgraph_lst']] = mstgraph_lst
init_val[['mu']] = mu
init_val[['cluster']] = cluster
init_val[['teams']] = teams
init_val[['mu_teams']] = mu_teams
init_val[['sigmasq_y']] = sigmasq_y

hyperpar = list()
hyperpar[['sigmasq_mu']] = (0.5/(2*sqrt(1)))^2 
hyperpar[['lambda_s']] = diag(1, p) # lambda_s
hyperpar[['nu']] = p #n-p+1
hyperpar[['M']] = M
hyperpar[['temp']] = temp
hyperpar[['k_max']] = k_max
hyperpar[['j_max']] = k_max

#-------------------------------------------------------------------------------

# MCMC parameters
# number of posterior samples = (MCMC - BURNIN) / THIN
MCMC = 3000 # MCMC iterations
BURNIN = 2000 # burnin period length
THIN = 5     # thinning intervals

#-------------------------------------------------------------------------------

SR_sim_SensAlpha <- list()

alpha_grid <- c(0.05, 0.1, 0.2, 0.5, 1)

## Run DPM for each initialization from PT
for (s in 1:length(alpha_grid)){

    hyperpar[['alpha']] = alpha_grid[s]

    SR_sim_DPM = DP.RST(Y_std, graph0, init_val, hyperpar,
                        MCMC, BURNIN, THIN,
                        PT = TRUE, PT_diff = 0.1, seed = seed)

    SR_sim_DPM$seed <- seed
    SR_sim_SensAlpha[[s]] <- SR_sim_DPM

    print("Loop is done")
}

## Save the results of the repetitive runs
save(SR_sim_SensAlpha, file = "Sensitivity/SR_sim_SensAlpha_HyperSetVersion2.RData")

#-------------------------------------------------------------------------------
##### EXPLORE THE RESULTS #####

load("~/Desktop/DP-RST-workload/Sensitivity/SR_sim_SensAlpha_HyperSetVersion2.RData")

z = sim_data$cluster

# Initialize lists to store plots for each simulation output
hist_list <- list()
mode_list <- list()

my_palette1 <- paletteer::paletteer_d("IslamicArt::samarqand2")

# Loop over each simulation output in the list
for(i in seq_along(SR_sim_SensAlpha)) {

  # Extract the output for the current alpha
  output <- SR_sim_SensAlpha[[i]]

  # ---------- Histogram Plot ----------
  # Compute frequency table of refined clusters (assumed stored as "j_teams_out")
  teams_freq <- as.data.frame(table(output[["j_teams_out"]]))

  # Create histogram plot for the number of refined clusters, labeling with the current alpha value
  hist_list[[i]] <- ggplot(teams_freq, aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity", color = "black", fill = "lightblue") +
    labs(title = bquote(alpha == .(alpha_grid[i])),
         x = "# of Refined Clusters", y = "Frequency") +
    theme_classic()

  # ---------- Compute Partitions and Accuracy ----------
  # Get the best partitions using two methods
  mode_based_partition <- partition(DP.RST_output = output, method = "mode_based", batch_size = 100)
  frobenius_partition  <- partition(DP.RST_output = output, method = "frobenius", batch_size = 100)

  # Compute the binary membership matrices and then team membership vectors
  # For mode-based method:
  X_mode <- table(seq_along(mode_based_partition$groups_partition), mode_based_partition$groups_partition)
  Z_mode <- table(seq_along(mode_based_partition$teams_partition), mode_based_partition$teams_partition)
  obs_in_teams_mode <- as.matrix(X_mode) %*% as.matrix(Z_mode)
  obs_in_teams_vec_mode <- obs_in_teams_mode %*% sort(unique(mode_based_partition$teams_partition))

  # For frobenius method:
  X_frob <- table(seq_along(frobenius_partition$groups_partition), frobenius_partition$groups_partition)
  Z_frob <- table(seq_along(frobenius_partition$teams_partition), frobenius_partition$teams_partition)
  obs_in_teams_frob <- as.matrix(X_frob) %*% as.matrix(Z_frob)
  obs_in_teams_vec_frob <- obs_in_teams_frob %*% sort(unique(frobenius_partition$teams_partition))

  # Create a data frame for ARI calculations
  Gut_df <- data.frame(
    res_groups_1 = mode_based_partition$groups_partition,
    res_teams_1  = as.vector(obs_in_teams_vec_mode),
    res_groups_2 = frobenius_partition$groups_partition,
    res_teams_2  = as.vector(obs_in_teams_vec_frob),
    z            = z
  )
  Gut_df <- subset(Gut_df, !is.na(z))

  # Calculate Adjusted Rand Index (ARI) for each refined partition
  accur_teams_1 <- mclust::adjustedRandIndex(Gut_df$res_teams_1, Gut_df$z)
  accur_teams_2 <- mclust::adjustedRandIndex(Gut_df$res_teams_2, Gut_df$z)

  # ---------- Create Refined Partition Plots ----------
  # Mode-based plot
  mode_list[[i]] <- teams_plot(coords, mode_based_partition$teams_partition,
                               mode_based_partition$groups_partition, bnd_scaled) +
    theme(legend.position = "none") +
    scale_color_manual(values = my_palette1) +
    ggtitle(paste("Mode-based partition, k =",
                  length(unique(mode_based_partition$teams_partition)),
                  ", ARI =", round(accur_teams_1, 3)))

  # Frobenius plot
  frobenius_list[[i]] <- teams_plot(coords, frobenius_partition$teams_partition,
                                    frobenius_partition$groups_partition, bnd_scaled) +
    theme(legend.position = "none") +
    scale_color_manual(values = my_palette1) +
    ggtitle(paste("Frobenius Partition, k =",
                  length(unique(frobenius_partition$teams_partition)),
                  ", ARI =", round(accur_teams_2, 3)))
}

# ---------- Arrange All Plots in a Grid ----------
# Combine the three lists into one vector of plots
all_plots <- c(hist_list, mode_list, frobenius_list)

# Arrange in a grid:
# - Number of columns = number of alpha values (i.e. length of SR_sim_SensAlpha)
# - 3 rows: row 1 (histograms), row 2 (mode-based plots), row 3 (frobenius plots)
grid.arrange(grobs = all_plots, ncol = length(SR_sim_SensAlpha), nrow = 3)

# Create the arranged grob
final_plot <- arrangeGrob(grobs = all_plots, ncol = length(SR_sim_SensAlpha), nrow = 3)

# Save as PDF in the "Sensitivity" folder
ggsave("Sensitivity/SR_sim_SensAlpha_HyperSetVersion2.pdf", final_plot, width = 20, height = 10)


#-------------------------------------------------------------------------------
##### PROCESS AND PLOT THE RESULTS FOR PAPER #####

load("~/Desktop/DP-RST-workload/Sensitivity/SR_sim_SensAlpha_HyperSetVersion2.RData")

z = sim_data$cluster
alpha_grid <- c(0.05, 0.1, 0.2, 0.5, 1)

# Initialize lists to store plots for each simulation output
hist_list <- list()
mode_list <- list()

my_palette1 <- paletteer::paletteer_d("IslamicArt::samarqand2")

# Loop over each simulation output in the list
for(i in seq_along(SR_sim_SensAlpha)) {
  
  # Extract the output for the current alpha
  output <- SR_sim_SensAlpha[[i]]
  
  # ---------- Histogram Plot ----------
  # Compute frequency table of refined clusters (assumed stored as "j_teams_out")
  teams_freq <- as.data.frame(table(output[["j_teams_out"]]))
  
  # Create histogram plot for the number of refined clusters, labeling with the current alpha value
  hist_list[[i]] <- ggplot(teams_freq, aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity", color = "black", fill = "lightblue") +
    labs(title = bquote(alpha == .(alpha_grid[i])),
         x = "# of Refined Clusters", y = "Frequency") +
    theme_classic()
  
  # ---------- Compute Partitions and Accuracy ----------
  # Get the best partitions using two methods
  mode_based_partition <- partition(DP.RST_output = output, method = "mode_based", batch_size = 100)
  
  # Compute the binary membership matrices and then team membership vectors
  # For mode-based method:
  X_mode <- table(seq_along(mode_based_partition$groups_partition), mode_based_partition$groups_partition)
  Z_mode <- table(seq_along(mode_based_partition$teams_partition), mode_based_partition$teams_partition)
  obs_in_teams_mode <- as.matrix(X_mode) %*% as.matrix(Z_mode)
  obs_in_teams_vec_mode <- obs_in_teams_mode %*% sort(unique(mode_based_partition$teams_partition))
  
  # Create a data frame for ARI calculations
  Gut_df <- data.frame(
    res_groups_1 = mode_based_partition$groups_partition,
    res_teams_1  = as.vector(obs_in_teams_vec_mode),
    z            = z
  )
  Gut_df <- subset(Gut_df, !is.na(z))
  
  # Calculate Adjusted Rand Index (ARI) for each refined partition
  accur_teams_1 <- mclust::adjustedRandIndex(Gut_df$res_teams_1, Gut_df$z)
  
  # ---------- Create Refined Partition Plots ----------
  # Mode-based plot
  mode_list[[i]] <- teams_plot(coords, mode_based_partition$teams_partition,
                               mode_based_partition$groups_partition, bnd_scaled) +
    theme_classic() +  # Remove grey background
    theme(legend.position = "none",
          plot.title = element_text(face = "italic", size = 12)) +  # Italicize title
    scale_color_manual(values = my_palette1) +
    ggtitle(paste("c =", length(unique(mode_based_partition$teams_partition)),
                  ", ARI =", round(accur_teams_1, 3)))
  
}

# ---------- Arrange All Plots in a Grid ----------
# Combine the three lists into one vector of plots
all_plots <- c(hist_list, mode_list)

# Arrange in a grid:
# - Number of columns = number of alpha values (i.e. length of SR_sim_SensAlpha)
# - 3 rows: row 1 (histograms), row 2 (mode-based plots), row 3 (frobenius plots)
grid.arrange(grobs = all_plots, ncol = length(SR_sim_SensAlpha), nrow = 2)

# Create the arranged grob
final_plot <- arrangeGrob(grobs = all_plots, ncol = length(SR_sim_SensAlpha), nrow = 2)

# Save as PDF in the "Sensitivity" folder
ggsave("Sensitivity/Sensetivity_Alpha.pdf", final_plot, width = 15, height = 6)


