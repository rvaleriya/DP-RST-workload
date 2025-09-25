library(DP.RST)
library(ggplot2)
library(igraph)
library(gridExtra)
library(paletteer)

# set seed
set.seed(9362)

# For debugging purposes
options(error=recover)

setwd("~/Desktop/DP-RST-workload")

#-------------------------------------------------------------------------------

### Load Data ###
load("Sensitivity/Swiss_Roll_sim_3k_data_1.RData")

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

# MCMC parameters
# number of posterior samples = (MCMC - BURNIN) / THIN
MCMC = 3000 # MCMC iterations
BURNIN = 2000 # burnin period length
THIN = 5     # thinning intervals

#-------------------------------------------------------------------------------

# Compute the geodesic distance matrix for the selected batch of nodes
geodesic_dist_matrix <- distances(graph0, v = V(graph0), to = V(graph0),
                                  weights = E(graph0)$weight)

# Replace Inf with a large value
geodesic_dist_matrix[is.infinite(geodesic_dist_matrix)] <- max(geodesic_dist_matrix[is.finite(geodesic_dist_matrix)], na.rm = TRUE) * 2

# Compute the Euclidean distance matrix for the PCA data
euclidean_dist_matrix <- as.matrix(dist(Y_std))

# Compute the weighted average of the geodesic and Euclidean distance matrices
alpha <- 1/2  # Adjust this value to control the weight of geodesic vs. Euclidean distance
combined_dist_matrix <- alpha * geodesic_dist_matrix + (1 - alpha) * euclidean_dist_matrix

# Perform hierarchical clustering on the combined distance matrix
combined_dist <- as.dist(combined_dist_matrix)
hc <- fastcluster::hclust(combined_dist, method = "ward.D2")

# Cut the tree to form clusters
num_clusters <- c(8, 12, 16, 20, 24)

clusters <- cutree(hc, k = num_clusters)

#-------------------------------------------------------------------------------

SR_sim_SensK <- list()

for (i in 1:length(num_clusters)){

  # Compute centroids for each cluster in the batch
  centroids_batch <- aggregate(Y_std, by = list(clusters[, i]), FUN = mean)

  # Store the centroids and cluster assignments for this batch
  centroids_batch <- as.matrix(centroids_batch[, -1])  # Remove cluster labels

  mu = list() # list for initial values of mu
  sigmasq_y = list()
  mstgraph_lst = list()  # initial spanning trees

  clusters_bast <- clusters[, i]
  clusters_bast_number <- length(unique(clusters_bast))
  k_max = clusters_bast_number

  cluster = matrix(rep(clusters_bast), nrow = n, ncol = M, byrow = F)

  # Initialize all the groups from one team
  teams = matrix(rep(1, k_max), ncol = M, nrow = k_max, byrow = F)
  mu_teams = list()
  # For each temperature, we have same initial cluster assignment
  for(m in 1:M) {
    mu[[m]] = centroids_batch #rep(0, p)
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

  hyperpar[['alpha']] = 0.1

  SR_sim_DPM = DP.RST(Y_std, graph0, init_val, hyperpar,
                      MCMC, BURNIN, THIN,
                      PT = TRUE, PT_diff = 0.1, seed = 141)

  SR_sim_SensK[[i]] <- SR_sim_DPM

  print("Loop is done")
}


## Save the results of the repetitive runs
save(SR_sim_SensK, file = "Sensitivity/SR_sim_SensK.RData")

#-------------------------------------------------------------------------------
##### EXPLORE THE RESULTS #####
z = sim_data$cluster

# Initialize lists to store plots for each simulation output
init_part_list <- list()
mode_list <- list()
frobenius_list <- list()

my_palette1 <- paletteer::paletteer_d("IslamicArt::samarqand2")

# Define the 25-color palette
custom_palette_25 <- c(
  "#1B1F8A", "#DA4167", "#2A9D8F", "#E76F51", "#264653",
  "#F4A261", "#8B5E3B", "#E9C46A", "#F94144", "#43AA8B",
  "#577590", "#FFB703", "#023047", "#9A031E", "#5F0F40",
  "#A7C957", "#D4A373", "#118AB2", "#8ECAE6", "#6A0572",
  "#E29578", "#2C6E49", "#E76F51", "#C77DFF", "#A52A2A"
)


# Loop over each simulation output in the list
for(i in seq_along(SR_sim_SensK)) {

  # Extract the output for the current alpha
  output <- SR_sim_SensK[[i]]

  # ---------- Initial Spatial Partition Plot ----------

  # Create histogram plot for the number of refined clusters, labeling with the current alpha value
  init_part_list[[i]] <- groups_plot(coords, output[["init_val"]][["cluster"]][,1], bnd_scaled) +
    theme(legend.position = "none") +
    scale_color_manual(values = custom_palette_25) +
    ggtitle(paste("Initial Spatial Partition, K =",
                  length(unique(output[["init_val"]][["cluster"]][,1]))))

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
    ggtitle(paste("Mode-based partition, C =",
                  length(unique(mode_based_partition$teams_partition)),
                  ", ARI =", round(accur_teams_1, 3)))

  # Frobenius plot
  frobenius_list[[i]] <- teams_plot(coords, frobenius_partition$teams_partition,
                                    frobenius_partition$groups_partition, bnd_scaled) +
    theme(legend.position = "none") +
    scale_color_manual(values = my_palette1) +
    ggtitle(paste("Frobenius Partition, C =",
                  length(unique(frobenius_partition$teams_partition)),
                  ", ARI =", round(accur_teams_2, 3)))
}

# ---------- Arrange All Plots in a Grid ----------
# Combine the three lists into one vector of plots
all_plots <- c(init_part_list, mode_list, frobenius_list)

# Arrange in a grid:
# - Number of columns = number of alpha values (i.e. length of SR_sim_SensK)
# - 3 rows: row 1 (initial partitions), row 2 (mode-based plots), row 3 (frobenius plots)
grid.arrange(grobs = all_plots, ncol = length(SR_sim_SensK), nrow = 3)

# Create the arranged grob
final_plot <- arrangeGrob(grobs = all_plots, ncol = length(SR_sim_SensK), nrow = 3)

# Save as PDF in the "Sensitivity" folder
ggsave("Sensitivity/SR_sim_SensK.pdf", final_plot, width = 20, height = 10)


#-------------------------------------------------------------------------------
##### PROCESS AND PLOT THE RESULTS FOR PAPER #####

load("./Sensitivity/SR_sim_SensK.RData")

z = sim_data$cluster

# Initialize lists to store plots for each simulation output
init_part_list <- list()
mode_list <- list()

my_palette1 <- paletteer::paletteer_d("IslamicArt::samarqand2")

# Define the 25-color palette
custom_palette_25 <- c(
  "#1B1F8A", "#DA4167", "#2A9D8F", "#E76F51", "#264653",
  "#F4A261", "#8B5E3B", "#E9C46A", "#F94144", "#43AA8B",
  "#577590", "#FFB703", "#023047", "#9A031E", "#5F0F40",
  "#A7C957", "#D4A373", "#118AB2", "#8ECAE6", "#6A0572",
  "#E29578", "#2C6E49", "#E76F51", "#C77DFF", "#A52A2A"
)


# Loop over each simulation output in the list
for(i in seq_along(SR_sim_SensK)) {
  
  # Extract the output for the current alpha
  output <- SR_sim_SensK[[i]]
  
  # ---------- Initial Spatial Partition Plot ----------
  
  # Create histogram plot for the number of refined clusters, labeling with the current alpha value
  init_part_list[[i]] <- groups_plot(coords, output[["init_val"]][["cluster"]][,1], bnd_scaled) +
    theme_classic() +  # Remove grey background
    theme(legend.position = "none", plot.title = element_text(face = "italic", size = 12)) +
    scale_color_manual(values = custom_palette_25) +
    ggtitle(bquote("Initial Spatial Partition, " * italic(k == .(length(unique(output[["init_val"]][["cluster"]][,1]))))))
  
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
    theme(legend.position = "none", plot.title = element_text(face = "italic", size = 12)) +
    scale_color_manual(values = my_palette1) +
    ggtitle(paste("k =", length(unique(output[["init_val"]][["cluster"]][,1])), ", c =", length(unique(mode_based_partition$teams_partition)),
                  ", ARI =", round(accur_teams_1, 3)))
  
}

# ---------- Arrange All Plots in a Grid ----------
# Combine the three lists into one vector of plots
# all_plots <- c(init_part_list, mode_list)
all_plots <- c(mode_list)

# Arrange in a grid:
# - Number of columns = number of alpha values (i.e. length of SR_sim_SensK)
# - 3 rows: row 1 (initial partitions), row 2 (mode-based plots), row 3 (frobenius plots)
grid.arrange(grobs = all_plots, ncol = length(SR_sim_SensK), nrow = 1)
# grid.arrange(grobs = all_plots, ncol = length(SR_sim_SensK), nrow = 2)

# Create the arranged grob
# final_plot <- arrangeGrob(grobs = all_plots, ncol = length(SR_sim_SensK), nrow = 2)
final_plot <- arrangeGrob(grobs = all_plots, ncol = length(SR_sim_SensK), nrow = 1)

# Save as PDF in the "Sensitivity" folder
# ggsave("Sensitivity/Sensetivity_SpatialClusters.pdf", final_plot, width = 15, height = 6)
ggsave("/Sensitivity/Sensetivity_SpatialClusters.pdf", final_plot, width = 15, height = 3)




