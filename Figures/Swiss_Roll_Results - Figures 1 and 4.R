# Load necessary libraries
library(ggplot2)
library(ggpubr)
library(patchwork)
library(pdftools)
library(magick)
library(mclust)
library(SingleCellExperiment)
library(BayesSpace)
library("SC.MEB")
library(mvtnorm)
library(GiRaF)
library("DR.SC")
library(Seurat)
library(readr)
library(scCustomize)
library(data.table)
library(png)
library(gridExtra)
library(grid)
library(paletteer)
library(rprojroot)

# Set working directory to the root of the R project
project_root <- find_root(is_rstudio_project)
setwd(project_root)
getwd()

# Set seed
set.seed(7303)

# For debugging purposes
options(error=recover)

# Load custom functions
source("./Codes/Custom Functions/ComplexDomainFun.R")
source("./Codes/Custom Functions/plotting_fun.R")

# Load pre-processed data
load("./Real Data/swiss_roll_wt_muscle_finaltouches1.RData")
load("./Real Data/swiss_roll_wt_muscle_boundary.RData")

# Extract first three principal components and spatial coordinates
Y_sample = swiss_roll_wt_muscle_finaltouches1[, 1:3]
loc = swiss_roll_wt_muscle_finaltouches1[, 4:5]
bnd = boundary
# Extract hystological labels
z = swiss_roll_wt_muscle_finaltouches1$z
z_man = swiss_roll_wt_muscle_finaltouches1$z_man

n = nrow(Y_sample) # number of observations
p = ncol(Y_sample) # number of variables in Y

# Standardize coordinates, Y values
coords <- apply(loc, 2, scale)
Y_std <- apply(Y_sample, 2, scale)

n = nrow(Y_sample) # number of observations
p = ncol(Y_sample) # number of variables in Y

# Standardize spatial coordinates and principal component values
coords <- apply(loc, 2, scale)
Y_std <- apply(Y_sample, 2, scale)

# Standardize the boundary to match the scaled coordinates
bnd_scaled = list()
bnd_scaled$x <- (bnd$x - mean(loc[,1]))/sd(loc[,1])
bnd_scaled$y <- (bnd$y - mean(loc[,2]))/sd(loc[,2])

# Add small random variable to the scaled boundary to avoid duplicates
bnd_scaled_r = list()
bnd_scaled_r$x <- bnd_scaled$x + c(0, runif(length(bnd_scaled$x)-2, min = 1e-04, max = 5e-04), 0)
bnd_scaled_r$y <- bnd_scaled$y + c(0, runif(length(bnd_scaled$y)-2, min = 1e-04, max = 5e-04),0)

# Create dataframe for plotting
df_subset <- data.frame(coords, Y_std, z, z_man)

true_hystology_labels <- ggplot() + 
  geom_boundary(bnd_scaled_r) +
  geom_point(aes(x = x, y = y, col = z_man), data = df_subset) +
  labs(x = 'Coordinate X', y = 'Coordinate Y', colour = "True label") + 
  ggtitle('Ground Truth Labels
          (hystology)') +
  theme(legend.position = "bottom")
true_hystology_labels

true_groupped_labels <- ggplot() + 
  geom_boundary(bnd_scaled_r) +
  geom_point(aes(x = x, y = y, col = z), data = df_subset) +
  labs(x = 'Coordinate X', y = 'Coordinate Y', colour = "True label") + 
  ggtitle('Ground Truth Labels
          (manually groupped)') +
  theme(legend.position = "bottom")
true_groupped_labels

# Create a reference vector with values from 1 to 9
reference_vector <- c(1:9)

# Define a function to reorder a cluster assignments vector based on a reference vector
reorder_based_on_reference <- function(cluster_vector, reference_vector) {
  # Create an empty vector to store the reordered values
  reordered_vector <- numeric(length(cluster_vector))
  # Loop through each element in the cluster_vector
  for (i in seq_along(cluster_vector)) {
    # Find the position of the current element in the reference_vector
    position <- match(cluster_vector[i], reference_vector)
    # Assign the position to the reordered_vector
    reordered_vector[i] <- position
  }
  return(reordered_vector)
}


### Load DP-RST output
load("./Real Data/Swiss_Roll_DPM_OutputOnly.RData")

results = Swiss_Roll_Robb_DPM

##### CHOOSE THE ITERATION #####
groups_assign = results$cluster_out
teams_assign = results$teams_out
teams_number = results$j_teams_out

# List to save co-clustering matrices
W_cum <- matrix(0, nrow = n, ncol = n)

# Get the mode of the team
tbl <- table(teams_number)
max_freq <- max(tbl)
modes <- as.numeric(names(tbl[tbl == max_freq]))

# Get the indices where teams_number equals the mode
mode_indices <- which(teams_number == modes)

# Run the loop to create co-clustering matrix for over the selected indices
for (i in 1:length(mode_indices)) {
  # Binary matrix for groups and teams
  X <- table(sequence(length(groups_assign[mode_indices[i], ])), groups_assign[mode_indices[i], ])
  Z <- table(sequence(length(teams_assign[[mode_indices[i]]])), teams_assign[[mode_indices[i]]])

  # Get the membership of observations in each team
  obs_in_teams <- X %*% Z

  # Compute W such that w_ij = 1 if observation i and observation j share same team
  W <- obs_in_teams %*% t(obs_in_teams)

  W_cum <- W_cum + W
}

mean_matrix <- W_cum/length(mode_indices)

# Calculate the vector with Frobenius norms
norm = c()

for (i in 1:length(mode_indices)) {
  # Binary matrix for groups and teams
  X <- table(sequence(length(groups_assign[mode_indices[i], ])), groups_assign[mode_indices[i], ])
  Z <- table(sequence(length(teams_assign[[mode_indices[i]]])), teams_assign[[mode_indices[i]]])

  # Get the membership of observations in each team
  obs_in_teams <- X %*% Z
  # Compute W such that w_ij = 1 if observation i and observation j share same team
  mat <- obs_in_teams %*% t(obs_in_teams)

  diff_matrix <- mat - mean_matrix
  norm[i] <- norm(diff_matrix, type = "F")
}

# Find the minimum Frobenius norm index
min_index <- mode_indices[which.min(norm)] 

groups_assign_out <- groups_assign[min_index, ] # Cluster assignments for the crude spatial partition
teams_assign_out <- teams_assign[[min_index]] # Cluster assignments for the refined partition

rm(W, W_cum, mean_matrix)


#### Plot the crude spatial and refined partitions #####
# Determine the number of unique groups in the crude spatial partition output
k_groups_out <- length(unique(groups_assign_out))

# Calculate the accuracy of the crude spatial partition by comparing it with the true labels (df_subset$z)
accur_groups = adjustedRandIndex(df_subset$z, groups_assign_out)
accur_groups # 0.2024632

### Plot the crude spatial partition
groups_plotted <- groups_plot(coords, groups_assign_out, bnd_scaled_r) +
  theme(legend.position = "none") +
  ggtitle(paste("Estimated Groups, k = ", k_groups_out, ". Spatial Clustering. 
                ARI (gr.) = ", 
                round(accur_groups, 3), sep = ""))
groups_plotted

# Plot the teams of observations (clusters in the refined partition)
k_teams_out <- length(unique(teams_assign_out))

#### Calculate the accuracy of the labeled observations #####
# Binary matrix for groups and teams
X <- table(sequence(length(groups_assign_out)), groups_assign_out)
Z <- table(sequence(length(teams_assign_out)), teams_assign_out)

# Get the membership of observations in each team
obs_in_teams <- X %*% Z
# Compute W such that w_ij = 1 if observation i and observation j share same team
obs_in_teams_vec <- obs_in_teams %*% sort(unique(teams_assign_out))

accur_teams = adjustedRandIndex(df_subset$z, obs_in_teams_vec) #
accur_teams # 0.285258

# Define a custom color palette for plotting
my_palette1 <- c("#FCDD23FF", "#F8B100FF","#CA697CFF", "#AB74CFFF", "#48439BFF", 
  "#D9A453FF", "#E3A8A3FF", "#6A47D8FF", "#D856A7FF")

# Reorder the partition based on a reference vector
DPM_partition <- reorder_based_on_reference(obs_in_teams_vec, reference_vector)

# Plot the refined partition (DP-RST partition)
DPM.RST_plot <- groups_plot(coords, DPM_partition, bnd_scaled_r, point_size = 3) +
  scale_color_manual(values = my_palette1) +
  geom_point(size = 6) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title.x = element_blank(),  # No x-axis title
    axis.title.y = element_blank(),  # No y-axis title
    plot.title = element_text(size = 35)
  ) + 
  ggtitle(expression(bold("DP-RST") ~ ", c = 9, " * bold("ARI = 0.285")))

DPM.RST_plot

# Save the refined partition plot as a PNG image with a transparent background (part of the Figure 4)
ggsave(filename = "./Plots/DPM-RST_plot.png", 
       plot = DPM.RST_plot, width = 8, height = 8, bg = "transparent")


##### Create density plots for each refined cluster #####
obs_teams_assign = DPM_partition
team_data <- as.data.frame(cbind(Y_sample, obs_teams_assign))

# Create density plot for PC1
plot_PC1 <- ggplot(team_data, aes(x = PC1, fill = as.factor(obs_teams_assign))) +
  scale_fill_manual(values = my_palette1) +
  geom_density(alpha = 0.5) +
  facet_grid(obs_teams_assign ~ ., scales = "free_y", switch = "y") +
  labs(y = "Cluster", x = "PC1", fill = "Cluster") +
  theme_minimal() +
  theme(strip.text.y.left = element_text(size = 12, angle = 0),
        axis.title.y = element_text(size = 14), # Add y-axis title
        axis.text.y = element_blank(),  # Remove y-axis numbers
        axis.ticks.y = element_blank(), # Remove y-axis ticks
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none")

# Create density plot for PC2
plot_PC2 <- ggplot(team_data, aes(x = PC2, fill = as.factor(obs_teams_assign))) +
  scale_fill_manual(values = my_palette1) +
  geom_density(alpha = 0.5) +
  facet_grid(obs_teams_assign ~ ., scales = "free_y", switch = "y") +
  labs(x = "PC2", fill = "Cluster") +
  theme_minimal() +
  theme(strip.text.y.left = element_blank(),
        axis.title.y = element_blank(), # Remove y-axis title
        axis.text.y = element_blank(),  # Remove y-axis numbers
        axis.ticks.y = element_blank(), # Remove y-axis ticks
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none")

plot_PC3 <- ggplot(team_data, aes(x = PC3, fill = as.factor(obs_teams_assign))) +
  scale_fill_manual(values = my_palette1) +
  geom_density(alpha = 0.5) +
  facet_grid(obs_teams_assign ~ ., scales = "free_y", switch = "y") +
  labs(x = "PC3", fill = "Cluster") +
  theme_minimal() +
  theme(strip.text.y.left = element_blank(),
        axis.title.y = element_blank(), # Remove y-axis title
        axis.text.y = element_blank(),  # Remove y-axis numbers
        axis.ticks.y = element_blank(), # Remove y-axis ticks
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none")

# Arrange plots side by side
hist_DPM <- grid.arrange(plot_PC1, plot_PC2, plot_PC3, ncol = 3)
hist_DPM

# Save the arranged density plots as a PNG image (part of the Figure 4)
ggsave("./Plots/hist_DPM.png",
       plot = hist_DPM, width = 7.97, height = 5.35,
       units = "in", dpi = 300, bg = "transparent")

##### Plot the crude spatial partition for Figure 1 #####
DPM_teams_plot <- groups_plot(coords, DPM_partition, bnd_scaled_r, point_size = 3) +
  scale_color_manual(values = my_palette1) +
  geom_point(size = 6) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title.x = element_blank(),  # No x-axis title
    axis.title.y = element_blank(),  # No y-axis title
    plot.title = element_text(size = 35, hjust = 0.5)
  ) + 
  ggtitle(expression(bold("Refined Partition")))
DPM_teams_plot

# Define palette for the crude spatial partition
my_palette2 <- c("#FCDD23FF", "#F8B100FF", "#CA697CFF", "#AB74CFFF", "#48439BFF", 
  "#D9A453FF", "#E3A8A3FF", "#6A47D8FF", "#D856A7FF", "#3B83BDFF", "#FF7F50FF", 
  "#9ACD32FF", "#D2691EFF", "#1E90FFFF"  )

DPM_groups_plot <- groups_plot(coords, groups_assign_out, bnd_scaled_r, point_size = 3) +
  scale_color_manual(values = my_palette2) +
  geom_point(size = 6) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title.x = element_blank(),  # No x-axis title
    axis.title.y = element_blank(),  # No y-axis title
    plot.title = element_text(size = 35, hjust = 0.5)
  ) + 
  ggtitle(expression(bold("Crude Spatial Partition")))
DPM_groups_plot

ggsave(filename = "./Plots/DPM_teams_plot_Figure1.png", 
       plot = DPM_teams_plot, width = 8, height = 8, bg = "transparent")

ggsave(filename = "./Plots/DPM_groups_plot_Figure1.png", 
       plot = DPM_groups_plot, width = 8, height = 8, bg = "transparent")

##### COMPARE WITH OTHER METHODS #####

### BayesSpace ###
data_dir <- './Real Data'
expression_matrices <- Read10X_h5_GEO(data_dir = data_dir)
names(loc) <- c("row", "col")
# To create object from single file
seurat_object = CreateSeuratObject(counts = expression_matrices[[1]], meta.data = loc)
# Convert to SingleCellExperiment object
sce <- as.SingleCellExperiment(seurat_object)

log_count <- log(sce@assays@data@listData$counts + 1) 
M1 <- as(log_count, "CsparseMatrix") #"dgCMatrix"

sce@assays@data@listData$logcounts <- M1
rm(log_count)

dec <- scran::modelGeneVar(sce)
top <- scran::getTopHVGs(dec, n = 2000)

sce <- scater::runPCA(sce, subset_row=top, ncomponents = 3)

## Add BayesSpace metadata
sce <- BayesSpace::spatialPreprocess(sce, platform="Visium", skip.PCA=TRUE)

pca_reduced <- sce@int_colData@listData$reducedDims@listData

sce <- spatialCluster(sce, q=5, d=3, platform='Visium', nrep=10000, gamma=3)

BS.df <- data.frame(x = sce@colData$row, y = sce@colData$col, BS_z = sce$spatial.cluster)

true_labels_df <- data.frame(loc, df_subset$z_man, df_subset$z)
names(true_labels_df)[1:2] <- c("x", "y")

BS_merged_df <- merge(BS.df, true_labels_df, 
                      by = c("x", "y"), 
                      all = TRUE)

BS_accur = adjustedRandIndex(BS_merged_df$df_subset.z, BS_merged_df$BS_z)
BS_accur #0.1181117

# Reorder the partition based on a reference vector
BS_partition <- reorder_based_on_reference(BS_merged_df$BS_z, reference_vector)

BS_plot <- groups_plot_wB(BS_merged_df[, 1:2], BS_partition, point_size = 3) +
  scale_color_manual(values = my_palette1) +
  theme_minimal(base_size = 20) +
  geom_point(size = 6) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title.x = element_blank(),  # No x-axis title
    axis.title.y = element_blank(),  # No y-axis title
    plot.title = element_text(size = 35)
  ) + 
  ggtitle(expression(bold("BayesSpace") ~ ", c = 5, " * bold("ARI = 0.118")))

BS_plot

# Save the BayesSpace partition plot as a PNG image (part of Figure 4)
ggsave(filename = "./Plots/BS_plot.png", 
       plot = BS_plot, width = 8, height = 8, bg = "transparent")


### SC.MEB ###
y = as.matrix(Y_sample[, 1:3])

### Set the parameters
#platform = "ST"
beta_grid = seq(0, 5, 0.2) #Vector specifying the smoothness of Random Markov Field
K_set= 2:10 #Numbers of mixture components
parallel=TRUE #Logical value specifing the run the model in parallel
num_core = 3
PX = TRUE #Logical value for paramter expansion in EM algorithm
maxIter_ICM = 10 #Maximum iteration of ICM algorithm
maxIter = 50 #Maximum iteration of EM algorithm

#### Calculating the neighborhood from coordinates
Adj_sp <- getneighborhood_fast(as.matrix(loc), cutoff = 1.2)

#### Run the SC-MEB in parallel
SC.MEB_fit = SC.MEB(y, Adj_sp, beta_grid = beta_grid, K_set = K_set, 
                    parallel = parallel, num_core = num_core, PX = PX, 
                    maxIter_ICM = maxIter_ICM, maxIter = maxIter)

#### Selecting the number of clusters using BIC
selectKPlot(SC.MEB_fit, K_set = K_set, criterion = "BIC")

#### Selecting the number of clusters using Modified BIC
selectKPlot(SC.MEB_fit, K_set = K_set, criterion = "MBIC")

# Plot the clustering result using the optimal number of clusters
out = selectK(SC.MEB_fit, K_set = K_set, criterion = "MBIC")
ClusterPlot(out, loc)

SC.MEB_df <- data.frame(loc, out[["best_K_label"]], df_subset$z, df_subset$z_man)

SC.MEB_accur = adjustedRandIndex(df_subset$z, SC.MEB_fit[,4]$x)
SC.MEB_accur # 0.1111712

# Create a dataframe with the clustering results for 5 clusters (the truth)
SC.MEB_df <- data.frame(loc, c(SC.MEB_fit[,4]$x), df_subset$z, df_subset$z_man)

# Reorder the partition based on a reference vector
SC.MEB_partition <- reorder_based_on_reference(SC.MEB_df$c.SC.MEB_fit...4..x., reference_vector)

SC.MEB_plot <- groups_plot_wB(SC.MEB_df[, 1:2], SC.MEB_partition, point_size = 3) +
  scale_color_manual(values = my_palette1) +
  theme_minimal(base_size = 20) +
  geom_point(size = 6) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title.x = element_blank(),  # No x-axis title
    axis.title.y = element_blank(),  # No y-axis title
    plot.title = element_text(size = 35)
  ) + 
  ggtitle(expression(bold("SC-MEB") ~ ", c = 5, " * bold("ARI = 0.112")))

SC.MEB_plot

# Save the SC-MEB partition plot as a PNG image (part of Figure 4)
ggsave(filename = "./Plots/SC.MEB_plot.png", 
       plot = SC.MEB_plot, width = 8, height = 8, bg = "transparent")


### DR.SC ###
### Create Seurat object ###
# Create count matrix
expression_matrix <- Read10X_h5_GEO(data_dir = data_dir)

colnames(loc) <- c("row", "col")

# Create Seurat object
seurat_object = CreateSeuratObject(counts = expression_matrix, meta.data = loc)
head(seurat_object)

# standard log-normalization
Swiss <- NormalizeData(seurat_object, verbose = F)
# choose 500 highly variable features
seu <- FindVariableFeatures(Swiss, nfeatures = 2000, verbose = F)
# Fit DR-SC model using 500 highly variable features
seu_fit <- DR.SC(seu, K=5, q = 3, platform = 'Other_SRT', verbose=T)
# plot the obtained partition
spatialPlotClusters(seu_fit)


DR.SC_labels <- data.frame(seu_fit$row, seu_fit$col, seu_fit$spatial.drsc.cluster)

names(DR.SC_labels)[1:3] <- c("x", "y", "z_DR.SC")
DR.SC_merged_df <- merge(DR.SC_labels, true_labels_df, 
                         by = c("x", "y"), 
                         all = TRUE)

DR.SC_accur = adjustedRandIndex(DR.SC_merged_df$df_subset.z, DR.SC_merged_df$z_DR.SC)
DR.SC_accur #0.2031965

# Reorder the partition based on a reference vector
DR.SC_partition <- reorder_based_on_reference(DR.SC_merged_df$z_DR.SC, reference_vector)

DR.SC_plot <- groups_plot_wB(DR.SC_merged_df[, 1:2], DR.SC_partition, point_size = 3) +
  scale_color_manual(values = my_palette1) +
  theme_minimal(base_size = 20) +
  geom_point(size = 6) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title.x = element_blank(),  # No x-axis title
    axis.title.y = element_blank(),  # No y-axis title
    plot.title = element_text(size = 35)
  ) + 
  ggtitle(expression(bold("DR-SC") ~ ", c = 5, " * bold("ARI = 0.203")))

DR.SC_plot

# Save the DR-SC partition plot as a PNG image (part of Figure 4)
ggsave(filename = "./Plots/DR.SC_plot.png", 
       plot = DR.SC_plot, width = 8, height = 8, bg = "transparent")


##### k-means#####
kmeans_res <- kmeans(df_subset[, 3:5], 5)$cluster
adjustedRandIndex(df_subset$z, kmeans_res) #0.09389801

# Reorder the partition based on a reference vector
kmeans_partition <- reorder_based_on_reference(kmeans_res, reference_vector)

kmeans_plot <- groups_plot_wB(coords, kmeans_partition, point_size = 3) +
  scale_color_manual(values = my_palette1) +
  theme_minimal(base_size = 20) +
  geom_point(size = 6) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title.x = element_blank(),  # No x-axis title
    axis.title.y = element_blank(),  # No y-axis title
    plot.title = element_text(size = 35)
  ) + 
  ggtitle(expression(bold("k-means") ~ ", c = 5, " * bold("ARI = 0.094")))

kmeans_plot

# Save the k-means partition plot as a PNG image (part of Figure 4)
ggsave(filename = "./Plots/kmeans_plot.png", 
       plot = kmeans_plot, width = 8, height = 8, bg = "transparent")


