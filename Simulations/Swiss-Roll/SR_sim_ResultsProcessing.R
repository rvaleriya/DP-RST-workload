# Load necessary libraries
library(DP.RST)
library(ggplot2)
library(patchwork)
library(mclust)
library(dplyr)
library(tidyr)
library(cowplot)

#-------------------------------------------------------------------------------
# Set working directory to the root of the R project
setwd("~/Desktop/DP-RST-workload")

# set seed
set.seed(9362)

# For debugging purposes
options(error=recover)

### Load Data ###
# SR_sim_data <- readr::read_csv("./New_Simulations/SR_Extra/SR_sim_data.csv")
load("/Swiss_Roll_sim_3k_data_1.RData")
SR_sim_data = sim_data
#-------------------------------------------------------------------------------

# Extract first three PCAs
Y_sample = SR_sim_data[, 5:7]
loc = SR_sim_data[, 1:2]
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
### Function to plot individual simulation results for Figure 3 in paper

# Define pallettes for plots
my_palette <- paletteer::paletteer_d("ochRe::galah")
my_palette1 <- paletteer::paletteer_d("IslamicArt::samarqand2")

plot_points_paper <- function(coords, values, bnd_scaled, title = NULL) {

  # Preparing data for points
  pointsdata <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    Value = values  # Converting scaled values to a factor for discrete coloring
  )
  
  # Create the ggplot object with boundary and points
  plot <- ggplot() +
    geom_point(aes(x = x, y = y, color = as.factor(Value)), data = pointsdata, size = 4) +
    scale_color_viridis_c() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.title.x = element_blank(),  # Remove x-axis title
          axis.title.y = element_blank(),  # Remove y-axis title
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank()
          # axis.ticks = element_line(color = "black"),
          # axis.line = element_line(colour = "black", size = 0.5,
          #                          arrow = arrow(type = "closed", length = unit(0.08, "inches")))
    )
  
  return(plot)
}

#-------------------------------------------------------------------------------

# Define reference vector to assign same clusters to the same colors
reference_vector <- c(1, 2, 3)

# Function to reorder cluster assignments based on the reference vector
reorder_based_on_reference <- function(cluster_vector, reference_vector) {
  map <- match(reference_vector, cluster_vector)
  reordered_vector <- map[cluster_vector]
  return(reordered_vector)
}

#-------------------------------------------------------------------------------

# Plot the ground truth for the U-shape example
true_clusters <- plot_points_paper(coords, SR_sim_data$cluster, bnd_scaled) +
  geom_boundary(bnd_scaled) +
  scale_color_manual(values = my_palette, labels = c(expression(mu == -1), expression(mu == 0), expression(mu == 1))) +
  ggtitle(expression(bold("Truth"))) +
  theme(plot.title = element_text(size = 20),
        legend.position = "left",
        legend.title = element_blank(),
        legend.text = element_text(size = 24))
true_clusters

# Add label to the plot for Figure 3
true_clusters_with_label <- ggdraw() +
  draw_plot(true_clusters) +
  draw_label("(A)", x = 0, y = 0.98, size = 20, hjust = 0, vjust = 1)
true_clusters_with_label

ggsave("./New_Simulations/SR_Extra/SR_sim_true_clusters.png", 
       plot = true_clusters_with_label, 
       width = 6.5, height = 5, units = "in", dpi = 300, bg = "transparent")

#-------------------------------------------------------------------------------
##### DP-RST results #####
load("./New_Simulations/SR_Extra/SR_sim_DPM_NewNu_30reps.RData")

### Calculate the accuracy for each run ###
DPM_partition <- list()
DPM_accur = c()
clust_number <- c()

for (j in 1:30) {
  mode_based_partition = partition(DP.RST_output = SR_sim_DPM_reps_NewNu[[j]], 
                                   method = "mode_based", batch_size = 100)
  DPM_partition [[j]] <- mode_based_partition
  
  ## Accuracy of the refined partition for the mode_based method
  # Binary membership matrix for groups and teams
  X <- table(sequence(length(mode_based_partition$groups_partition)),
             mode_based_partition$groups_partition)
  Z <- table(sequence(length(mode_based_partition$teams_partition)),
             mode_based_partition$teams_partition)
  
  # Get the membership of observations in each team
  obs_in_teams <- X %*% Z
  # Compute W such that w_ij = 1 if observation i and observation j share same team
  obs_in_teams_vec_mode_based <- obs_in_teams %*% sort(unique(mode_based_partition$teams_partition))
  
  DPM_accur[j] = adjustedRandIndex(SR_sim_data$cluster, obs_in_teams_vec_mode_based)
  clust_number[j] = length(unique(obs_in_teams_vec_mode_based))
}

# Plot the accuracy over runs
plot(DPM_accur, type = "l")
summary(DPM_accur)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6512  0.6898  0.7282  0.7306  0.7622  0.8522 

table(clust_number)
# clust_number
# 3 
# 30 

## Get the cluster assignments of the highest accuracy
original_index = which.max(DPM_accur)

X <- table(sequence(length(DPM_partition[[original_index]]$groups_partition)), DPM_partition[[original_index]]$groups_partition)
Z <- table(sequence(length(DPM_partition[[original_index]]$teams_partition)), DPM_partition[[original_index]]$teams_partition)
obs_in_teams <- X %*% Z
DPM_partition <- obs_in_teams %*% sort(unique(DPM_partition[[original_index]]$teams_partition))

DPM_partition <- reorder_based_on_reference(DPM_partition, reference_vector)
se_DPM <- sd(DPM_accur) / sqrt(length(DPM_accur))

## Plot the partition
DPM_results <- plot_points_paper(coords, DPM_partition, bnd_scaled) +
  geom_boundary(bnd_scaled) +
  scale_color_manual(values = my_palette1) +
  ggtitle(bquote(bold("DP-RST, ") ~ italic("ARI = ") ~ .(round(max(DPM_accur), 3)) ~ " (" ~ .(round(se_DPM, 3)) ~ ")")) + 
  theme(plot.title = element_text(size = 20))
DPM_results

DPM_results_with_label <- ggdraw() +
  draw_plot(DPM_results) +
  draw_label("(i)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
DPM_results_with_label

ggsave("./New_Simulations/SR_Extra/SR_sim_DPM.png", 
       plot = DPM_results_with_label, 
       width = 5, height = 5, units = "in", dpi = 300, bg = "transparent")

#-------------------------------------------------------------------------------
### Make a histogram with the number of teams chosen over runs (for the Figure 1 in Supplementary Materials)
# Convert to a data frame
clust_number <- data.frame(clust_number = clust_number)

# Create the histogram
SR_hist_plot <- ggplot(clust_number, aes(x = clust_number)) +
  geom_histogram(binwidth = 0.1, fill = "#F0A8C0FF", color = "black") +
  theme_minimal() +
  labs(title = "Cluster Number Estimation in Swiss-Roll Simulation", 
       x = "Estimated number of clusters", 
       y = "Frequency") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  xlim(1, 5)  # Extend the x-axis from 1 to 5
SR_hist_plot

ggsave("./New_Simulations/SR_Extra/SR_hist_plot.png", SR_hist_plot, width = 6, height = 4)

#-------------------------------------------------------------------------------
##### BayesSpace results #####
load("/SR_sim_BayesSpace_30reps.RData")

BayesSpace_acc <- c()
for (i in 1:30){
  BayesSpace_acc[i] <- adjustedRandIndex(SR_sim_data$cluster, SR_sim_BayesSpace_reps[[i]]$spatial.cluster)
}

plot(BayesSpace_acc, type = "l")
summary(BayesSpace_acc)

BayesSpace_index <- which.max(BayesSpace_acc)

BS_partition <- reorder_based_on_reference(SR_sim_BayesSpace_reps[[BayesSpace_index]]$spatial.cluster, reference_vector)
se_BS <- sd(BayesSpace_acc) / sqrt(length(BayesSpace_acc))

BayesSpace_results <- plot_points_paper(coords, BS_partition, bnd_scaled) +
  geom_boundary(bnd_scaled) +
  scale_color_manual(values = my_palette1) +
  ggtitle(bquote(bold("BayesSpace, ") ~ italic("ARI = ") ~ .(round(max(BayesSpace_acc), 3)) ~ " (" ~ .(round(se_BS, 3)) ~ ")"))+
  theme(plot.title = element_text(size = 20))
BayesSpace_results

BayesSpace_results_with_label <- ggdraw() +
  draw_plot(BayesSpace_results) +
  draw_label("(ii)", x = -0.004, y = 0.99, size = 20, hjust = 0, vjust = 1)
BayesSpace_results_with_label

ggsave("./New_Simulations/SR_Extra/SR_sim_BayesSpace.png", 
       plot = BayesSpace_results_with_label, 
       width = 5.1, height = 5.1, units = "in", dpi = 300, bg = "transparent")

#-------------------------------------------------------------------------------

##### SC-MEB results #####
load("~/SR_sim_SC.MEB_30reps.RData")

SC.MEB_acc <- c()
for (i in 1:30){
  SC.MEB_acc[i] <- adjustedRandIndex(SR_sim_data$cluster, SR_sim_SC.MEB_reps[[i]][[10]])
}

plot(SC.MEB_acc, type = "l")
summary(SC.MEB_acc)

SC.MEB_index <- which.max(SC.MEB_acc)

SC.MEB_partition <- reorder_based_on_reference(SR_sim_SC.MEB_reps[[SC.MEB_index]][[10]], reference_vector)
se_SC.MEB <- sd(SC.MEB_acc) / sqrt(length(SC.MEB_acc))

SC.MEB_results <- plot_points_paper(coords, SC.MEB_partition, bnd_scaled) +
  geom_boundary(bnd_scaled) +
  scale_color_manual(values = my_palette1) +
  ggtitle(bquote(bold("SC-MEB, ") ~ italic("ARI = ") ~ .(round(max(SC.MEB_acc), 3)) ~ " (" ~ .(round(se_SC.MEB, 3)) ~ ")")) +
  theme(plot.title = element_text(size = 20))
SC.MEB_results

SC.MEB_results_with_label <- ggdraw() +
  draw_plot(SC.MEB_results) +
  draw_label("(iii)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
SC.MEB_results_with_label

ggsave("./New_Simulations/SR_Extra/SR_sim_SC.MEB.png", 
       plot = SC.MEB_results_with_label, 
       width = 5, height = 5, units = "in", dpi = 300, bg = "transparent")

#-------------------------------------------------------------------------------

##### DR-RC results #####
load("~/SR_sim_DR.SC_30reps.RData")

DR.SC_acc <- c()
for (i in 1:30){
  DR.SC_acc[i] <- adjustedRandIndex(SR_sim_data$cluster, SR_sim_DR.SC_reps[[i]][["Objdrsc"]][[1]][["cluster"]])
}

plot(DR.SC_acc, type = "l")
summary(DR.SC_acc)

DR.SC_index <- which.max(DR.SC_acc)

DR.SC_partition <- reorder_based_on_reference(SR_sim_DR.SC_reps[[DR.SC_index]][["Objdrsc"]][[1]][["cluster"]], reference_vector)
se_DR.SC <- sd(DR.SC_acc) / sqrt(length(DR.SC_acc))

DR.SC_results <- plot_points_paper(coords, DR.SC_partition, bnd_scaled) +
  geom_boundary(bnd_scaled) +
  scale_color_manual(values = my_palette1) +
  ggtitle(bquote(bold("DR-SC, ") ~ italic("ARI = ") ~ .(round(max(DR.SC_acc), 3)) ~ " (" ~ .(round(se_DR.SC, 3)) ~ ")")) +
  theme(plot.title = element_text(size = 20))
DR.SC_results

DR.SC_results_with_label <- ggdraw() +
  draw_plot(DR.SC_results) +
  draw_label("(iv)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
DR.SC_results_with_label

ggsave("./New_Simulations/SR_Extra/SR_sim_DR.SC.png", 
       plot = DR.SC_results_with_label, 
       width = 5, height = 5, units = "in", dpi = 300, bg = "transparent")

#-------------------------------------------------------------------------------

##### k-means results #####
load("/SR_sim_kmeans_30reps.RData")

kmeans_acc <- c()
for (i in 1:30){
  kmeans_acc[i] <- adjustedRandIndex(SR_sim_data$cluster, SR_sim_kmeans_reps[[i]])
}

plot(kmeans_acc, type = "l")
summary(kmeans_acc)

kmeans_index <- which.max(kmeans_acc)

kmeans_partition <- reorder_based_on_reference(SR_sim_kmeans_reps[[kmeans_index]], reference_vector)
se_kmeans <- sd(kmeans_acc) / sqrt(length(kmeans_acc))

kmeans_results <- plot_points_paper(coords, kmeans_partition, bnd_scaled) +
  geom_boundary(bnd_scaled) +
  scale_color_manual(values = my_palette1) +
  ggtitle(bquote(bold("k-means, ") ~ italic("ARI = ") ~ .(round(max(kmeans_acc), 3)) ~ " (" ~ .(round(se_kmeans, 3)) ~ ")")) +
  theme(plot.title = element_text(size = 20))
kmeans_results

kmeans_results_with_label <- ggdraw() +
  draw_plot(kmeans_results) +
  draw_label("(v)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
kmeans_results_with_label

ggsave("./New_Simulations/SR_Extra/SR_sim_kmeans.png", 
       plot = kmeans_results_with_label, 
       width = 5, height = 5, units = "in", dpi = 300, bg = "transparent")

#-------------------------------------------------------------------------------
##### Load GraphST data #####
GraphST_sim_results <- readr::read_csv("./New_Simulations/SR_Extra/graphst_swiss_roll_multiseed_results.csv")

# Extract ARIs from all GraphST columns
graphst_partition_matrix <- GraphST_sim_results %>%
  select(starts_with("refined_mclust_seed_")) %>%
  as.matrix()

graphst_ari <- apply(
  graphst_partition_matrix, 
  2, 
  function(col) adjustedRandIndex(GraphST_sim_results$original_cluster, col)
)

graphst_best_index <- which.max(graphst_ari)
graphst_best_ari <- graphst_ari[graphst_best_index]
graphst_se <- sd(graphst_ari) / sqrt(length(graphst_ari))

# Best partition
graphst_best_partition <- reorder_based_on_reference(
  graphst_partition_matrix[, graphst_best_index],
  reference_vector
)

# Plot
graphst_plot <- plot_points_paper(
  coords,
  graphst_best_partition,
  bnd_scaled
) +
  geom_boundary(bnd_scaled) +
  scale_color_manual(values = my_palette1) +
  ggtitle(bquote(bold("GraphST, ") ~ italic("ARI = ") ~ .(round(graphst_best_ari, 3)) ~ 
                   " (" ~ .(round(graphst_se, 3)) ~ ")")) +
  theme(plot.title = element_text(size = 20))
graphst_plot

graphst_plot_with_label <- ggdraw() +
  draw_plot(graphst_plot) +
  draw_label("(vi)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)

# Save the figure
ggsave("./New_Simulations/SR_Extra/SR_sim_GraphST.png",
       plot = graphst_plot_with_label,
       width = 5, height = 5, units = "in", dpi = 300, bg = "transparent")

#-------------------------------------------------------------------------------
##### Load SEDR data #####
SEDR_sim_results <- readr::read_csv("./New_Simulations/SR_Extra/sedr_swiss_roll_multiseed_k3_results.csv")

# Extract ARIs from all SEDR columns
sedr_partition_matrix <- SEDR_sim_results %>%
  select(starts_with("mclust_seed_")) %>%
  as.matrix()

sedr_ari <- apply(
  sedr_partition_matrix, 
  2, 
  function(col) adjustedRandIndex(SEDR_sim_results$original_cluster, col)
)

sedr_best_index <- which.max(sedr_ari)
sedr_best_ari <- sedr_ari[sedr_best_index]
sedr_se <- sd(sedr_ari) / sqrt(length(sedr_ari))

# Best partition
sedr_best_partition <- reorder_based_on_reference(
  sedr_partition_matrix[, sedr_best_index],
  reference_vector
)

# Plot
sedr_plot <- plot_points_paper(
  coords,
  sedr_best_partition,
  bnd_scaled
) +
  geom_boundary(bnd_scaled) +
  scale_color_manual(values = my_palette1) +
  ggtitle(bquote(bold("SEDR, ") ~ italic("ARI = ") ~ .(round(sedr_best_ari, 3)) ~ 
                   " (" ~ .(round(sedr_se, 3)) ~ ")")) +
  theme(plot.title = element_text(size = 20))
sedr_plot

sedr_plot_with_label <- ggdraw() +
  draw_plot(sedr_plot) +
  draw_label("(vii)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)

# Save the figure
ggsave("./New_Simulations/SR_Extra/SR_sim_SEDR.png",
       plot = sedr_plot_with_label,
       width = 5, height = 5, units = "in", dpi = 300, bg = "transparent")

#-------------------------------------------------------------------------------
##### Load STAGATE data #####
STAGATE_sim_results <- readr::read_csv("./New_Simulations/SR_Extra/stagate_simdata_multiseed_3clusters_results.csv")

# Extract ARIs from all STAGATE columns
stagate_partition_matrix <- STAGATE_sim_results %>%
  select(starts_with("stagate_seed_")) %>%
  as.matrix()

stagate_ari <- apply(
  stagate_partition_matrix, 
  2, 
  function(col) adjustedRandIndex(STAGATE_sim_results$cluster, col)
)

stagate_best_index <- which.max(stagate_ari)
stagate_best_ari <- stagate_ari[stagate_best_index]
stagate_se <- sd(stagate_ari) / sqrt(length(stagate_ari))

# Best partition
stagate_best_partition <- reorder_based_on_reference(
  stagate_partition_matrix[, stagate_best_index],
  reference_vector
)

# Plot
stagate_plot <- plot_points_paper(
  coords,
  stagate_best_partition,
  bnd_scaled
) +
  geom_boundary(bnd_scaled) +
  scale_color_manual(values = my_palette1) +
  ggtitle(bquote(bold("STAGATE, ") ~ italic("ARI = ") ~ .(round(stagate_best_ari, 3)) ~ 
                   " (" ~ .(round(stagate_se, 3)) ~ ")")) +
  theme(plot.title = element_text(size = 20))
stagate_plot

stagate_plot_with_label <- ggdraw() +
  draw_plot(stagate_plot) +
  draw_label("(viii)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)

# Save the figure
ggsave("./New_Simulations/SR_Extra/SR_sim_STAGATE.png",
       plot = stagate_plot_with_label,
       width = 5, height = 5, units = "in", dpi = 300, bg = "transparent")

#-------------------------------------------------------------------------------
##### Load SpaGCN data #####
SpaGCN_sim_results <- readr::read_csv("./New_Simulations/SR_Extra/SpaGCN_simulated_clustering_results_multi_seed.csv")

# Extract ARIs from all SpaGCN columns
spagcn_partition_matrix <- SpaGCN_sim_results %>%
  select(starts_with("refined_pred_seed_")) %>%
  as.matrix()

spagcn_ari <- apply(
  spagcn_partition_matrix, 
  2, 
  function(col) adjustedRandIndex(SpaGCN_sim_results$ground_truth_cluster, col)
)

spagcn_best_index <- which.max(spagcn_ari)
spagcn_best_ari <- spagcn_ari[spagcn_best_index]
spagcn_se <- sd(spagcn_ari) / sqrt(length(spagcn_ari))

# Best partition
spagcn_best_partition <- reorder_based_on_reference(
  spagcn_partition_matrix[, spagcn_best_index]+1,
  reference_vector
)

# Plot
spagcn_plot <- plot_points_paper(
  coords,
  spagcn_best_partition,
  bnd_scaled
) +
  geom_boundary(bnd_scaled) +
  scale_color_manual(values = my_palette1) +
  ggtitle(bquote(bold("SpaGCN, ") ~ italic("ARI = ") ~ .(round(spagcn_best_ari, 3)) ~ 
                   " (" ~ .(round(spagcn_se, 3)) ~ ")")) +
  theme(plot.title = element_text(size = 20))
spagcn_plot

spagcn_plot_with_label <- ggdraw() +
  draw_plot(spagcn_plot) +
  draw_label("(ix)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)

# Save the figure
ggsave("./New_Simulations/SR_Extra/SR_sim_SpaGCN.png",
       plot = spagcn_plot_with_label,
       width = 5, height = 5, units = "in", dpi = 300, bg = "transparent")


#-------------------------------------------------------------------------------
##### Create box-plot for accuracy results for all methods#####
# Combine the vectors into a data frame
data <- data.frame(
  Accuracy = c(DPM_accur, BayesSpace_acc, SC.MEB_acc, DR.SC_acc, kmeans_acc, 
               graphst_ari, sedr_ari, stagate_ari, spagcn_ari),
  Method = factor(rep(c("DP-RST", "Bayes\nSpace", "DR.SC", "SC.MEB", "k-means", 
                        "GraphST", "SEDR", "STAGATE", "SpaGCN"), each = 30), 
                  levels = c("DP-RST", "Bayes\nSpace", "DR.SC", "SC.MEB", "k-means", 
                             "GraphST", "SEDR", "STAGATE", "SpaGCN"))
)

# Create the box plot
boxplot <- ggplot(data, aes(x = Method, y = Accuracy, fill = Method)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = c("DP-RST" = "#54C4E3FF", "Bayes\nSpace" = "#D84860FF", 
                               "DR.SC" = "#903030FF", "SC.MEB" = "#FFFFFFFF", 
                               "k-means" = "#A8A8A8FF", "GraphST" = "#FFFFFFFF",
                               "SEDR" = "#F0A8C0FF", "STAGATE" = "#F0A8C0FF", 
                               "SpaGCN" = "#E0B2CDFF")) +
  labs(title = "Compare ARI",
       x = NULL,  # Remove the x-axis label
       y = "ARI") +
  theme(legend.position = "none",
        plot.title = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.margin = unit(c(2.5, 2.5, 2.5, 2.5), "cm"))  # Add margins around the plot
boxplot

boxplot_with_label <- ggdraw() +
  draw_plot(boxplot) +
  draw_label("(B)", x = 0.1, y = 0.98, size = 20, hjust = 0, vjust = 1)
boxplot_with_label

ggsave("./New_Simulations/SR_Extra/SR_sim_boxplot.png", 
       plot = boxplot_with_label, 
       width = 12, height = 6, units = "in", dpi = 300, bg = "transparent")
