##### LOAD NECESSARY LIBRARIES #####

library(DP.RST)
library(mclust)
library(ggplot2)
library(readr)
library(dplyr)

#-------------------------------------------------------------------------------
# Set working directory
setwd("~/Desktop/DP-RST-workload/Datasets Results/Gut")
getwd()

# Set seed
set.seed(7303)

# For debugging purposes
options(error=recover)
#-------------------------------------------------------------------------------
##### LOAD PRE-PROCESSED INTESTINE DATA #####

# Load pre-processed data
gut_df_wt_muscle <- readRDS("./gut_df_wt_muscle.rds")
load("swiss_roll_wt_muscle_boundary.RData")

# Extract spatial coordinates and PCs
Y_sample = gut_df_wt_muscle[, 1:3]
loc = gut_df_wt_muscle[, c("x", "y")]
bnd = boundary

# Extract hystological labels
z = gut_df_wt_muscle$z

# Standardize coordinates, Y values
coords <- apply(loc, 2, scale)
Y_std <- apply(Y_sample, 2, scale)

# Standardize the boundary to match the scaled coordinates
bnd_scaled = list()
bnd_scaled$x <- (bnd$x - mean(loc[,1]))/sd(loc[,1])
bnd_scaled$y <- (bnd$y - mean(loc[,2]))/sd(loc[,2])

# Create dataframe for plotting
df_subset <- data.frame(loc, z)

true_hystology_labels <- ggplot() +
  geom_boundary(bnd_scaled) +
  geom_point(aes(x = x, y = y, col = z), data = data.frame(coords, z)) +
  labs(x = 'Coordinate X', y = 'Coordinate Y', colour = "True label") +
  ggtitle('Ground Truth Labels
          (manually groupped)') +
  theme(legend.position = "bottom")
true_hystology_labels

#-------------------------------------------------------------------------------
##### FUNCTION TO REORDER CLUSTER LABELS BASED ON THE REFERENCE VECTOR #####

# Create a reference vector with values from 1 to 9
reference_vector <- c(1:5)

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

#-------------------------------------------------------------------------------
##### PLOTTING FUNCTION WITHOUT BOUNDARY #####

groups_plot_wB <- function(coords, group_assign, point_size = 2) {
  
  # Create a dataframe of coordinates and groups' assignment
  group_data <- as.data.frame(cbind(coords, group_assign))
  group_data$group_assign <- as.factor(group_data$group_assign)
  
  # Make a plot
  group_pred <- ggplot() + 
    geom_point(aes(x = group_data[,1], y = group_data[,2], colour = as.factor(group_data[,3])), 
               data = group_data, size = point_size) +
    # geom_point(aes(x = group_data[,1], y = group_data[,2], colour = group_assign), 
    #            data = group_data, size = point_size) +
    # scale_colour_manual(values = colors, name = 'Groups', guide = guide_legend(ncol = 2)) +
    # scale_colour_manual(values = getPalette(colourCount), name = 'groups', guide = guide_legend(ncol = 2)) +
    labs(x = 'Coordinate X', y = 'Coordinate Y') + 
    ggtitle('Groups of Observations') 
  
  return(group_pred)
}
#-------------------------------------------------------------------------------
##### LOAD AND PROCESS DP-RST OUTPUT #####

load("./Results/Gut_DP.RST_FromNewBastPT_p3_Version2_OutputOnly.RData")
output = Gut_DP.RST_FromNewBastPT_Version2

# Get the best partition
mode_based_partition = partition(DP.RST_output = output, method = "mode_based", 
                                 batch_size = 100)

#### Plot the crude spatial and refined partitions #####
# Determine the number of unique groups in the crude spatial partition output
k_groups_out <- length(unique(mode_based_partition$groups_partition)) # 20

# Calculate the accuracy of the crude spatial partition by comparing it with the true labels (df_subset$z)
accur_groups = adjustedRandIndex(z, mode_based_partition$groups_partition)
accur_groups # 0.14534

### Plot the crude spatial partition
groups_plotted <- groups_plot(coords, mode_based_partition$groups_partition, bnd_scaled) +
  theme(legend.position = "none") +
  ggtitle(paste("Estimated Groups, k = ", k_groups_out, ". Spatial Clustering.
                ARI (gr.) = ",
                round(accur_groups, 3), sep = ""))
groups_plotted

# Plot the teams of observations (clusters in the refined partition)
k_teams_out <- length(unique(mode_based_partition$teams_partition)) # 5

#### Calculate the accuracy of the labeled observations #####
# Binary matrix for groups and teams
X <- table(sequence(length(mode_based_partition$groups_partition)), 
           mode_based_partition$groups_partition)
Z <- table(sequence(length(mode_based_partition$teams_partition)), 
           mode_based_partition$teams_partition)

# Get the membership of observations in each team
obs_in_teams <- X %*% Z
# Compute W such that w_ij = 1 if observation i and observation j share same team
obs_in_teams_vec <- obs_in_teams %*% sort(unique(mode_based_partition$teams_partition))

accur_teams = adjustedRandIndex(df_subset$z, obs_in_teams_vec) #
accur_teams # 0.4348815

#-------------------------------------------------------------------------------
##### PLOT DP-RST RESULTS ###

# Define a custom color palette for plotting
my_palette1 <- c("#FCDD23FF", "#AD98F1FF", "#D856A7FF", "#F8B100FF", "#48439BFF")

# Reorder the partition based on a reference vector
DPM_partition <- reorder_based_on_reference(obs_in_teams_vec, reference_vector)

# Plot the refined partition (DP-RST partition)
DPM.RST_plot <- groups_plot(coords, DPM_partition, bnd_scaled, point_size = 3) +
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
  ggtitle(bquote(bold("DP-RST") * ", c = " * .(k_teams_out) * ", " * bold("ARI = ") * bold(.(format(round(accur_teams, 3), nsmall = 3)))))

DPM.RST_plot

# Save the refined partition plot as a PNG image with a transparent background (part of the Figure 4)
ggsave(filename = "~/Desktop/DP-RST-workload/Figures/figure3/DPM-RST_plot.png",
       plot = DPM.RST_plot, width = 8, height = 8, bg = "transparent")

#-------------------------------------------------------------------------------
##### DENSITY PLOTS FOR THE DP-RST RESULTS #####

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
hist_DPM <- gridExtra::grid.arrange(plot_PC1, plot_PC2, plot_PC3, ncol = 3)
hist_DPM

# # Save the arranged density plots as a PNG image (part of the Figure 4)
# ggsave("./Figures/figure4/hist_DP-RST.png",
#        plot = hist_DPM, width = 7.97, height = 5.35,
#        units = "in", dpi = 300, bg = "transparent")

#-------------------------------------------------------------------------------
##### HISTOGRAM OF REFINED CLUSTERS (MCMC ITERATIONS) #####

# Extract the frequency of team counts across all MCMC iterations
# output$j_teams_out contains the number of teams for each iteration
teams_freq <- as.data.frame(table(output$j_teams_out))

# Plot the histogram
cluster_count_hist <- ggplot(teams_freq, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", color = "#C497A7", fill = "#F8E6EC", width = 0.8) +
  labs(title = "Posterior Distribution of Refined Clusters",
       # subtitle = "Number of Teams over MCMC Iterations",
       x = "Number of Refined Clusters", 
       y = "Frequency") +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank(), # Remove vertical grid lines for cleaner look
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# Display the plot
print(cluster_count_hist)

# Save the histogram as a PNG image (part of Figure 4)
ggsave("~/Desktop/DP-RST-workload/Figures/figure3/hist_MCMC_clusters.png",
       plot = cluster_count_hist, width = 8, height = 6,
       units = "in", dpi = 300, bg = "transparent")
#-------------------------------------------------------------------------------
##### CREATE PLOTS FOR THE FIGURE 1 #####

##### Plot the crude spatial partition for Figure 1 #####
DPM_teams_plot <- groups_plot(coords, DPM_partition, bnd_scaled, point_size = 3) +
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
my_palette2 <- c(
  "#FCDD23FF", "#F8B100FF", "#CA697CFF", "#AB74CFFF", "#48439BFF",
  "#D9A453FF", "#E3A8A3FF", "#6A47D8FF", "#D856A7FF", "#3B83BDFF",
  "#FF7F50FF", "#9ACD32FF", "#D2691EFF", "#1E90FFFF", "#228B22FF", 
  "#40E0D0FF", "#DC143CFF", "#87CEFAFF", "#6A5ACDFF", "#C71585FF" )

DPM_groups_plot <- groups_plot(coords, mode_based_partition$groups_partition, 
                               bnd_scaled, point_size = 3) +
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

ggsave(filename = "./Figures/figure4/DP-RST_refined_Figure1.png",
       plot = DPM_teams_plot, width = 8, height = 8, bg = "transparent")

ggsave(filename = "./Figures/figure4/DP-RST_crude_Figure1.png",
       plot = DPM_groups_plot, width = 8, height = 8, bg = "transparent")

#-------------------------------------------------------------------------------
##### LOAD BAYES-SPACE RESULTS #####

BayesSpace_Gut_clustering_results <- readRDS("./Results/BayesSpace_Gut_clustering_results.rds")

BS.df <- data.frame(x = BayesSpace_Gut_clustering_results$imagerow, 
                    y = BayesSpace_Gut_clustering_results$imagecol, 
                    BS_z = BayesSpace_Gut_clustering_results$cluster_q5_pc3)

BS_merged_df <- merge(BS.df, df_subset,
                      by = c("x", "y"),
                      all = TRUE)

BS_accur = adjustedRandIndex(BS_merged_df$z, BS_merged_df$BS_z)
BS_accur # 0.1280336

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
  ggtitle(bquote(bold("BayesSpace") * ", c = 5, " * bold("ARI = ") * bold(.(format(round(BS_accur, 3), nsmall = 3)))))

BS_plot

# Save the BayesSpace partition plot as a PNG image (part of Figure 4)
ggsave(filename = "~/Desktop/DP-RST-workload/Figures/figure3/BS_plot.png",
       plot = BS_plot, width = 8, height = 8, bg = "transparent")

#-------------------------------------------------------------------------------
##### LOAD SC-MEB RESULTS #####

SC_MEB_Gut_clustering_results <- readRDS("./Results/SC-MEB_Gut_clustering_results.rds")

SC_MEB.df <- data.frame(x = SC_MEB_Gut_clustering_results$x, 
                    y = SC_MEB_Gut_clustering_results$y, 
                    SC_z = SC_MEB_Gut_clustering_results$pc3_k5)

SC_MEB_merged_df <- merge(SC_MEB.df, df_subset,
                      by = c("x", "y"),
                      all = TRUE)

SC.MEB_accur = adjustedRandIndex(SC_MEB_merged_df$z, SC_MEB_merged_df$SC_z)
SC.MEB_accur # 0.1124159

# Reorder the partition based on a reference vector
SC.MEB_partition <- reorder_based_on_reference(SC_MEB_merged_df$SC_z, reference_vector)

SC.MEB_plot <- groups_plot_wB(SC_MEB_merged_df[, 1:2], SC.MEB_partition, point_size = 3) +
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
  ggtitle(bquote(bold("SC-MEB") * ", c = 5, " * bold("ARI = ") * bold(.(format(round(SC.MEB_accur, 3), nsmall = 3)))))

SC.MEB_plot

# Save the SC-MEB partition plot as a PNG image (part of Figure 4)
ggsave(filename = "./Figures/figure4/SC.MEB_plot.png",
       plot = SC.MEB_plot, width = 8, height = 8, bg = "transparent")

#-------------------------------------------------------------------------------
##### LOAD DR-SC RESULTS #####

DR_SC_Gut_clustering_results <- readRDS("./Results/DR-SC_Gut_clustering_results.rds")

DR_SC.df <- data.frame(x = DR_SC_Gut_clustering_results$imagerow, 
                        y = DR_SC_Gut_clustering_results$imagecol, 
                        DR_z = DR_SC_Gut_clustering_results$cluster_q5_pc3)

DR_SC_merged_df <- merge(DR_SC.df, df_subset,
                          by = c("x", "y"),
                          all = TRUE)

DR_SC_accur = adjustedRandIndex(DR_SC_merged_df$z, DR_SC_merged_df$DR_z)
DR_SC_accur # 0.1530579

# Reorder the partition based on a reference vector
DR.SC_partition <- reorder_based_on_reference(DR_SC_merged_df$DR_z, reference_vector)

DR.SC_plot <- groups_plot_wB(DR_SC_merged_df[, 1:2], DR.SC_partition, point_size = 3) +
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
  ggtitle(bquote(bold("DR-SC") * ", c = 5, " * bold("ARI = ") * bold(.(format(round(DR_SC_accur, 3), nsmall = 3)))))

DR.SC_plot

# Save the DR-SC partition plot as a PNG image (part of Figure 4)
ggsave(filename = "./Figures/figure4/DR.SC_plot.png",
       plot = DR.SC_plot, width = 8, height = 8, bg = "transparent")


#-------------------------------------------------------------------------------
##### LOAD kmeans RESULTS #####

kmeans_Gut_clustering_results <- readRDS("./Results/kmeans_Gut_clustering_results.rds")

kmeans.df <- data.frame(x = kmeans_Gut_clustering_results$x, 
                       y = kmeans_Gut_clustering_results$y, 
                       kmeans_z = kmeans_Gut_clustering_results$pc3_k5)

kmeans_merged_df <- merge(kmeans.df, df_subset,
                         by = c("x", "y"),
                         all = TRUE)

kmeans_accur = adjustedRandIndex(kmeans_merged_df$z, kmeans_merged_df$kmeans_z)
kmeans_accur # 0.02973941

# Reorder the partition based on a reference vector
kmeans_partition <- reorder_based_on_reference(kmeans_merged_df$kmeans_z, reference_vector)

kmeans_plot <- groups_plot_wB(kmeans_merged_df[, 1:2], kmeans_partition, point_size = 3) +
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
  ggtitle(bquote(bold("k-means") * ", c = 5, " * bold("ARI = ") * bold(.(format(round(kmeans_accur, 3), nsmall = 3)))))

kmeans_plot

# Save the k-means partition plot as a PNG image (part of Figure 4)
ggsave(filename = "./Figures/figure4/kmeans_plot.png",
       plot = kmeans_plot, width = 8, height = 8, bg = "transparent")

#-------------------------------------------------------------------------------
##### LOAD BASS RESULTS #####

BASS_Gut_clustering_results <- readRDS("./Results/BASS_Intestine_3PCs.rds")

bass.df <- data.frame(x = BASS_Gut_clustering_results$x, 
                      y = BASS_Gut_clustering_results$y, 
                      bass_z = BASS_Gut_clustering_results$z)

bass_merged_df <- merge(bass.df, df_subset,
                        by = c("x", "y"),
                        all = TRUE)

bass_accur = adjustedRandIndex(bass_merged_df$z, bass_merged_df$bass_z)
bass_accur # 0.2790147

# Reorder the partition based on a reference vector
bass_partition <- reorder_based_on_reference(bass_merged_df$bass_z, reference_vector)

bass_plot <- groups_plot_wB(bass_merged_df[, 1:2], bass_partition, point_size = 3) +
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
  ggtitle(bquote(bold("BASS") * ", c = 5, " * bold("ARI = ") * bold(.(format(round(bass_accur, 3), nsmall = 3)))))

bass_plot

# Save the BASS partition plot as a PNG image (part of Figure 4)
ggsave(filename = "~/Desktop/DP-RST-workload/Figures/figure3/BASS_plot.png",
       plot = bass_plot, width = 8, height = 8, bg = "transparent")

#-------------------------------------------------------------------------------
##### LOAD SEDR RESULTS #####

sedr_gut_results <- read_csv("./Competing_Methods/SEDR_runs/res_realdata_csv/sedr_gut_reduced_pca3_results.csv")

sedr.df <- data.frame(y = sedr_gut_results$x, 
                      x = sedr_gut_results$y, 
                      sedr_z = sedr_gut_results$mclust_pca3_5)

sedr_merged_df <- merge(sedr.df, df_subset,
                          by = c("x", "y"),
                          all = TRUE)

sedr_accur = adjustedRandIndex(sedr_merged_df$z, sedr_merged_df$sedr_z)
sedr_accur # 0.0481976

# Reorder the partition based on a reference vector
sedr_partition <- reorder_based_on_reference(sedr_merged_df$sedr_z, reference_vector)

sedr_plot <- groups_plot_wB(sedr_merged_df[, 1:2], sedr_partition, point_size = 3) +
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
  ggtitle(bquote(bold("SEDR") * ", c = 5, " * bold("ARI = ") * bold(.(format(round(sedr_accur, 3), nsmall = 3)))))

sedr_plot

# Save the k-means partition plot as a PNG image (part of Figure 4)
ggsave(filename = "./Figures/figure4/sedr_plot.png",
       plot = sedr_plot, width = 8, height = 8, bg = "transparent")


#-------------------------------------------------------------------------------
##### LOAD GRAPHST RESULTS #####

graphst_gut_results <- read_csv("./Competing_Methods/GraphST_runs/res_realdata_csv/graphst_gut_reduced_results.csv")

graphst.df <- data.frame(y = graphst_gut_results$x, 
                      x = graphst_gut_results$y, 
                      graphst_z = graphst_gut_results$mclust_5)

graphst_merged_df <- merge(graphst.df, df_subset,
                        by = c("x", "y"),
                        all = TRUE)

graphst_accur = adjustedRandIndex(graphst_merged_df$z, graphst_merged_df$graphst_z)
graphst_accur # 0.1419926

# Reorder the partition based on a reference vector
graphst_partition <- reorder_based_on_reference(graphst_merged_df$graphst_z, reference_vector)

graphst_plot <- groups_plot_wB(graphst_merged_df[, 1:2], graphst_partition, point_size = 3) +
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
  ggtitle(bquote(bold("GraphST") * ", c = 5, " * bold("ARI = ") * bold(.(format(round(graphst_accur, 3), nsmall = 3)))))

graphst_plot

# Save the k-means partition plot as a PNG image (part of Figure 4)
ggsave(filename = "./Figures/figure4/graphst_plot.png",
       plot = graphst_plot, width = 8, height = 8, bg = "transparent")


#-------------------------------------------------------------------------------
##### LOAD STAGATE RESULTS #####

stagate_gut_results <- read_csv("./Competing_Methods/STAGATE_runs/res_realdata_csv/stagate_gut_reduced_results.csv")

stagate.df <- data.frame(y = stagate_gut_results$x, 
                         x = stagate_gut_results$y, 
                         stagate_z = stagate_gut_results$mclust_5)

stagate_merged_df <- merge(stagate.df, df_subset,
                           by = c("x", "y"),
                           all = TRUE)

stagate_accur = adjustedRandIndex(stagate_merged_df$z, stagate_merged_df$stagate_z)
stagate_accur # 0.2669045

# Reorder the partition based on a reference vector
stagate_partition <- reorder_based_on_reference(stagate_merged_df$stagate_z, reference_vector)

stagate_plot <- groups_plot_wB(stagate_merged_df[, 1:2], stagate_partition, point_size = 3) +
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
  ggtitle(bquote(bold("STAGATE") * ", c = 5, " * bold("ARI = ") * bold(.(format(round(stagate_accur, 3), nsmall = 3)))))

stagate_plot

# Save the k-means partition plot as a PNG image (part of Figure 4)
ggsave(filename = "./Figures/figure4/stagate_plot.png",
       plot = stagate_plot, width = 8, height = 8, bg = "transparent")


#-------------------------------------------------------------------------------
##### LOAD SpaGCN RESULTS #####

spagcn_gut_results <- read_csv("Competing_Methods/SpaGCN_runs/res_realdata_csv/all_clustering_results_Gut_reduced.csv")

spagcn.df <- data.frame(spots = spagcn_gut_results$...1, 
                        x = spagcn_gut_results$array_row, 
                        y = spagcn_gut_results$array_col, 
                        spagcn_z_with = (spagcn_gut_results$refined_pred_5clusters_3pcs_with_histology+1),
                        spagcn_z_without = (spagcn_gut_results$refined_pred_5clusters_3pcs_without_histology+1))

df_subset$spots <- rownames(df_subset)
spagcn_merged_df <- merge(df_subset, spagcn.df, by = "spots")

spagcn_accur_with = adjustedRandIndex(spagcn_merged_df$spagcn_z_with, spagcn_merged_df$z)
spagcn_accur_with # 0.1304065

# Reorder the partition based on a reference vector
spagcn_partition_with <- reorder_based_on_reference(spagcn_merged_df$spagcn_z_with, reference_vector)

spagcn_plot_with <- groups_plot_wB(spagcn_merged_df[, c("x.x", "y.x")], spagcn_partition_with, point_size = 3) +
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
  ggtitle(bquote(bold("SpaGCN (w/)") * ", c = 5, " * bold("ARI = ") * bold(.(format(round(spagcn_accur_with, 3), nsmall = 3)))))

spagcn_plot_with

# Save the k-means partition plot as a PNG image (part of Figure 4)
ggsave(filename = "./Figures/figure4/spagcn_with_plot.png",
       plot = spagcn_plot_with, width = 8, height = 8, bg = "transparent")


spagcn_accur_without = adjustedRandIndex(spagcn_merged_df$spagcn_z_without, spagcn_merged_df$z)
spagcn_accur_without # 0.1990036

# Reorder the partition based on a reference vector
spagcn_partition_without <- reorder_based_on_reference(spagcn_merged_df$spagcn_z_without, reference_vector)

spagcn_plot_without <- groups_plot_wB(spagcn_merged_df[, c("x.x", "y.x")], spagcn_partition_without, point_size = 3) +
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
  ggtitle(bquote(bold("SpaGCN (w/o)") * ", c = 5, " * bold("ARI = ") * bold(.(format(round(spagcn_accur_without, 3), nsmall = 3)))))

spagcn_plot_without

# Save the k-means partition plot as a PNG image (part of Figure 4)
ggsave(filename = "./Figures/figure4/spagcn_without_plot.png",
       plot = spagcn_plot_without, width = 8, height = 8, bg = "transparent")


#-------------------------------------------------------------------------------
##### LOAD stLearn RESULTS #####

stlearn_gut_results_with <- read_csv("./Competing_Methods/stLearn_runs/res_realdata_csv/Gut_reduced_clustering_3PCs_5clusters_with_image.csv")
stlearn_gut_results_without <- read_csv("./Competing_Methods/stLearn_runs/res_realdata_csv/Gut_reduced_clustering_3PCs_5clusters_no_image.csv")

stlearn.df <- data.frame(y = stlearn_gut_results_with$x, 
                         x = stlearn_gut_results_with$y, 
                         stlearn_z_with = (stlearn_gut_results_with$cluster_5+1),
                         stlearn_z_without = (stlearn_gut_results_without$cluster_5+1))

stlearn_merged_df <- merge(stlearn.df, df_subset,
                           by = c("x", "y"),
                           all = TRUE)

stlearn_accur_with = adjustedRandIndex(stlearn_merged_df$stlearn_z_with, stlearn_merged_df$z)
stlearn_accur_with # 0.2753501

# Reorder the partition based on a reference vector
stlearn_partition_with <- reorder_based_on_reference(stlearn_merged_df$stlearn_z_with, reference_vector)

stlearn_plot_with <- groups_plot_wB(stlearn_merged_df[, c("x", "y")], stlearn_partition_with, point_size = 3) +
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
  ggtitle(bquote(bold("stLearn (w/)") * ", c = 5, " * bold("ARI = ") * bold(.(format(round(stlearn_accur_with, 3), nsmall = 3)))))

stlearn_plot_with

# Save the k-means partition plot as a PNG image (part of Figure 4)
ggsave(filename = "./Figures/figure4/stlearn_with_plot.png",
       plot = stlearn_plot_with, width = 8, height = 8, bg = "transparent")


stlearn_accur_without = adjustedRandIndex(stlearn_merged_df$stlearn_z_without, stlearn_merged_df$z)
stlearn_accur_without # 0.1507263

# Reorder the partition based on a reference vector
stlearn_partition_without <- reorder_based_on_reference(stlearn_merged_df$stlearn_z_without, reference_vector)

stlearn_plot_without <- groups_plot_wB(stlearn_merged_df[, c("x", "y")], stlearn_partition_without, point_size = 3) +
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
  ggtitle(bquote(bold("stLearn (w/o)") * ", c = 5, " * bold("ARI = ") * bold(.(format(round(stlearn_accur_without, 3), nsmall = 3)))))

stlearn_plot_without

# Save the k-means partition plot as a PNG image (part of Figure 4)
ggsave(filename = "./Figures/figure4/stlearn_without_plot.png",
       plot = stlearn_plot_without, width = 8, height = 8, bg = "transparent")

#-------------------------------------------------------------------------------
##### LOAD SpaMask RESULTS #####

SpaMask_gut_results <- read_csv("~/Desktop/DP-RST-workload/Datasets Results/Gut/Results/SpaMask_Gut_3PCs.csv")

SpaMask.df <- data.frame(y = SpaMask_gut_results$x, 
                         x = SpaMask_gut_results$y, 
                         SpaMask_z = SpaMask_gut_results$kmeans)
# Fix SpaMask.df so columns follow (x, y) convention
SpaMask.df.fixed <- SpaMask.df %>%
  dplyr::rename(
    x = y,
    y = x
  )

SpaMask_merged_df <- merge(
  SpaMask.df.fixed,
  df_subset,
  by = c("x", "y"),
  all = TRUE
)

SpaMask_accur = adjustedRandIndex(SpaMask_merged_df$z, SpaMask_merged_df$SpaMask_z)
SpaMask_accur # 0.2000117

# Reorder the partition based on a reference vector
SpaMask_partition <- reorder_based_on_reference(SpaMask_merged_df$SpaMask_z + 1, reference_vector)

SpaMask_plot <- groups_plot_wB(SpaMask_merged_df[, 1:2], SpaMask_partition, point_size = 3) +
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
  ggtitle(bquote(bold("SpaMask") * ", c = 5, " * bold("ARI = ") * bold(.(format(round(SpaMask_accur, 3), nsmall = 3)))))

SpaMask_plot

# Save the k-means partition plot as a PNG image (part of Figure 4)
ggsave(filename = "~/Desktop/DP-RST-workload/Figures/figure3/SpaMask_plot.png",
       plot = SpaMask_plot, width = 8, height = 8, bg = "transparent")
