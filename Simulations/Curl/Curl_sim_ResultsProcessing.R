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
Curl_sim_data <- readr::read_csv("./Simulations/Curl/Curl_sim_data.csv")

### Load Boundary ###
# Unzip the folder
zip_file <- "./Simulations/Curl/Curl_contours_csv.zip"
unzip(zip_file, exdir = "contours_csv")

# List all CSV files
csv_files <- list.files("contours_csv", pattern = "\\.csv$", full.names = TRUE)

# Read and combine all contours into one data frame
all_contours <- lapply(seq_along(csv_files), function(i) {
  df <- read.csv(csv_files[i])
  df$contour_id <- paste0("Contour_", i)  # Add contour ID
  return(df)
}) %>% bind_rows()
#-------------------------------------------------------------------------------

# Extract first three PCAs
Y_sample = Curl_sim_data[, 5:14]
loc = Curl_sim_data[, 1:2]

z = Curl_sim_data$super_cluster

n = nrow(Y_sample) # number of observations
p = ncol(Y_sample) # number of variables in Y

# Standardize coordinates, Y values
coords <- apply(loc, 2, scale)
Y_std <- apply(Y_sample, 2, scale)

# Create dataframe for plotting
df_subset <- data.frame(coords, Y_std)

#-------------------------------------------------------------------------------
### Function to plot individual simulation results for Figure 3 in paper

# Define pallettes for plots
my_palette_curl <- c(
  "#4BA3C7", # bright sky blue
  "#F4A259", # bright apricot orange
  "#80C342", # bright fresh green
  "#EF476F", # lively coral pink
  "#A07BEF", # vivid lavender
  "#FFD166", # warm sunflower yellow
  "#06D6A0", # aqua mint
  "#118AB2", # deep cerulean blue
  "#FF7C43", # vibrant tangerine
  "#C77DFF"  # playful soft purple
)

my_palette_curl_truth <- c(
  "#89CFF0", # baby blue
  "#F08080", # light coral
  "#77DD77", # pastel green
  "#FFB347", # soft orange
  "#B19CD9", # light purple
  "#FFDAC1", # peach
  "#AEC6CF", # pastel cyan
  "#CFCFC4", # light stone
  "#E6A8D7", # orchid pink
  "#C5E384"  # pale lime
)
#-------------------------------------------------------------------------------

plot_points_paper <- function(coords, values, boundary_df, title = NULL, palette = my_palette_curl) {
  
  # Prepare point data
  pointsdata <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    Value = as.factor(values)  # Ensure factor for discrete colors
  )
  
  # Build ggplot
  plot <- ggplot() +
    # Add boundary paths
    geom_path(data = boundary_df, aes(x = x, y = y, group = contour_id),
              color = "black", size = 1) +
    
    # Add points
    geom_point(data = pointsdata, aes(x = x, y = y, color = Value), size = 2) +
    
    # Color palette
    scale_color_manual(values = palette) +
    
    # Clean theme
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    ) +
    
    # Title
    labs(title = title)
  
  return(plot)
}

#-------------------------------------------------------------------------------

# Define reference vector to assign same clusters to the same colors
reference_vector <- seq(1, 10, 1)

# Function to reorder cluster assignments based on the reference vector
reorder_based_on_reference <- function(cluster_vector, reference_vector) {
  map <- match(reference_vector, cluster_vector)
  reordered_vector <- map[cluster_vector]
  return(reordered_vector)
}

#-------------------------------------------------------------------------------

# Plot the ground truth for the Curl example
true_clusters <- plot_points_paper(
  coords = as.matrix(Curl_sim_data[, c("X", "Y")]),
  values = Curl_sim_data$super_cluster,
  boundary_df = all_contours,
  palette = my_palette_curl_truth,
  title = bquote(bold("Truth, ") ~ italic("c = 10"))
  # title = expression(bold("Truth"))
) +
  scale_y_reverse() +  # Optional flip Y-axis
  theme(
    plot.title = element_text(size = 20),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 24)
  )

# Display
true_clusters

# Add label to the plot for Figure 3
true_clusters_with_label <- ggdraw() +
  draw_plot(true_clusters) +
  draw_label("(A)", x = 0, y = 0.98, size = 20, hjust = 0, vjust = 1)
true_clusters_with_label

ggsave("./Simulations/Curl/Curl_plots/Curl_sim_true_clusters.png", 
       plot = true_clusters_with_label, 
       width = 5, height = 5, units = "in", dpi = 300, bg = "transparent")
#-------------------------------------------------------------------------------
##### DP-RST results for p=3 #####

load("./Simulations/Curl/Curl_results/split_DP-RST/Curl_sim_DP.RST_p3_FirstTen_reps.RData")
Curl_sim_DP.RST_p3_FirstTen_reps <- Curl_sim_DP.RST_p3_reps
load("./Simulations/Curl/Curl_results/split_DP-RST/Curl_sim_DP.RST_p3_SecondTen_reps.RData")
Curl_sim_DP.RST_p3_SecondTen_reps <- Curl_sim_DP.RST_p3_reps[11:20]
load("./Simulations/Curl/Curl_results/split_DP-RST/Curl_sim_DP.RST_p3_ThirdTen_reps.RData")
Curl_sim_DP.RST_p3_ThirdTen_reps <- Curl_sim_DP.RST_p3_reps[21:30]

# Combine all into one list
output_p3 <- c(Curl_sim_DP.RST_p3_FirstTen_reps, 
               Curl_sim_DP.RST_p3_SecondTen_reps, 
               Curl_sim_DP.RST_p3_ThirdTen_reps)

### Calculate the accuracy for each run ###
DPM_partition_p3 <- list()
DPM_accur_p3 = c()
clust_number_p3 <- c()

for (j in 1:30) {
  mode_based_partition = partition(DP.RST_output = output_p3[[j]], 
                                   method = "mode_based", batch_size = 100)
  DPM_partition_p3[[j]] <- mode_based_partition
  
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
  
  DPM_accur_p3[j] = adjustedRandIndex(Curl_sim_data$super_cluster, obs_in_teams_vec_mode_based)
  clust_number_p3[j] = length(unique(obs_in_teams_vec_mode_based))
}

# Plot the accuracy over runs
plot(DPM_accur_p3, type = "l")

table(clust_number_p3)
# clust_number_p3
# 5 
# 30 

summary(DPM_accur_p3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3310  0.3659  0.3756  0.3759  0.3847  0.4155 

## Get the cluster assignments of the highest accuracy
# which.max(DPM_accur_p3) # 24 for 5 clusters

original_index_p3 = which.max(DPM_accur_p3)

X <- table(sequence(length(DPM_partition_p3[[original_index_p3]]$groups_partition)), DPM_partition_p3[[original_index_p3]]$groups_partition)
Z <- table(sequence(length(DPM_partition_p3[[original_index_p3]]$teams_partition)), DPM_partition_p3[[original_index_p3]]$teams_partition)
obs_in_teams <- X %*% Z
DPM_partition_p3 <- obs_in_teams %*% sort(unique(DPM_partition_p3[[original_index_p3]]$teams_partition))

DPM_partition_p3 <- reorder_based_on_reference(DPM_partition_p3, reference_vector)
se_DPM_p3 <- sd(DPM_accur_p3) / sqrt(length(DPM_accur_p3))

ari_val <- round(max(DPM_accur_p3[which.max(DPM_accur_p3)]), 3)
se_val <- round(se_DPM_p3, 3)
c_val <- length(unique(DPM_partition_p3))

title_expr <- substitute(
  atop(
    bold("p = 3,") * italic(paste(" ARI = ", ARI)) * paste(" (", SE, ")", sep = ""),
    italic(paste("c = ", C))
  ),
  list(ARI = ari_val, SE = se_val, C = c_val)
)

DPM_results_p3 <- plot_points_paper(
  coords = as.matrix(Curl_sim_data[, c("X", "Y")]),
  values = DPM_partition_p3,
  boundary_df = all_contours,
  palette = my_palette_curl,
  title = title_expr
) +
  scale_y_reverse() +
  theme(plot.title = element_text(size = 20, lineheight = 1.2))

DPM_results_p3

#-------------------------------------------------------------------------------
##### DP-RST results for p=10 #####

load("./Simulations/Curl/Curl_results/split_DP-RST/Curl_sim_DP.RST_p10_FirstTen_30reps.RData")
Curl_sim_DP.RST_p10_FirstTen_reps <- Curl_sim_DP.RST_p10_reps
load("./Simulations/Curl/Curl_results/split_DP-RST/Curl_sim_DP.RST_p10_SecondTen_30reps.RData")
Curl_sim_DP.RST_p10_SecondTen_reps <- Curl_sim_DP.RST_p10_reps[11:20]
load("./Simulations/Curl/Curl_results/split_DP-RST/Curl_sim_DP.RST_p10_ThirdTen_30reps.RData")
Curl_sim_DP.RST_p10_ThirdTen_reps <- Curl_sim_DP.RST_p10_reps[21:30]

# Combine all into one list
output_p10 <- c(Curl_sim_DP.RST_p10_FirstTen_reps, 
               Curl_sim_DP.RST_p10_SecondTen_reps, 
               Curl_sim_DP.RST_p10_ThirdTen_reps)

### Calculate the accuracy for each run ###
DPM_partition_p10 <- list()
DPM_accur_p10 = c()
clust_number_p10 <- c()

for (j in 1:30) {
  mode_based_partition = partition(DP.RST_output = output_p10[[j]], method = "mode_based", batch_size = 100)
  DPM_partition_p10[[j]] <- mode_based_partition
  
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
  
  DPM_accur_p10[j] = adjustedRandIndex(Curl_sim_data$super_cluster, obs_in_teams_vec_mode_based)
  clust_number_p10[j] = length(unique(obs_in_teams_vec_mode_based))
}

# Plot the accuracy over runs
plot(DPM_accur_p10, type = "l")
summary(DPM_accur_p10)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5476  0.5883  0.6078  0.6096  0.6298  0.6706 

table(clust_number_p10)
# clust_number_p10
# 7  8  9 
# 17 12  1 

## Get the cluster assignments of the highest accuracy
# which.max(DPM_accur_p10) # 24 for 5 clusters

original_index_p10 <- which.max(DPM_accur_p10)

X <- table(sequence(length(DPM_partition_p10[[original_index_p10]]$groups_partition)), DPM_partition_p10[[original_index_p10]]$groups_partition)
Z <- table(sequence(length(DPM_partition_p10[[original_index_p10]]$teams_partition)), DPM_partition_p10[[original_index_p10]]$teams_partition)
obs_in_teams <- X %*% Z
DPM_partition_p10 <- obs_in_teams %*% sort(unique(DPM_partition_p10[[original_index_p10]]$teams_partition))

DPM_partition_p10 <- reorder_based_on_reference(DPM_partition_p10, reference_vector)
se_DPM_p10 <- sd(DPM_accur_p10) / sqrt(length(DPM_accur_p10))

ari_val <- round(max(DPM_accur_p10[which.max(DPM_accur_p10)]), 3)
se_val <- round(se_DPM_p10, 3)
c_val <- length(unique(DPM_partition_p10))

title_expr <- substitute(
  atop(
    bold("p = 10,") * italic(paste(" ARI = ", ARI)) * paste(" (", SE, ")", sep = ""),
    italic(paste("c = ", C))
  ),
  list(ARI = ari_val, SE = se_val, C = c_val)
)

DPM_results_p10 <- plot_points_paper(
  coords = as.matrix(Curl_sim_data[, c("X", "Y")]),
  values = DPM_partition_p10,
  boundary_df = all_contours,
  palette = my_palette_curl,
  title = title_expr
) +
  scale_y_reverse() +
  theme(plot.title = element_text(size = 20, lineheight = 1.2))

DPM_results_p10


#-------------------------------------------------------------------------------

combined_plot <- DPM_results_p3 + DPM_results_p10 +
  plot_annotation(title = bquote(bold("DP-RST"))) &
  theme(
    plot.title = element_text(size = 24, hjust = 0.5)  # Center title
  )

combined_plot

DPM_results_with_label <- ggdraw() +
  draw_plot(combined_plot) +
  draw_label("(i)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
DPM_results_with_label

DPM_results_with_label_MainPaper <- ggdraw() +
  draw_plot(combined_plot) +
  draw_label("(B)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
DPM_results_with_label_MainPaper

ggsave("./Simulations/Curl/Curl_plots/Curl_sim_DPM.png", 
       plot = DPM_results_with_label, 
       width = 10, height = 5.5, units = "in", dpi = 300, bg = "transparent")

# ggsave("./Simulations/Curl/Curl_plots/Curl_sim_DPM_new_MainPaper.png", 
#        plot = DPM_results_with_label_MainPaper, 
#        width = 10, height = 5.5, units = "in", dpi = 300, bg = "transparent")

#-------------------------------------------------------------------------------

### Make a histogram with the number of teams chosen over runs (for the Figure 1 in Supplementary Materials)
clust_number_p3 <- data.frame(clust_number = clust_number_p3)
clust_number_p3$clust_number <- factor(clust_number_p3$clust_number, levels = 1:10)

Curl_hist_plot_p3 <- ggplot(clust_number_p3, aes(x = clust_number)) +
  geom_bar(fill = "#C77DFF", color = "black") +
  theme_minimal() +
  labs(
    title = "Cluster Number Estimation in Curl Simulation
    p = 3",
    x = "Estimated number of clusters",
    y = "Frequency"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

Curl_hist_plot_p3


clust_number_p10 <- data.frame(clust_number = clust_number_p10)
clust_number_p10$clust_number <- factor(clust_number_p10$clust_number, levels = 1:11)

Curl_hist_plot_p10 <- ggplot(clust_number_p10, aes(x = clust_number)) +
  geom_bar(fill = "#EF476F", color = "black") +
  theme_minimal() +
  labs(
    title = "Cluster Number Estimation in Curl Simulation
    p = 10",
    x = "Estimated number of clusters",
    y = "Frequency"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

Curl_hist_plot_p10

combined_plot_hist <- Curl_hist_plot_p3 + Curl_hist_plot_p10
combined_plot_hist

ggsave("./Simulations/Curl/Curl_plots/Curl_hist_combined_plot.png", combined_plot_hist, width = 10, height = 4)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
##### BayesSpace results #####
load("./Simulations/Curl/Curl_results/Curl_sim_BayesSpace_p3_30reps.RData")

BayesSpace_acc_p3 <- c()
for (i in 1:30){
  BayesSpace_acc_p3[i] <- adjustedRandIndex(Curl_sim_data$super_cluster, Curl_sim_BayesSpace_p3_reps[[i]]$spatial.cluster)
  print(table(Curl_sim_BayesSpace_p3_reps[[i]]$spatial.cluster))
}

plot(BayesSpace_acc_p3, type = "l")
summary(BayesSpace_acc_p3)

BayesSpace_index_p3 <- which.max(BayesSpace_acc_p3)

BS_partition_p3 <- reorder_based_on_reference(Curl_sim_BayesSpace_p3_reps[[BayesSpace_index_p3]]$spatial.cluster, reference_vector)
se_BS_p3 <- sd(BayesSpace_acc_p3) / sqrt(length(BayesSpace_acc_p3))

BayesSpace_results_p3 <- plot_points_paper(
  coords = as.matrix(Curl_sim_data[, c("X", "Y")]),
  values = BS_partition_p3,
  boundary_df = all_contours,
  palette = my_palette_curl,
  title = bquote(bold("p = 3, ") ~ italic("ARI = ") ~ .(round(max(BayesSpace_acc_p3[BayesSpace_index_p3]), 3)) ~ " (" ~ .(round(se_BS_p3, 3)) ~ ")")
) +
  scale_y_reverse() +
  theme(plot.title = element_text(size = 20))
  
BayesSpace_results_p3

#-------------------------------------------------------------------------------
load("./Simulations/Curl/Curl_results/Curl_sim_BayesSpace_p10_30reps.RData")

BayesSpace_acc_p10 <- c()
for (i in 1:30){
  BayesSpace_acc_p10[i] <- adjustedRandIndex(Curl_sim_data$super_cluster, Curl_sim_BayesSpace_p10_reps[[i]]$spatial.cluster)
  print(table(Curl_sim_BayesSpace_p10_reps[[i]]$spatial.cluster))
}

plot(BayesSpace_acc_p10, type = "l")
summary(BayesSpace_acc_p10)

BayesSpace_index_p10 <- which.max(BayesSpace_acc_p10)

BS_partition_p10 <- reorder_based_on_reference(Curl_sim_BayesSpace_p10_reps[[BayesSpace_index_p10]]$spatial.cluster, reference_vector)
se_BS_p10 <- sd(BayesSpace_acc_p10) / sqrt(length(BayesSpace_acc_p10))

BayesSpace_results_p10 <- plot_points_paper(
  coords = as.matrix(Curl_sim_data[, c("X", "Y")]),
  values = BS_partition_p10,
  boundary_df = all_contours,
  palette = my_palette_curl,
  title = bquote(bold("p = 10, ") ~ italic("ARI = ") ~ .(round(max(BayesSpace_acc_p10[BayesSpace_index_p10]), 3)) ~ " (" ~ .(round(se_BS_p10, 3)) ~ ")")
) +
  scale_y_reverse() +
  theme(plot.title = element_text(size = 20))
BayesSpace_results_p10

#-------------------------------------------------------------------------------

combined_plot_BS <- BayesSpace_results_p3 + BayesSpace_results_p10 +
  plot_annotation(title = bquote(bold("BayesSpace"))) &
  theme(
    plot.title = element_text(size = 24, hjust = 0.5)  # Center title
  )

combined_plot_BS

BayesSpace_results_with_label <- ggdraw() +
  draw_plot(combined_plot_BS) +
  draw_label("(ii)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
BayesSpace_results_with_label


ggsave("./Simulations/Curl/Curl_plots/Curl_sim_BayesSpace.png", 
       plot = BayesSpace_results_with_label, 
       width = 10, height = 5, units = "in", dpi = 300, bg = "transparent")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
##### SC-MEB results #####
load("./Simulations/Curl/Curl_results/Curl_sim_SC.MEB_p3_30reps.RData")

SC.MEB_acc_p3 <- c()
for (i in 1:30){
  SC.MEB_acc_p3[i] <- adjustedRandIndex(Curl_sim_data$super_cluster, Curl_sim_SC.MEB_p3_reps[[i]][[1]])
}

plot(SC.MEB_acc_p3, type = "l")
summary(SC.MEB_acc_p3)

SC.MEB_index_p3 <- which.max(SC.MEB_acc_p3)

SC.MEB_partition_p3 <- reorder_based_on_reference(Curl_sim_SC.MEB_p3_reps[[SC.MEB_index_p3]][[1]], reference_vector)
se_SC.MEB_p3 <- sd(SC.MEB_acc_p3) / sqrt(length(SC.MEB_acc_p3))

SC.MEB_results_p3 <- plot_points_paper(
  coords = as.matrix(Curl_sim_data[, c("X", "Y")]),
  values = SC.MEB_partition_p3,
  boundary_df = all_contours,
  palette = my_palette_curl,
  title = bquote(bold("p = 3, ") ~ italic("ARI = ") ~ .(round(max(SC.MEB_acc_p3[SC.MEB_index_p3]), 3)) ~ " (" ~ .(round(se_SC.MEB_p3, 3)) ~ ")")
) +
  scale_y_reverse() +
  theme(plot.title = element_text(size = 20))

SC.MEB_results_p3

#-------------------------------------------------------------------------------
load("./Simulations/Curl/Curl_results/Curl_sim_SC.MEB_p10_30reps.RData")

SC.MEB_acc_p10 <- c()
for (i in 1:30){
  SC.MEB_acc_p10[i] <- adjustedRandIndex(Curl_sim_data$super_cluster, Curl_sim_SC.MEB_p10_reps[[i]][[1]])
}

plot(SC.MEB_acc_p10, type = "l")
summary(SC.MEB_acc_p10)

SC.MEB_index_p10 <- which.max(SC.MEB_acc_p10)

SC.MEB_partition_p10 <- reorder_based_on_reference(Curl_sim_SC.MEB_p10_reps[[SC.MEB_index_p10]][[1]], reference_vector)
se_SC.MEB_p10 <- sd(SC.MEB_acc_p10) / sqrt(length(SC.MEB_acc_p10))

SC.MEB_results_p10 <- plot_points_paper(
  coords = as.matrix(Curl_sim_data[, c("X", "Y")]),
  values = SC.MEB_partition_p10,
  boundary_df = all_contours,
  palette = my_palette_curl,
  title = bquote(bold("p = 10, ") ~ italic("ARI = ") ~ .(round(max(SC.MEB_acc_p10[SC.MEB_index_p10]), 3)) ~ " (" ~ .(round(se_SC.MEB_p10, 3)) ~ ")")
) +
  scale_y_reverse() +
  theme(plot.title = element_text(size = 20))

SC.MEB_results_p10

#-------------------------------------------------------------------------------

combined_plot_SC.MEB <- SC.MEB_results_p3 + SC.MEB_results_p10 +
  plot_annotation(title = bquote(bold("SC-MEB"))) &
  theme(
    plot.title = element_text(size = 24, hjust = 0.5)  # Center title
  )

combined_plot_SC.MEB

SC.MEB_results_with_label <- ggdraw() +
  draw_plot(combined_plot_SC.MEB) +
  draw_label("(iii)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
SC.MEB_results_with_label


ggsave("./Simulations/Curl/Curl_plots/Curl_sim_SC.MEB.png", 
       plot = SC.MEB_results_with_label, 
       width = 10, height = 5, units = "in", dpi = 300, bg = "transparent")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
##### DR-RC results #####
load("./Simulations/Curl/Curl_results/Curl_sim_DR.SC_p3_30reps.RData")

DR.SC_acc_p3 <- c()
for (i in 1:30){
  DR.SC_acc_p3[i] <- adjustedRandIndex(Curl_sim_data$super_cluster, Curl_sim_DR.SC_p3_reps[[i]][["Objdrsc"]][[1]][["cluster"]])
}

plot(DR.SC_acc_p3, type = "l")
summary(DR.SC_acc_p3)

DR.SC_index_p3 <- which.max(DR.SC_acc_p3)

DR.SC_partition_p3 <- reorder_based_on_reference(Curl_sim_DR.SC_p3_reps[[DR.SC_index_p3]][["Objdrsc"]][[1]][["cluster"]], reference_vector)
se_DR.SC_p3 <- sd(DR.SC_acc_p3) / sqrt(length(DR.SC_acc_p3))

DR.SC_results_p3 <- plot_points_paper(
  coords = as.matrix(Curl_sim_data[, c("X", "Y")]),
  values = DR.SC_partition_p3,
  boundary_df = all_contours,
  palette = my_palette_curl,
  title = bquote(bold("p = 3, ") ~ italic("ARI = ") ~ .(round(max(DR.SC_acc_p3[DR.SC_index_p3]), 3)) ~ " (" ~ .(round(se_DR.SC_p3, 3)) ~ ")")
) +
  scale_y_reverse() +
  theme(plot.title = element_text(size = 20))

DR.SC_results_p3

#-------------------------------------------------------------------------------
load("./Simulations/Curl/Curl_results/Curl_sim_DR.SC_p10_30reps.RData")

DR.SC_acc_p10 <- c()
for (i in 1:30){
  DR.SC_acc_p10[i] <- adjustedRandIndex(Curl_sim_data$super_cluster, Curl_sim_DR.SC_p10_reps[[i]][["Objdrsc"]][[1]][["cluster"]])
}

plot(DR.SC_acc_p10, type = "l")
summary(DR.SC_acc_p10)

DR.SC_index_p10 <- which.max(DR.SC_acc_p10)

DR.SC_partition_p10 <- reorder_based_on_reference(Curl_sim_DR.SC_p10_reps[[DR.SC_index_p10]][["Objdrsc"]][[1]][["cluster"]], reference_vector)
se_DR.SC_p10 <- sd(DR.SC_acc_p10) / sqrt(length(DR.SC_acc_p10))

DR.SC_results_p10 <- plot_points_paper(
  coords = as.matrix(Curl_sim_data[, c("X", "Y")]),
  values = DR.SC_partition_p10,
  boundary_df = all_contours,
  palette = my_palette_curl,
  title = bquote(bold("p = 10, ") ~ italic("ARI = ") ~ .(round(max(DR.SC_acc_p10[DR.SC_index_p10]), 3)) ~ " (" ~ .(round(se_DR.SC_p10, 3)) ~ ")")
) +
  scale_y_reverse() +
  theme(plot.title = element_text(size = 20))

DR.SC_results_p10

#-------------------------------------------------------------------------------

combined_plot_DR.SC <- DR.SC_results_p3 + DR.SC_results_p10 +
  plot_annotation(title = bquote(bold("DR-SC"))) &
  theme(
    plot.title = element_text(size = 24, hjust = 0.5)  # Center title
  )

combined_plot_DR.SC

DR.SC_results_with_label <- ggdraw() +
  draw_plot(combined_plot_DR.SC) +
  draw_label("(iv)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
DR.SC_results_with_label


ggsave("./Simulations/Curl/Curl_plots/Curl_sim_DR.SC.png", 
       plot = DR.SC_results_with_label, 
       width = 10, height = 5, units = "in", dpi = 300, bg = "transparent")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
##### k-means results #####
load("./Simulations/Curl/Curl_results/Curl_sim_kmeans_3p_30reps.RData")

kmeans_acc_p3 <- c()
for (i in 1:30){
  kmeans_acc_p3[i] <- adjustedRandIndex(Curl_sim_data$super_cluster, Curl_sim_kmeans_3p_reps[[i]])
}

plot(kmeans_acc_p3, type = "l")
summary(kmeans_acc_p3)

kmeans_index_p3 <- which.max(kmeans_acc_p3)

kmeans_partition_p3 <- reorder_based_on_reference(Curl_sim_kmeans_3p_reps[[kmeans_index_p3]], reference_vector)
se_kmeans_p3 <- sd(kmeans_acc_p3) / sqrt(length(kmeans_acc_p3))

kmeans_results_p3 <- plot_points_paper(
  coords = as.matrix(Curl_sim_data[, c("X", "Y")]),
  values = kmeans_partition_p3,
  boundary_df = all_contours,
  palette = my_palette_curl,
  title = bquote(bold("p = 3, ") ~ italic("ARI = ") ~ .(round(max(kmeans_acc_p3[kmeans_index_p3]), 3)) ~ " (" ~ .(round(se_kmeans_p3, 3)) ~ ")")
) +
  scale_y_reverse() +
  theme(plot.title = element_text(size = 20))

kmeans_results_p3

#-------------------------------------------------------------------------------
load("./Simulations/Curl/Curl_results/Curl_sim_kmeans_10p_30reps.RData")

kmeans_acc_p10 <- c()
for (i in 1:30){
  kmeans_acc_p10[i] <- adjustedRandIndex(Curl_sim_data$super_cluster, Curl_sim_kmeans_10p_reps[[i]])
}

plot(kmeans_acc_p10, type = "l")
summary(kmeans_acc_p10)

kmeans_index_p10 <- which.max(kmeans_acc_p10)

kmeans_partition_p10 <- reorder_based_on_reference(Curl_sim_kmeans_10p_reps[[kmeans_index_p10]], reference_vector)
se_kmeans_p10 <- sd(kmeans_acc_p10) / sqrt(length(kmeans_acc_p10))

kmeans_results_p10 <- plot_points_paper(
  coords = as.matrix(Curl_sim_data[, c("X", "Y")]),
  values = kmeans_partition_p10,
  boundary_df = all_contours,
  palette = my_palette_curl,
  title = bquote(bold("p = 10, ") ~ italic("ARI = ") ~ .(round(max(kmeans_acc_p10[kmeans_index_p10]), 3)) ~ " (" ~ .(round(se_kmeans_p10, 3)) ~ ")")
) +
  scale_y_reverse() +
  theme(plot.title = element_text(size = 20))

kmeans_results_p10

#-------------------------------------------------------------------------------

combined_plot_kmeans <- kmeans_results_p3 + kmeans_results_p10 +
  plot_annotation(title = bquote(bold("k-means"))) &
  theme(
    plot.title = element_text(size = 24, hjust = 0.5)  # Center title
  )

combined_plot_kmeans

kmeans_results_with_label <- ggdraw() +
  draw_plot(combined_plot_kmeans) +
  draw_label("(v)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
kmeans_results_with_label


ggsave("./Simulations/Curl/Curl_plots/Curl_sim_kmeans.png", 
       plot = kmeans_results_with_label, 
       width = 10, height = 5, units = "in", dpi = 300, bg = "transparent")


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
##### Load GraphST data #####
GraphST_sim_results_p3 <- readr::read_csv("./Simulations/Curl/Curl_results/graphst_curl_3p_multiseed_k10_radius50_results.csv")

# Extract ARIs from all GraphST columns
graphst_partition_matrix_p3 <- GraphST_sim_results_p3 %>%
  select(starts_with("refined_mclust_seed_")) %>%
  as.matrix()

graphst_ari_p3 <- apply(
  graphst_partition_matrix_p3, 
  2, 
  function(col) adjustedRandIndex(GraphST_sim_results_p3$original_cluster, col)
)

graphst_best_index_p3 <- which.max(graphst_ari_p3)
graphst_best_ari_p3 <- graphst_ari_p3[graphst_best_index_p3]
graphst_se_p3 <- sd(graphst_ari_p3) / sqrt(length(graphst_ari_p3))

# Best partition
graphst_best_partition_p3 <- reorder_based_on_reference(
  graphst_partition_matrix_p3[, graphst_best_index_p3],
  reference_vector
)

graphst_results_p3 <- plot_points_paper(
  coords = as.matrix(Curl_sim_data[, c("X", "Y")]),
  values = graphst_best_partition_p3,
  boundary_df = all_contours,
  palette = my_palette_curl,
  title = bquote(bold("p = 3, ") ~ italic("ARI = ") ~ .(round(graphst_best_ari_p3, 3)) ~ " (" ~ .(round(graphst_se_p3, 3)) ~ ")")
) +
  scale_y_reverse() +
  theme(plot.title = element_text(size = 20))

graphst_results_p3
#-------------------------------------------------------------------------------

GraphST_sim_results_p10 <- readr::read_csv("./Simulations/Curl/Curl_results/graphst_curl_10p_multiseed_k10_radius50_results.csv")

# Extract ARIs from all GraphST columns
graphst_partition_matrix_p10 <- GraphST_sim_results_p10 %>%
  select(starts_with("refined_mclust_seed_")) %>%
  as.matrix()

graphst_ari_p10 <- apply(
  graphst_partition_matrix_p10, 
  2, 
  function(col) adjustedRandIndex(GraphST_sim_results_p10$original_cluster, col)
)

graphst_best_index_p10 <- which.max(graphst_ari_p10)
graphst_best_ari_p10 <- graphst_ari_p10[graphst_best_index_p10]
graphst_se_p10 <- sd(graphst_ari_p10) / sqrt(length(graphst_ari_p10))

# Best partition
graphst_best_partition_p10 <- reorder_based_on_reference(
  graphst_partition_matrix_p10[, graphst_best_index_p10],
  reference_vector
)

graphst_results_p10 <- plot_points_paper(
  coords = as.matrix(Curl_sim_data[, c("X", "Y")]),
  values = graphst_best_partition_p10,
  boundary_df = all_contours,
  palette = my_palette_curl,
  title = bquote(bold("p = 10, ") ~ italic("ARI = ") ~ .(round(graphst_best_ari_p10, 3)) ~ " (" ~ .(round(graphst_se_p10, 3)) ~ ")")
) +
  scale_y_reverse() +
  theme(plot.title = element_text(size = 20))

graphst_results_p10

#-------------------------------------------------------------------------------

combined_plot_graphst <- graphst_results_p3 + graphst_results_p10 +
  plot_annotation(title = bquote(bold("GraphST"))) &
  theme(
    plot.title = element_text(size = 24, hjust = 0.5)  # Center title
  )

combined_plot_graphst

graphst_results_with_label <- ggdraw() +
  draw_plot(combined_plot_graphst) +
  draw_label("(vi)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
graphst_results_with_label


ggsave("./Simulations/Curl/Curl_plots/Curl_sim_GraphST.png", 
       plot = graphst_results_with_label, 
       width = 10, height = 5, units = "in", dpi = 300, bg = "transparent")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
##### Load SEDR data #####
SEDR_sim_results_p3 <- readr::read_csv("./Simulations/Curl/Curl_results/sedr_curl_3p_multiseed_k10_results_BestResult.csv")

# Extract ARIs from all SEDR columns
sedr_partition_matrix_p3 <- SEDR_sim_results_p3 %>%
  select(starts_with("mclust_seed_")) %>%
  as.matrix()

sedr_ari_p3 <- apply(
  sedr_partition_matrix_p3, 
  2, 
  function(col) adjustedRandIndex(SEDR_sim_results_p3$original_cluster, col)
)

sedr_best_index_p3 <- which.max(sedr_ari_p3)
sedr_best_ari_p3 <- sedr_ari_p3[sedr_best_index_p3]
sedr_se_p3 <- sd(sedr_ari_p3) / sqrt(length(sedr_ari_p3))

# Best partition
sedr_best_partition_p3 <- reorder_based_on_reference(
  sedr_partition_matrix_p3[, sedr_best_index_p3],
  reference_vector
)

sedr_results_p3 <- plot_points_paper(
  coords = as.matrix(Curl_sim_data[, c("X", "Y")]),
  values = sedr_best_partition_p3,
  boundary_df = all_contours,
  palette = my_palette_curl,
  title = bquote(bold("p = 3, ") ~ italic("ARI = ") ~ .(round(sedr_best_ari_p3, 3)) ~ " (" ~ .(round(sedr_se_p3, 3)) ~ ")")
) +
  scale_y_reverse() +
  theme(plot.title = element_text(size = 20))

sedr_results_p3
#-------------------------------------------------------------------------------

SEDR_sim_results_p10 <- readr::read_csv("./Simulations/Curl/Curl_results/sedr_curl_10p_multiseed_k10_results_BestResult.csv")

# Extract ARIs from all SEDR columns
sedr_partition_matrix_p10 <- SEDR_sim_results_p10 %>%
  select(starts_with("mclust_seed_")) %>%
  as.matrix()

sedr_ari_p10 <- apply(
  sedr_partition_matrix_p10, 
  2, 
  function(col) adjustedRandIndex(SEDR_sim_results_p10$original_cluster, col)
)

sedr_best_index_p10 <- which.max(sedr_ari_p10)
sedr_best_ari_p10 <- sedr_ari_p10[sedr_best_index_p10]
sedr_se_p10 <- sd(sedr_ari_p10) / sqrt(length(sedr_ari_p10))

# Best partition
sedr_best_partition_p10 <- reorder_based_on_reference(
  sedr_partition_matrix_p10[, sedr_best_index_p10],
  reference_vector
)

sedr_results_p10 <- plot_points_paper(
  coords = as.matrix(Curl_sim_data[, c("X", "Y")]),
  values = sedr_best_partition_p10,
  boundary_df = all_contours,
  palette = my_palette_curl,
  title = bquote(bold("p = 10, ") ~ italic("ARI = ") ~ .(round(sedr_best_ari_p10, 3)) ~ " (" ~ .(round(sedr_se_p10, 3)) ~ ")")
) +
  scale_y_reverse() +
  theme(plot.title = element_text(size = 20))

sedr_results_p10

#-------------------------------------------------------------------------------

combined_plot_sedr <- sedr_results_p3 + sedr_results_p10 +
  plot_annotation(title = bquote(bold("SEDR"))) &
  theme(
    plot.title = element_text(size = 24, hjust = 0.5)  # Center title
  )

combined_plot_sedr

sedr_results_with_label <- ggdraw() +
  draw_plot(combined_plot_sedr) +
  draw_label("(vii)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
sedr_results_with_label


ggsave("./Simulations/Curl/Curl_plots/Curl_sim_SEDR.png", 
       plot = sedr_results_with_label, 
       width = 10, height = 5, units = "in", dpi = 300, bg = "transparent")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
##### Load STAGATE data #####
STAGATE_sim_results_p3 <- readr::read_csv("./Simulations/Curl/Curl_results/stagate_simdata_multiseed_10clusters_3p_results.csv")

# Extract ARIs from all STAGATE columns
stagate_partition_matrix_p3 <- STAGATE_sim_results_p3 %>%
  select(starts_with("stagate_seed_")) %>%
  as.matrix()

stagate_ari_p3 <- apply(
  stagate_partition_matrix_p3, 
  2, 
  function(col) adjustedRandIndex(STAGATE_sim_results_p3$cluster, col)
)

stagate_best_index_p3 <- which.max(stagate_ari_p3)
stagate_best_ari_p3 <- stagate_ari_p3[stagate_best_index_p3]
stagate_se_p3 <- sd(stagate_ari_p3) / sqrt(length(stagate_ari_p3))

# Best partition
stagate_best_partition_p3 <- reorder_based_on_reference(
  stagate_partition_matrix_p3[, stagate_best_index_p3],
  reference_vector
)

stagate_results_p3 <- plot_points_paper(
  coords = as.matrix(Curl_sim_data[, c("X", "Y")]),
  values = stagate_best_partition_p3,
  boundary_df = all_contours,
  palette = my_palette_curl,
  title = bquote(bold("p = 3, ") ~ italic("ARI = ") ~ .(round(stagate_best_ari_p3, 3)) ~ " (" ~ .(round(stagate_se_p3, 3)) ~ ")")
) +
  scale_y_reverse() +
  theme(plot.title = element_text(size = 20))

stagate_results_p3
#-------------------------------------------------------------------------------

STAGATE_sim_results_p10 <- readr::read_csv("./Simulations/Curl/Curl_results/stagate_simdata_multiseed_10clusters_10p_results.csv")

# Extract ARIs from all STAGATE columns
stagate_partition_matrix_p10 <- STAGATE_sim_results_p10 %>%
  select(starts_with("stagate_seed_")) %>%
  as.matrix()

stagate_ari_p10 <- apply(
  stagate_partition_matrix_p10, 
  2, 
  function(col) adjustedRandIndex(STAGATE_sim_results_p10$cluster, col)
)

stagate_best_index_p10 <- which.max(stagate_ari_p10)
stagate_best_ari_p10 <- stagate_ari_p10[stagate_best_index_p10]
stagate_se_p10 <- sd(stagate_ari_p10) / sqrt(length(stagate_ari_p10))

# Best partition
stagate_best_partition_p10 <- reorder_based_on_reference(
  stagate_partition_matrix_p10[, stagate_best_index_p10],
  reference_vector
)

stagate_results_p10 <- plot_points_paper(
  coords = as.matrix(Curl_sim_data[, c("X", "Y")]),
  values = stagate_best_partition_p10,
  boundary_df = all_contours,
  palette = my_palette_curl,
  title = bquote(bold("p = 10, ") ~ italic("ARI = ") ~ .(round(stagate_best_ari_p10, 3)) ~ " (" ~ .(round(stagate_se_p10, 3)) ~ ")")
) +
  scale_y_reverse() +
  theme(plot.title = element_text(size = 20))

stagate_results_p10

#-------------------------------------------------------------------------------

combined_plot_stagate <- stagate_results_p3 + stagate_results_p10 +
  plot_annotation(title = bquote(bold("STAGATE"))) &
  theme(
    plot.title = element_text(size = 24, hjust = 0.5)  # Center title
  )

combined_plot_stagate

stagate_results_with_label <- ggdraw() +
  draw_plot(combined_plot_stagate) +
  draw_label("(viii)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
stagate_results_with_label


ggsave("./Simulations/Curl/Curl_plots/Curl_sim_STAGATE.png", 
       plot = stagate_results_with_label, 
       width = 10, height = 5, units = "in", dpi = 300, bg = "transparent")


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
##### Load SpaGCN data #####
SpaGCN_sim_results_p3 <- readr::read_csv("./Simulations/Curl/Curl_results/curl_spagcn_clustering_results_3p_multi_seed.csv")

# Extract ARIs from all SpaGCN columns
spagcn_partition_matrix_p3 <- SpaGCN_sim_results_p3 %>%
  select(starts_with("refined_pred_seed_")) %>%
  as.matrix()

spagcn_ari_p3 <- apply(
  spagcn_partition_matrix_p3, 
  2, 
  function(col) adjustedRandIndex(SpaGCN_sim_results_p3$ground_truth_cluster, col)
)

spagcn_best_index_p3 <- which.max(spagcn_ari_p3)
spagcn_best_ari_p3 <- spagcn_ari_p3[spagcn_best_index_p3]
spagcn_se_p3 <- sd(spagcn_ari_p3) / sqrt(length(spagcn_ari_p3))

# Best partition
spagcn_best_partition_p3 <- reorder_based_on_reference(
  spagcn_partition_matrix_p3[, spagcn_best_index_p3]+1,
  reference_vector
)

spagcn_results_p3 <- plot_points_paper(
  coords = as.matrix(Curl_sim_data[, c("X", "Y")]),
  values = spagcn_best_partition_p3,
  boundary_df = all_contours,
  palette = my_palette_curl,
  title = bquote(bold("p = 3, ") ~ italic("ARI = ") ~ .(round(spagcn_best_ari_p3, 3)) ~ " (" ~ .(round(spagcn_se_p3, 3)) ~ ")")
) +
  scale_y_reverse() +
  theme(plot.title = element_text(size = 20))

spagcn_results_p3

#-------------------------------------------------------------------------------
SpaGCN_sim_results_p10 <- readr::read_csv("./Simulations/Curl/Curl_results/curl_spagcn_clustering_results_10p_multi_seed.csv")

# Extract ARIs from all SpaGCN columns
spagcn_partition_matrix_p10 <- SpaGCN_sim_results_p10 %>%
  select(starts_with("refined_pred_seed_")) %>%
  as.matrix()

spagcn_ari_p10 <- apply(
  spagcn_partition_matrix_p10, 
  2, 
  function(col) adjustedRandIndex(SpaGCN_sim_results_p10$ground_truth_cluster, col)
)

spagcn_best_index_p10 <- which.max(spagcn_ari_p10)
spagcn_best_ari_p10 <- spagcn_ari_p10[spagcn_best_index_p10]
spagcn_se_p10 <- sd(spagcn_ari_p10) / sqrt(length(spagcn_ari_p10))

# Best partition
spagcn_best_partition_p10 <- reorder_based_on_reference(
  spagcn_partition_matrix_p10[, spagcn_best_index_p10]+1,
  reference_vector
)

spagcn_results_p10 <- plot_points_paper(
  coords = as.matrix(Curl_sim_data[, c("X", "Y")]),
  values = spagcn_best_partition_p10,
  boundary_df = all_contours,
  palette = my_palette_curl,
  title = bquote(bold("p = 10, ") ~ italic("ARI = ") ~ .(round(spagcn_best_ari_p10, 3)) ~ " (" ~ .(round(spagcn_se_p10, 3)) ~ ")")
) +
  scale_y_reverse() +
  theme(plot.title = element_text(size = 20))

spagcn_results_p10

#-------------------------------------------------------------------------------

combined_plot_spagcn <- spagcn_results_p3 + spagcn_results_p10 +
  plot_annotation(title = bquote(bold("SpaGCN"))) &
  theme(
    plot.title = element_text(size = 24, hjust = 0.5)  # Center title
  )

combined_plot_spagcn

spagcn_results_with_label <- ggdraw() +
  draw_plot(combined_plot_spagcn) +
  draw_label("(ix)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
spagcn_results_with_label


ggsave("./Simulations/Curl/Curl_plots/Curl_sim_SpaGCN.png", 
       plot = spagcn_results_with_label, 
       width = 10, height = 5, units = "in", dpi = 300, bg = "transparent")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
##### Create box-plot for accuracy results for all methods#####
# Combine the vectors into a data frame
# Combine all accuracies into one data frame
data <- data.frame(
  Accuracy = c(DPM_accur_p3, DPM_accur_p10,
               BayesSpace_acc_p3, BayesSpace_acc_p10,
               DR.SC_acc_p3, DR.SC_acc_p10,
               SC.MEB_acc_p3, SC.MEB_acc_p10,
               kmeans_acc_p3, kmeans_acc_p10,
               graphst_ari_p3, graphst_ari_p10,
               sedr_ari_p3, sedr_ari_p10,
               stagate_ari_p3, stagate_ari_p10,
               spagcn_ari_p3, spagcn_ari_p10),
  Method = factor(rep(c("DP-RST", "Bayes\nSpace", "DR.SC", "SC.MEB", "k-means", 
                        "GraphST", "SEDR", "STAGATE", "SpaGCN"), each = 60),
                  levels = c("DP-RST", "Bayes\nSpace", "DR.SC", "SC.MEB", "k-means",
                             "GraphST", "SEDR", "STAGATE", "SpaGCN")),
  P = factor(rep(rep(c("p = 3", "p = 10"), each = 30), times = 9),
             levels = c("p = 3", "p = 10"))
)

# Create the box plot
boxplot <- ggplot(data, aes(x = Method, y = Accuracy, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8), aes(group = interaction(Method, P))) +
  facet_wrap(~P, ncol = 1) +  # Optional: if you want facets. If you prefer grouped bars, remove.
  scale_fill_manual(values = c(
    "DP-RST" = "#C5E384",
    "Bayes\nSpace" = "#E6A8D7",
    "DR.SC" = "#89CFF0",
    "SC.MEB" = "#FFB347",
    "k-means" = "#FFDAC1",
    "GraphST"     = "#A3A0FB",  # periwinkle purple
    "SEDR"        = "#77DD77",  # pastel green
    "STAGATE"     = "#FF6961",  # soft coral red
    "SpaGCN"      = "#779ECB"   # muted denim blue
  )) +
  theme_minimal() +
  labs(x = NULL,
       y = "ARI") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 20),
    plot.margin = unit(c(2.5, 2.5, 2.5, 2.5), "cm")
  )

boxplot

boxplot_with_label <- ggdraw() +
  draw_label("Compare ARI", x = 0.5, y = 0.98, size = 20, fontface = "bold", hjust = 0.5, vjust = 1) +
  draw_plot(boxplot) +
  draw_label("(B)", x = 0.1, y = 0.98, size = 20, hjust = 0, vjust = 1)
boxplot_with_label

ggsave("./Simulations/Curl/Curl_plots/Curl_sim_boxplot_BestSEDR.png", 
       plot = boxplot_with_label, 
       width = 12, height = 8, units = "in", dpi = 300, bg = "transparent")

