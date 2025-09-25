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
Ellipse_sim_data <- readr::read_csv("./New_Simulations/Ellipse/Ellipse_sim_data.csv")

#-------------------------------------------------------------------------------

# Extract first three PCAs
Y_sample = Ellipse_sim_data[, 5:14]
loc = Ellipse_sim_data[, 1:2]

z = Ellipse_sim_data$super_cluster

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
my_palette_ellipse <- c(
  "#A6D75B",  # spring green
  "#7BDFF2",  # sky cyan
  "#FFADAD",  # rose pink
  "#FFE066",  # lemon yellow
  "#B28DFF"   # soft lilac
)

my_palette_ellipse_truth <- c(
  "#FF6B6B",  # vivid coral red
  "#4D96FF",  # bright azure blue
  "#6A4C93",  # royal violet
  "#00B86B",  # fresh green teal
  "#FFA41B"   # mango orange
)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
### Function to generate the boundary of the elipses
generate_ellipse <- function(center, a, b, angle = 0, n = 100) {
  t <- seq(0, 2*pi, length.out = n)
  x <- a * cos(t)
  y <- b * sin(t)
  
  # Apply rotation
  x_rot <- x * cos(angle) - y * sin(angle)
  y_rot <- x * sin(angle) + y * cos(angle)
  
  data.frame(
    x = x_rot + center[1],
    y = y_rot + center[2]
  )
}
#-------------------------------------------------------------------------------
# Parameters for the two ellipses
# Ellipse 1: centered at (0, 0), a=3, b=2, rotated by pi/6
center1 <- c(0, 0)
a1 <- 3; b1 <- 2; angle1 <- pi/6

# Ellipse 2: centered at (8, 0), a=3, b=2, rotated by -pi/8
center2 <- c(8, 0)
a2 <- 3; b2 <- 2; angle2 <- -pi/8

ellipse1 <- generate_ellipse(center1, a1, b1, angle1)
ellipse2 <- generate_ellipse(center2, a2, b2, angle2)

#-------------------------------------------------------------------------------

ellipse1$contour_id <- 1
ellipse2$contour_id <- 2
boundary_df <- rbind(ellipse1, ellipse2)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

plot_points_paper <- function(coords, values, boundary_df, title = NULL, palette = my_palette_ellipse) {
  
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
    geom_point(data = pointsdata, aes(x = x, y = y, color = Value), size = 2.5) +
    
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
reference_vector <- seq(1, 5, 1)

# Function to reorder cluster assignments based on the reference vector
reorder_based_on_reference <- function(cluster_vector, reference_vector) {
  map <- match(reference_vector, cluster_vector)
  reordered_vector <- map[cluster_vector]
  return(reordered_vector)
}

#-------------------------------------------------------------------------------

# Plot the ground truth for the Ellipse example
true_clusters <- plot_points_paper(
  coords = as.matrix(Ellipse_sim_data[, c("x", "y")]),
  values = Ellipse_sim_data$super_cluster,
  boundary_df = boundary_df,
  palette = my_palette_ellipse_truth,
  title = bquote(bold("Truth, ") ~ italic("c = 5"))
) +
  # scale_y_reverse() +  # Optional flip Y-axis
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

ggsave("./New_Simulations/Ellipse/Ellipse_results/Ellipse_sim_true_clusters_MainPaper.png", 
       plot = true_clusters_with_label, 
       width = 5, height = 5, units = "in", dpi = 300, bg = "transparent")

#-------------------------------------------------------------------------------
##### DP-RST results for p=3 #####
load("./New_Simulations/Ellipse/Ellipse_results/Ellipse_sim_DP.RST_p3_30reps.RData")
output_p3 = Ellipse_sim_DP.RST_p3_reps

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
  
  DPM_accur_p3[j] = adjustedRandIndex(Ellipse_sim_data$super_cluster, obs_in_teams_vec_mode_based)
  clust_number_p3[j] = length(unique(obs_in_teams_vec_mode_based))
}

# Plot the accuracy over runs
plot(DPM_accur_p3, type = "l")

table(clust_number_p3)
# clust_number_p3
# 4  5 
# 28  2 

summary(DPM_accur_p3)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5160  0.5640  0.6244  0.6225  0.6882  0.7778

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
  coords = as.matrix(Ellipse_sim_data[, c("x", "y")]),
  values = DPM_partition_p3,
  boundary_df = boundary_df,
  palette = my_palette_ellipse,
  title = title_expr
) +
  scale_y_reverse() +
  theme(plot.title = element_text(size = 20, lineheight = 1.2))

DPM_results_p3


#-------------------------------------------------------------------------------
##### DP-RST results for p=10 #####
load("./New_Simulations/Ellipse/Ellipse_results/Ellipse_sim_DP.RST_p10_30reps.RData")
output_p10 = Ellipse_sim_DP.RST_p10_reps

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
  
  DPM_accur_p10[j] = adjustedRandIndex(Ellipse_sim_data$super_cluster, obs_in_teams_vec_mode_based)
  clust_number_p10[j] = length(unique(obs_in_teams_vec_mode_based))
}

# Plot the accuracy over runs
plot(DPM_accur_p10, type = "l")
summary(DPM_accur_p10)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7091  0.8162  0.8345  0.8281  0.8568  0.9053 

table(clust_number_p10)
# 5  6  7 
# 26  3  1

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
  coords = as.matrix(Ellipse_sim_data[, c("x", "y")]),
  values = DPM_partition_p10,
  boundary_df = boundary_df,
  palette = my_palette_ellipse,
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


ggsave("./New_Simulations/Ellipse/Ellipse_results/Ellipse_sim_DPM.png", 
       plot = DPM_results_with_label, 
       width = 11, height = 6, units = "in", dpi = 300, bg = "transparent")


ggsave("./New_Simulations/Ellipse/Ellipse_results/Ellipse_sim_DPM_MainPaper.png", 
       plot = DPM_results_with_label_MainPaper, 
       width = 11, height = 6, units = "in", dpi = 300, bg = "transparent")


#-------------------------------------------------------------------------------

### Make a histogram with the number of teams chosen over runs (for the Figure 1 in Supplementary Materials)
clust_number_p3 <- data.frame(clust_number = clust_number_p3)
clust_number_p3$clust_number <- factor(clust_number_p3$clust_number, levels = 1:5)

Ellipse_hist_plot_p3 <- ggplot(clust_number_p3, aes(x = clust_number)) +
  geom_bar(fill = "#FFADAD", color = "black") +
  theme_minimal() +
  labs(
    title = "Cluster Number Estimation in Ellipse Simulation
    p = 3",
    x = "Estimated number of clusters",
    y = "Frequency"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

Ellipse_hist_plot_p3


clust_number_p10 <- data.frame(clust_number = clust_number_p10)
clust_number_p10$clust_number <- factor(clust_number_p10$clust_number, levels = 1:7)

Ellipse_hist_plot_p10 <- ggplot(clust_number_p10, aes(x = clust_number)) +
  geom_bar(fill = "#B28DFF", color = "black") +
  theme_minimal() +
  labs(
    title = "Cluster Number Estimation in Ellipse Simulation
    p = 10",
    x = "Estimated number of clusters",
    y = "Frequency"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

Ellipse_hist_plot_p10

combined_plot_hist <- Ellipse_hist_plot_p3 + Ellipse_hist_plot_p10
combined_plot_hist

ggsave("./New_Simulations/Ellipse/Ellipse_hist_combined_plot.png", combined_plot_hist, width = 10, height = 4)

#-------------------------------------------------------------------------------
##### DENSITY PLOTS OF EACH PC PER CLUSTER #####

# Melt the data
df_melted <- reshape2::melt(
  Ellipse_sim_data,
  id.vars = "super_cluster",
  measure.vars = paste0("PC", 1:10),
  variable.name = "Variable",
  value.name = "Value"
)

# Make sure super_cluster is a factor
df_melted$super_cluster <- as.factor(df_melted$super_cluster)

# Compute means per group for vertical lines
means_df <- df_melted %>%
  group_by(super_cluster, Variable) %>%
  summarise(mean_value = mean(Value), .groups = "drop")

# Plot with super-cluster color mapping
density_plots_cluster <- ggplot(df_melted, aes(x = Value, fill = super_cluster)) +
  geom_density(alpha = 0.6, adjust = 1.5) +
  # Add vertical mean lines
  geom_vline(data = means_df, aes(xintercept = mean_value, color = super_cluster),
             linetype = "dashed", size = 0.5, show.legend = FALSE) +
  facet_grid(rows = vars(super_cluster), cols = vars(Variable), scales = "free") +
  scale_fill_manual(values = my_palette_ellipse_truth) +
  scale_color_manual(values = my_palette_ellipse_truth) +
  labs(
    # title = "Density Plots for Each Variable by Super-Cluster",
    x = "value",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    legend.position = "none",
    legend.key.size = unit(0.4, "cm"),
    strip.text = element_text(size = 16)
  )

density_plots_cluster

# Add label to the plot for Figure 3
density_with_label <- ggdraw() +
  draw_plot(density_plots_cluster) +
  draw_label("(B)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
density_with_label


# Save the combined plot
ggsave("./New_Simulations/Ellipse/Ellipse_results/Ellips_DensityFeatures.png", 
       plot = density_with_label, 
       width = 17, height = 7, units = "in", dpi = 300, bg = "transparent")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
##### BayesSpace results #####
load("./New_Simulations/Ellipse/Ellipse_results/Ellipse_sim_BayesSpace_p3_30reps.RData")

BayesSpace_acc_p3 <- c()
for (i in 1:30){
  BayesSpace_acc_p3[i] <- adjustedRandIndex(Ellipse_sim_data$super_cluster, Ellipse_sim_BayesSpace_p3_reps[[i]]$spatial.cluster)
  print(table(Ellipse_sim_BayesSpace_p3_reps[[i]]$spatial.cluster))
}

plot(BayesSpace_acc_p3, type = "l")
summary(BayesSpace_acc_p3)

BayesSpace_index_p3 <- which.max(BayesSpace_acc_p3)

BS_partition_p3 <- reorder_based_on_reference(Ellipse_sim_BayesSpace_p3_reps[[BayesSpace_index_p3]]$spatial.cluster, reference_vector)
se_BS_p3 <- sd(BayesSpace_acc_p3) / sqrt(length(BayesSpace_acc_p3))

BayesSpace_results_p3 <- plot_points_paper(
  coords = as.matrix(Ellipse_sim_data[, c("x", "y")]),
  values = BS_partition_p3,
  boundary_df = boundary_df,
  palette = my_palette_ellipse,
  title = bquote(bold("p = 3, ") ~ italic("ARI = ") ~ .(round(BayesSpace_acc_p3[BayesSpace_index_p3], 3)) ~ " (" ~ .(round(se_BS_p3, 3)) ~ ")")
) +
  # scale_y_reverse() +
  theme(plot.title = element_text(size = 20))
  
BayesSpace_results_p3

#-------------------------------------------------------------------------------
load("./New_Simulations/Ellipse/Ellipse_results/Ellipse_sim_BayesSpace_p10_30reps.RData")

BayesSpace_acc_p10 <- c()
for (i in 1:30){
  BayesSpace_acc_p10[i] <- adjustedRandIndex(Ellipse_sim_data$super_cluster, Ellipse_sim_BayesSpace_p10_reps[[i]]$spatial.cluster)
  print(table(Ellipse_sim_BayesSpace_p10_reps[[i]]$spatial.cluster))
}

plot(BayesSpace_acc_p10, type = "l")
summary(BayesSpace_acc_p10)

BayesSpace_index_p10 <- which.max(BayesSpace_acc_p10)

BS_partition_p10 <- reorder_based_on_reference(Ellipse_sim_BayesSpace_p10_reps[[BayesSpace_index_p10]]$spatial.cluster, reference_vector)
se_BS_p10 <- sd(BayesSpace_acc_p10) / sqrt(length(BayesSpace_acc_p10))

BayesSpace_results_p10 <- plot_points_paper(
  coords = as.matrix(Ellipse_sim_data[, c("x", "y")]),
  values = BS_partition_p10,
  boundary_df = boundary_df,
  palette = my_palette_ellipse,
  title = bquote(bold("p = 10, ") ~ italic("ARI = ") ~ .(round(BayesSpace_acc_p10[BayesSpace_index_p10], 3)) ~ " (" ~ .(round(se_BS_p10, 3)) ~ ")")
) +
  # scale_y_reverse() +
  theme(plot.title = element_text(size = 20))
BayesSpace_results_p10

#-------------------------------------------------------------------------------

# Title row with label and method name
title_row_BayesSpace <- plot_grid(
  ggdraw() + draw_label("(ii)", fontface = "bold", size = 20, hjust = 0, vjust = 1),
  ggdraw() + draw_label("BayesSpace", fontface = "bold", size = 24, hjust = 0.5, vjust = 1),
  ncol = 2,
  rel_widths = c(0.1, 0.9)
)

# Combined plots vertically
combined_body_BayesSpace <- BayesSpace_results_p3 / BayesSpace_results_p10 +
  plot_layout(heights = c(1, 1)) &
  theme(plot.margin = margin(t = 10, r = 5, b = 5, l = 5))

# Final layout
final_plot_BayesSpace <- plot_grid(
  title_row_BayesSpace,
  combined_body_BayesSpace,
  ncol = 1,
  rel_heights = c(0.1, 1.9)
)
final_plot_BayesSpace

# Save
ggsave(
  "./New_Simulations/Ellipse/Ellipse_results/Ellipse_sim_BayesSpace.png", 
  plot = final_plot_BayesSpace, 
  width = 5, height = 10, units = "in", dpi = 300, bg = "transparent"
)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
##### SC-MEB results #####
load("./New_Simulations/Ellipse/Ellipse_results/Ellipse_sim_SC.MEB_p3_30reps.RData")

SC.MEB_acc_p3 <- c()
for (i in 1:30){
  SC.MEB_acc_p3[i] <- adjustedRandIndex(Ellipse_sim_data$super_cluster, Ellipse_sim_SC.MEB_p3_reps[[i]][[1]])
}

plot(SC.MEB_acc_p3, type = "l")
summary(SC.MEB_acc_p3)

SC.MEB_index_p3 <- which.max(SC.MEB_acc_p3)

SC.MEB_partition_p3 <- reorder_based_on_reference(Ellipse_sim_SC.MEB_p3_reps[[SC.MEB_index_p3]][[1]], reference_vector)
se_SC.MEB_p3 <- sd(SC.MEB_acc_p3) / sqrt(length(SC.MEB_acc_p3))

SC.MEB_results_p3 <- plot_points_paper(
  coords = as.matrix(Ellipse_sim_data[, c("x", "y")]),
  values = SC.MEB_partition_p3,
  boundary_df = boundary_df,
  palette = my_palette_ellipse,
  title = bquote(bold("p = 3, ") ~ italic("ARI = ") ~ .(round(SC.MEB_acc_p3[SC.MEB_index_p3], 3)) ~ " (" ~ .(round(se_SC.MEB_p3, 3)) ~ ")")
) +
  # scale_y_reverse() +
  theme(plot.title = element_text(size = 20))

SC.MEB_results_p3

#-------------------------------------------------------------------------------
load("./New_Simulations/Ellipse/Ellipse_results/Ellipse_sim_SC.MEB_p10_30reps.RData")

SC.MEB_acc_p10 <- c()
for (i in 1:30){
  SC.MEB_acc_p10[i] <- adjustedRandIndex(Ellipse_sim_data$super_cluster, Ellipse_sim_SC.MEB_p10_reps[[i]][[1]])
}

plot(SC.MEB_acc_p10, type = "l")
summary(SC.MEB_acc_p10)

SC.MEB_index_p10 <- which.max(SC.MEB_acc_p10)

SC.MEB_partition_p10 <- reorder_based_on_reference(Ellipse_sim_SC.MEB_p10_reps[[SC.MEB_index_p10]][[1]], reference_vector)
se_SC.MEB_p10 <- sd(SC.MEB_acc_p10) / sqrt(length(SC.MEB_acc_p10))

SC.MEB_results_p10 <- plot_points_paper(
  coords = as.matrix(Ellipse_sim_data[, c("x", "y")]),
  values = SC.MEB_partition_p10,
  boundary_df = boundary_df,
  palette = my_palette_ellipse,
  title = bquote(bold("p = 10, ") ~ italic("ARI = ") ~ .(round(SC.MEB_acc_p10[SC.MEB_index_p10], 3)) ~ " (" ~ .(round(se_SC.MEB_p10, 3)) ~ ")")
) +
  # scale_y_reverse() +
  theme(plot.title = element_text(size = 20))

SC.MEB_results_p10

#-------------------------------------------------------------------------------

# Title row with label and method name
title_row_SC.MEB <- plot_grid(
  ggdraw() + draw_label("(iii)", fontface = "bold", size = 20, hjust = 0, vjust = 1),
  ggdraw() + draw_label("SC-MEB", fontface = "bold", size = 24, hjust = 0.5, vjust = 1),
  ncol = 2,
  rel_widths = c(0.1, 0.9)
)

# Combined plots vertically
combined_body_SC.MEB <- SC.MEB_results_p3 / SC.MEB_results_p10 +
  plot_layout(heights = c(1, 1)) &
  theme(plot.margin = margin(t = 10, r = 5, b = 5, l = 5))

# Final layout
final_plot_SC.MEB <- plot_grid(
  title_row_SC.MEB,
  combined_body_SC.MEB,
  ncol = 1,
  rel_heights = c(0.1, 1.9)
)
final_plot_SC.MEB

# Save
ggsave(
  "./New_Simulations/Ellipse/Ellipse_results/Ellipse_sim_SC.MEB.png", 
  plot = final_plot_SC.MEB, 
  width = 5, height = 10, units = "in", dpi = 300, bg = "transparent"
)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
##### DR-RC results #####
load("./New_Simulations/Ellipse/Ellipse_results/Ellipse_sim_DR.SC_p3_30reps.RData")

DR.SC_acc_p3 <- c()
for (i in 1:30){
  DR.SC_acc_p3[i] <- adjustedRandIndex(Ellipse_sim_data$super_cluster, Ellipse_sim_DR.SC_p3_reps[[i]][["Objdrsc"]][[1]][["cluster"]])
}

plot(DR.SC_acc_p3, type = "l")
summary(DR.SC_acc_p3)

DR.SC_index_p3 <- which.max(DR.SC_acc_p3)

DR.SC_partition_p3 <- reorder_based_on_reference(Ellipse_sim_DR.SC_p3_reps[[DR.SC_index_p3]][["Objdrsc"]][[1]][["cluster"]], reference_vector)
se_DR.SC_p3 <- sd(DR.SC_acc_p3) / sqrt(length(DR.SC_acc_p3))

DR.SC_results_p3 <- plot_points_paper(
  coords = as.matrix(Ellipse_sim_data[, c("x", "y")]),
  values = DR.SC_partition_p3,
  boundary_df = boundary_df,
  palette = my_palette_ellipse,
  title = bquote(bold("p = 3, ") ~ italic("ARI = ") ~ .(round(DR.SC_acc_p3[DR.SC_index_p3], 3)) ~ " (" ~ .(round(se_DR.SC_p3, 3)) ~ ")")
) +
  # scale_y_reverse() +
  theme(plot.title = element_text(size = 20))

DR.SC_results_p3

#-------------------------------------------------------------------------------
load("./New_Simulations/Ellipse/Ellipse_results/Ellipse_sim_DR.SC_p10_30reps.RData")

DR.SC_acc_p10 <- c()
for (i in 1:30){
  DR.SC_acc_p10[i] <- adjustedRandIndex(Ellipse_sim_data$super_cluster, Ellipse_sim_DR.SC_p10_reps[[i]][["Objdrsc"]][[1]][["cluster"]])
}

plot(DR.SC_acc_p10, type = "l")
summary(DR.SC_acc_p10)

DR.SC_index_p10 <- which.max(DR.SC_acc_p10)

DR.SC_partition_p10 <- reorder_based_on_reference(Ellipse_sim_DR.SC_p10_reps[[DR.SC_index_p10]][["Objdrsc"]][[1]][["cluster"]], reference_vector)
se_DR.SC_p10 <- sd(DR.SC_acc_p10) / sqrt(length(DR.SC_acc_p10))

DR.SC_results_p10 <- plot_points_paper(
  coords = as.matrix(Ellipse_sim_data[, c("x", "y")]),
  values = DR.SC_partition_p10,
  boundary_df = boundary_df,
  palette = my_palette_ellipse,
  title = bquote(bold("p = 10, ") ~ italic("ARI = ") ~ .(round(DR.SC_acc_p10[DR.SC_index_p10], 3)) ~ " (" ~ .(round(se_DR.SC_p10, 3)) ~ ")")
) +
  # scale_y_reverse() +
  theme(plot.title = element_text(size = 20))

DR.SC_results_p10

#-------------------------------------------------------------------------------

# Title row with label and method name
title_row_DR.SC <- plot_grid(
  ggdraw() + draw_label("(iv)", fontface = "bold", size = 20, hjust = 0, vjust = 1),
  ggdraw() + draw_label("DR-SC", fontface = "bold", size = 24, hjust = 0.5, vjust = 1),
  ncol = 2,
  rel_widths = c(0.1, 0.9)
)

# Combined plots vertically
combined_body_DR.SC <- DR.SC_results_p3 / DR.SC_results_p10 +
  plot_layout(heights = c(1, 1)) &
  theme(plot.margin = margin(t = 10, r = 5, b = 5, l = 5))

# Final layout
final_plot_DR.SC <- plot_grid(
  title_row_DR.SC,
  combined_body_DR.SC,
  ncol = 1,
  rel_heights = c(0.1, 1.9)
)
final_plot_DR.SC

# Save
ggsave(
  "./New_Simulations/Ellipse/Ellipse_results/Ellipse_sim_DR.SC.png", 
  plot = final_plot_DR.SC, 
  width = 5, height = 10, units = "in", dpi = 300, bg = "transparent"
)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
##### k-means results #####
load("./New_Simulations/Ellipse/Ellipse_results/Ellipse_sim_kmeans_3p_30reps.RData")

kmeans_acc_p3 <- c()
for (i in 1:30){
  kmeans_acc_p3[i] <- adjustedRandIndex(Ellipse_sim_data$super_cluster, Ellips_sim_kmeans_3p_reps[[i]])
}

plot(kmeans_acc_p3, type = "l")
summary(kmeans_acc_p3)

kmeans_index_p3 <- which.max(kmeans_acc_p3)

kmeans_partition_p3 <- reorder_based_on_reference(Ellips_sim_kmeans_3p_reps[[kmeans_index_p3]], reference_vector)
se_kmeans_p3 <- sd(kmeans_acc_p3) / sqrt(length(kmeans_acc_p3))

kmeans_results_p3 <- plot_points_paper(
  coords = as.matrix(Ellipse_sim_data[, c("x", "y")]),
  values = kmeans_partition_p3,
  boundary_df = boundary_df,
  palette = my_palette_ellipse,
  title = bquote(bold("p = 3, ") ~ italic("ARI = ") ~ .(round(kmeans_acc_p3[kmeans_index_p3], 3)) ~ " (" ~ .(round(se_kmeans_p3, 3)) ~ ")")
) +
  # scale_y_reverse() +
  theme(plot.title = element_text(size = 20))

kmeans_results_p3

#-------------------------------------------------------------------------------
load("./New_Simulations/Ellipse/Ellipse_results/Ellipse_sim_kmeans_10p_30reps.RData")

kmeans_acc_p10 <- c()
for (i in 1:30){
  kmeans_acc_p10[i] <- adjustedRandIndex(Ellipse_sim_data$super_cluster, Ellips_sim_kmeans_10p_reps[[i]])
}

plot(kmeans_acc_p10, type = "l")
summary(kmeans_acc_p10)

kmeans_index_p10 <- which.max(kmeans_acc_p10)

kmeans_partition_p10 <- reorder_based_on_reference(Ellips_sim_kmeans_10p_reps[[kmeans_index_p10]], reference_vector)
se_kmeans_p10 <- sd(kmeans_acc_p10) / sqrt(length(kmeans_acc_p10))

kmeans_results_p10 <- plot_points_paper(
  coords = as.matrix(Ellipse_sim_data[, c("x", "y")]),
  values = kmeans_partition_p10,
  boundary_df = boundary_df,
  palette = my_palette_ellipse,
  title = bquote(bold("p = 10, ") ~ italic("ARI = ") ~ .(round(max(kmeans_acc_p10[kmeans_index_p10]), 3)) ~ " (" ~ .(round(se_kmeans_p10, 3)) ~ ")")
) +
  # scale_y_reverse() +
  theme(plot.title = element_text(size = 20))

kmeans_results_p10

#-------------------------------------------------------------------------------

# Title row with label and method name
title_row_kmeans <- plot_grid(
  ggdraw() + draw_label("(v)", fontface = "bold", size = 20, hjust = 0, vjust = 1),
  ggdraw() + draw_label("k-means", fontface = "bold", size = 24, hjust = 0.5, vjust = 1),
  ncol = 2,
  rel_widths = c(0.1, 0.9)
)

# Combine plots vertically
combined_body_kmeans <- kmeans_results_p3 / kmeans_results_p10 +
  plot_layout(heights = c(1, 1)) &
  theme(plot.margin = margin(t = 10, r = 5, b = 5, l = 5))

# Final stacked layout
final_plot_kmeans <- plot_grid(
  title_row_kmeans,
  combined_body_kmeans,
  ncol = 1,
  rel_heights = c(0.1, 1.9)
)
final_plot_kmeans

# Save the plot
ggsave(
  "./New_Simulations/Ellipse/Ellipse_results/Ellipse_sim_kmeans.png",
  plot = final_plot_kmeans,
  width = 5, height = 10, units = "in", dpi = 300, bg = "transparent"
)
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
               kmeans_acc_p3, kmeans_acc_p10),
  Method = factor(rep(c("DP-RST", "Bayes\nSpace", "DR.SC", "SC.MEB", "k-means"), each = 60),
                  levels = c("DP-RST", "Bayes\nSpace", "DR.SC", "SC.MEB", "k-means")),
  P = factor(rep(rep(c("p = 3", "p = 10"), each = 30), times = 5),
             levels = c("p = 3", "p = 10"))
)


# Create the box plot
boxplot <- ggplot(data, aes(x = Method, y = Accuracy, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8), aes(group = interaction(Method, P))) +
  facet_wrap(~P, ncol = 2) +  # Optional: if you want facets. If you prefer grouped bars, remove.
  scale_fill_manual(values = c(
    "DP-RST" = "#B28DFF",
    "Bayes\nSpace" = "#FFE066",
    "DR.SC" = "#903030FF",
    "SC.MEB" = "#FFFFFFFF",
    "k-means" = "#7BDFF2"
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
  draw_label("(C)", x = 0.1, y = 0.98, size = 20, hjust = 0, vjust = 1)
boxplot_with_label

ggsave("./New_Simulations/Ellipse/Ellipse_results/Ellipse_sim_boxplot.png", 
       plot = boxplot_with_label, 
       width = 12, height = 5, units = "in", dpi = 300, bg = "transparent")
