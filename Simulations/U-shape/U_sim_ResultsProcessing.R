# Load necessary libraries
library(DP.RST)
library(ggplot2)
library(patchwork)
library(mclust)
library(dplyr)
library(tidyr)
library(magick)
library(cowplot)

#-------------------------------------------------------------------------------
# Set working directory to the root of the R project
setwd("~/Desktop/DP-RST-workload")

# set seed
set.seed(9362)

# For debugging purposes
options(error=recover)

### Load Data ###
load("~/UshapeSim_coords_boundary1.RData")

#-------------------------------------------------------------------------------

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

# my_palette <- paletteer::paletteer_d("rcartocolor::Fall")

my_palette <- c(
  "#698C50",   # slightly more green than original
  "#95CC5EFF",   # brighter yellow-green
  "#F3E96B",   # deeper golden yellow vs pale cream
  "#EDBB8AFF",   # more saturated orange-beige
  "#C83E4D",   # brighter and redder coral
  "#8E1B10"    # deeper burnt orange
)

# my_palette1 <- paletteer::paletteer_d("ggsci::legacy_tron")

my_palette1 <- c(
  "#FF3C00",   # punchier red-orange
  "#F7D02C",   # vivid goldenrod
  "#AA1D2F",   # saturated slate blue
  "#3EC4ED",   # deeper cyan (less white)
  "#5E4FA2",   # more vivid tangerine
  "#2CA02C"   # brighter grass green
)

plot_points_paper <- function(coords, values, bnd_scaled, title = NULL, angle = 0) {
  require(ggplot2)
  
  # Rotation matrix
  rotation_matrix <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), nrow = 2)
  
  # Rotate the points coordinates
  rotated_coords <- as.matrix(coords) %*% rotation_matrix
  
  # Convert boundary coordinates to a matrix and rotate
  bnd_matrix <- cbind(bnd_scaled$x, bnd_scaled$y)
  rotated_bnd <- bnd_matrix %*% rotation_matrix
  
  # Preparing data for points
  pointsdata <- data.frame(
    x = rotated_coords[, 1],
    y = rotated_coords[, 2],
    Value = values  # Converting scaled values to a factor for discrete coloring
  )
  
  # Preparing data for boundary
  bnd_data <- data.frame(
    x = rotated_bnd[, 1],
    y = rotated_bnd[, 2]
  )
  
  # Create the ggplot object with boundary and points
  plot <- ggplot() +
    geom_point(aes(x = x, y = y, color = as.factor(Value)), data = pointsdata, size = 4) +
    geom_path(data = bnd_data, aes(x = x, y = y), color = "black", size = 1) +  # Adding boundary
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
    )
  
  return(plot)
}


# Apply a scaling transformation to the coordinates
scaled_coords <- coords
scaled_coords[, 1] <- coords[, 1] * 2  # Scale the x-coordinates by a factor (e.g., 2)

# Apply the same scaling transformation to the boundary coordinates
scaled_bnd <- bnd_scaled
scaled_bnd$x <- bnd_scaled$x * 2  # Scale the x-coordinates of the boundary

reference_vector <- c(1, 2, 3, 4, 5, 6)

reorder_based_on_reference <- function(cluster_vector, reference_vector) {
  map <- match(reference_vector, cluster_vector)
  reordered_vector <- map[cluster_vector]
  return(reordered_vector)
}

library(grid)
true_clusters <- plot_points_paper(scaled_coords, sim_data$cluster, scaled_bnd, angle = pi/4) +
  scale_color_manual(values = rev(my_palette), 
                     labels = c(expression(mu == -2), 
                                expression(mu == -1.2), 
                                expression(mu == -0.4),
                                expression(mu == 0.4),
                                expression(mu == 1.2),
                                expression(mu == 2))) +
  ggtitle(expression(bold("Truth"))) +
  theme(plot.title = element_text(size = 20),
        legend.position = "left",
        legend.title = element_blank(),
        legend.text = element_text(size = 24))

true_clusters

true_clusters_with_label <- ggdraw() +
  draw_plot(true_clusters) +
  draw_label("(A)", x = 0, y = 0.98, size = 20, hjust = 0, vjust = 1)
true_clusters_with_label


ggsave("./New_Simulations/Ushape_Extra/U_sim_true_clusters.png",
       plot = true_clusters_with_label,
       width = 6.5, height = 5, units = "in", dpi = 300, bg = "transparent")


##### DPM results #####
load("/U_sim_DPM_30reps.RData")

### Choose the iteration from each run ###
DPM_spatial <- list()
DPM_refinment <- list()

for (m in 1:30){
  ##### CHOOSE THE ITERATION #####
  groups_assign = U_sim_DPM_reps[[m]][["cluster_out"]]
  teams_assign = U_sim_DPM_reps[[m]][["teams_out"]]
  teams_number = U_sim_DPM_reps[[m]][["j_teams_out"]]
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
  
  # Fins the minimum Frobenius norm index
  min_index <- mode_indices[which.min(norm)]
  
  DPM_spatial[[m]] <- groups_assign[min_index, ]
  DPM_refinment[[m]] <- teams_assign[[min_index]]
}


### Calculate the accuracy for each run ###
DPM_accur = c()
k_numbers = c()

for (j in 1:30) {
  X <- table(sequence(length(DPM_spatial[[j]])), DPM_spatial[[j]])
  Z <- table(sequence(length(DPM_refinment[[j]])), DPM_refinment[[j]])
  
  obs_in_teams <- X %*% Z
  
  obs_in_teams_vec <- obs_in_teams %*% sort(unique(DPM_refinment[[j]]))
  DPM_accur[j] = adjustedRandIndex(sim_data$cluster, obs_in_teams_vec)
  # first value is 0.3944193
  k_numbers[j] <- length(unique(DPM_refinment[[j]]))
}

plot(DPM_accur, type = "l")
summary(DPM_accur)

## Get the cluster assignments of the highest accuracy
indices <- which(k_numbers == 6)
max_accuracy_index <- indices[which.max(DPM_accur[indices])] #30

X <- table(sequence(length(DPM_spatial[[30]])), DPM_spatial[[30]])
Z <- table(sequence(length(DPM_refinment[[30]])), DPM_refinment[[30]])
obs_in_teams <- X %*% Z
DPM_partition <- obs_in_teams %*% sort(unique(DPM_refinment[[30]]))

DPM_partition <- reorder_based_on_reference(DPM_partition, reference_vector)
se_DPM <- sd(DPM_accur) / sqrt(length(DPM_accur))

DPM_results <- plot_points_paper(scaled_coords, DPM_partition, scaled_bnd, angle = pi/4) +
  scale_color_manual(values = my_palette1) +
  ggtitle(bquote(bold("DP-RST, ") ~ italic("ARI = ") ~ .(round(max(DPM_accur), 3)) ~ " (" ~ .(round(se_DPM, 3)) ~ ")")) + 
  theme(plot.title = element_text(size = 20),
        legend.position = "none")
DPM_results

DPM_results_with_label <- ggdraw() +
  draw_plot(DPM_results) +
  draw_label("(i)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
DPM_results_with_label


DPM_results_MainPaper <- plot_points_paper(scaled_coords, DPM_partition, scaled_bnd, angle = pi/4) +
  scale_color_manual(values = my_palette1) +
  ggtitle(
    bquote(bold("p = 3, ") * italic("ARI") * " = " * .(round(max(DPM_accur), 3)) * 
             " (" * .(round(se_DPM, 3)) * ")")
  ) +
  theme(plot.title = element_text(size = 12),
        legend.position = "none")

DPM_results_MainPaper

combined_plot_DPM_results_MainPaper <- DPM_results_MainPaper +
  plot_annotation(title = bquote(bold("DP-RST"))) &
  theme(
    plot.title = element_text(size = 20, hjust = 0.5)  # Center title
  )
combined_plot_DPM_results_MainPaper

DPM_results_with_label_MainPaper <- ggdraw() +
  draw_plot(combined_plot_DPM_results_MainPaper) +
  draw_label("(B)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
DPM_results_with_label_MainPaper

ggsave("./New_Simulations/Ushape_Extra/U_sim_DPM.png", 
       plot = DPM_results_with_label, 
       width = 5, height = 5, units = "in", dpi = 300, bg = "transparent")

ggsave("./New_Simulations/Ushape_Extra/U_sim_DPM_MainPaper.png", 
       plot = DPM_results_with_label_MainPaper, 
       width = 5, height = 5, units = "in", dpi = 300, bg = "transparent")


# Convert to a data frame
clust_number <- data.frame(clust_number = k_numbers)

# Create the histogram
U_hist_plot <- ggplot(clust_number, aes(x = clust_number)) +
  geom_histogram(binwidth = 0.5, fill = "#778868FF", color = "black") +
  theme_minimal() +
  labs(title = "Cluster Number Estimation in U-shape Simulation", 
       x = "Estimated number of clusters", 
       y = "Frequency") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  xlim(4, 8)  # Extend the x-axis from 1 to 5

ggsave("./New_Simulations/Ushape_Extra/U_hist_plot.png", U_hist_plot, width = 6, height = 4)


##### BayesSpace results #####
load("/U_sim_BayesSpace_30reps.RData")

BayesSpace_acc <- c()
clust_6 <- c()

for (i in 1:30){
  BayesSpace_acc[i] <- adjustedRandIndex(sim_data$cluster, U_sim_BayesSpace_reps[[i]]$spatial.cluster)
  if (length(unique(U_sim_BayesSpace_reps[[i]]$spatial.cluster)) == 6) {
    clust_6[i] <- i
  }
}
# first value is 0.2166038
plot(BayesSpace_acc, type = "l")
summary(BayesSpace_acc[clust_6])

which.max(BayesSpace_acc[clust_6]) #19

BS_partition <- reorder_based_on_reference(U_sim_BayesSpace_reps[[19]]$spatial.cluster, reference_vector)
se_BS <- sd(BayesSpace_acc) / sqrt(length(BayesSpace_acc))

BayesSpace_results <- plot_points_paper(scaled_coords, 
                                        BS_partition, 
                                        scaled_bnd, angle = pi/4) +
  scale_color_manual(values = my_palette1) +
  ggtitle(bquote(bold("BS, ") ~ italic("ARI = ") ~ .(round(max(BayesSpace_acc), 3)) ~ " (" ~ .(round(se_BS, 3)) ~ ")"))+
  theme(plot.title = element_text(size = 20))
BayesSpace_results

BayesSpace_results_with_label <- ggdraw() +
  draw_plot(BayesSpace_results) +
  draw_label("(ii)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
BayesSpace_results_with_label

ggsave("./New_Simulations/Ushape_Extra/U_sim_BayesSpace.png", 
       plot = BayesSpace_results_with_label, 
       width = 5, height = 5, units = "in", dpi = 300, bg = "transparent")


##### SC-MEB results #####
load("/U_sim_SC.MEB_30reps.RData")

SC.MEB_acc <- c()

for (i in 1:30){
  SC.MEB_acc[i] <- adjustedRandIndex(sim_data$cluster, U_sim_SC.MEB_reps[[i]][[37]])
}
# first is 0.2875501
plot(SC.MEB_acc, type = "l")
summary(SC.MEB_acc)

which.max(SC.MEB_acc) #1

SC.MEB_partition <- reorder_based_on_reference(U_sim_SC.MEB_reps[[1]][[37]], reference_vector)
se_SC.MEB <- sd(SC.MEB_acc) / sqrt(length(SC.MEB_acc))

SC.MEB_results <- plot_points_paper(scaled_coords, SC.MEB_partition, 
                                    scaled_bnd, angle = pi/4) +
  scale_color_manual(values = my_palette1) +
  ggtitle(bquote(bold("SC-MEB, ") ~ italic("ARI = ") ~ .(round(max(SC.MEB_acc), 3)) ~ " (" ~ .(round(se_SC.MEB, 3)) ~ ")")) +
  theme(plot.title = element_text(size = 20))
SC.MEB_results

SC.MEB_results_with_label <- ggdraw() +
  draw_plot(SC.MEB_results) +
  draw_label("(iii)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
SC.MEB_results_with_label

ggsave("./New_Simulations/Ushape_Extra/U_sim_SC.MEB.png", 
       plot = SC.MEB_results_with_label, 
       width = 5, height = 5, units = "in", dpi = 300, bg = "transparent")


##### DR-RC results #####
load("/U_sim_DR.SC_30reps.RData")

DR.SC_acc <- c()

for (i in 1:30){
  DR.SC_acc[i] <- adjustedRandIndex(sim_data$cluster, 
                                    U_sim_DR.SC_reps[[i]][["Objdrsc"]][[1]][["cluster"]])
}
# first is 0.209915
plot(DR.SC_acc, type = "l")
summary(DR.SC_acc)

which.max(DR.SC_acc) #1

DR.SC_partition <- reorder_based_on_reference(U_sim_DR.SC_reps[[i]][["Objdrsc"]][[1]][["cluster"]], reference_vector)
se_DR.SC <- sd(DR.SC_acc) / sqrt(length(DR.SC_acc))

DR.SC_results <- plot_points_paper(scaled_coords, 
                                   DR.SC_partition, 
                                   scaled_bnd, angle = pi/4) +
  scale_color_manual(values = my_palette1) +
  ggtitle(bquote(bold("DR-SC, ") ~ italic("ARI = ") ~ .(round(max(DR.SC_acc), 3)) ~ " (" ~ .(round(se_DR.SC, 3)) ~ ")")) +
  theme(plot.title = element_text(size = 20))
DR.SC_results

DR.SC_results_with_label <- ggdraw() +
  draw_plot(DR.SC_results) +
  draw_label("(iv)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
DR.SC_results_with_label

ggsave("./New_Simulations/Ushape_Extra/U_sim_DR.SC.png", 
       plot = DR.SC_results_with_label, 
       width = 5, height = 5, units = "in", dpi = 300, bg = "transparent")



##### k-means results #####
load("/U_sim_kmeans_30reps.RData")

kmeans_acc <- c()

for (i in 1:30){
  kmeans_acc[i] <- adjustedRandIndex(sim_data$cluster, U_sim_kmeans_reps[[i]])
}
# first is 0.2575392
plot(kmeans_acc, type = "l")
summary(kmeans_acc)

which.max(kmeans_acc) #5

kmeans_partition <- reorder_based_on_reference(U_sim_kmeans_reps[[5]], reference_vector)
se_kmeans <- sd(kmeans_acc) / sqrt(length(kmeans_acc))

kmeans_results <- plot_points_paper(scaled_coords, kmeans_partition, 
                                    scaled_bnd, angle = pi/4) +
  scale_color_manual(values = my_palette1) +
  ggtitle(bquote(bold("k-means, ") ~ italic("ARI = ") ~ .(round(max(kmeans_acc), 3)) ~ " (" ~ .(round(se_kmeans, 3)) ~ ")")) +
  theme(plot.title = element_text(size = 20))
kmeans_results

kmeans_results_with_label <- ggdraw() +
  draw_plot(kmeans_results) +
  draw_label("(v)", x = 0, y = 0.99, size = 20, hjust = 0, vjust = 1)
kmeans_results_with_label

ggsave("./New_Simulations/Ushape_Extra/U_sim_kmeans.png", 
       plot = kmeans_results_with_label, 
       width = 5, height = 5, units = "in", dpi = 300, bg = "transparent")


##### Create box-plot for accuracy #####
# Combine the vectors into a data frame
# Combine the vectors into a data frame
data <- data.frame(
  Accuracy = c(DPM_accur, BayesSpace_acc, DR.SC_acc, SC.MEB_acc, kmeans_acc),
  Method = factor(rep(c("DP-RST", "BayesSpace", "DR.SC", "SC.MEB", "k-means"), each = 30), 
                  levels = c("DP-RST", "BayesSpace", "DR.SC", "SC.MEB", "k-means"))
)


# Create the box plot
boxplot <- ggplot(data, aes(x = Method, y = Accuracy, fill = Method)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = c("DP-RST" = "#95CC5EFF", "BayesSpace" = "#EDBB8AFF", 
                               "DR.SC" = "#CA562CFF", "SC.MEB" = "#3D5941FF", 
                               "k-means" = "#F6EDBDFF")) +
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

boxplot_with_label_MainPaper <- ggdraw() +
  draw_label("Compare ARI", x = 0.5, y = 0.98, size = 20, fontface = "bold", hjust = 0.5, vjust = 1) +
  draw_plot(boxplot) +
  draw_label("(C)", x = 0.1, y = 0.98, size = 20, hjust = 0, vjust = 1)
boxplot_with_label_MainPaper


ggsave("./New_Simulations/Ushape_Extra/U_sim_boxplot.png", 
       plot = boxplot_with_label, 
       width = 8.5, height = 5, units = "in", dpi = 300, bg = "transparent")

ggsave("./New_Simulations/Ushape_Extra/U_sim_boxplot_MainPaper.png", 
       plot = boxplot_with_label_MainPaper, 
       width = 8.5, height = 5, units = "in", dpi = 300, bg = "transparent")


##### Read the plots from .png files and arrange them for the paper #####

# Read the images
image_true <- image_read("./New_Simulations/Ushape_Extra/U_sim_true_clusters.png")
image_DPM <- image_read("./New_Simulations/Ushape_Extra/U_sim_DPM.png")
image_BS <- image_read("./New_Simulations/Ushape_Extra/U_sim_BayesSpace.png")
image_DR.SC <- image_read("./New_Simulations/Ushape_Extra/U_sim_DR.SC.png")
image_SC.MEB <- image_read("./New_Simulations/Ushape_Extra/U_sim_SC.MEB.png")
image_kmeans <- image_read("./New_Simulations/Ushape_Extra/U_sim_kmeans.png")
image_boxplot <- image_read("./New_Simulations/Ushape_Extra/U_sim_boxplot.png")

# Create a collage
# Append images side by side (horizontal append)
row1 <- image_append(c(image_true, image_DPM, image_boxplot))
row2 <- image_append(c(image_BS, image_DR.SC, image_SC.MEB, image_kmeans))

# Stack the rows (vertical append)
U_sim_collage <- image_append(c(row1, row2), stack = TRUE)

# Save the collage
image_write(U_sim_collage, "./New_Simulations/Ushape_Extra/U_sim_collage.png")
