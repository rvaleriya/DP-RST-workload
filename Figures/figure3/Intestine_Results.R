##### LOAD NECESSARY LIBRARIES #####
library(DP.RST)
library(mclust)
library(ggplot2)
library(readr)
library(dplyr)
library(gridExtra)

#-------------------------------------------------------------------------------
# SETUP & GLOBAL VARIABLES
#-------------------------------------------------------------------------------
setwd("~/Desktop/DP-RST-workload/Datasets Results/Gut")
OUTPUT_DIR <- "~/Desktop/DP-RST-workload/Figures/figure3/"
set.seed(7303)

# Palettes
my_palette1 <- c("#FCDD23FF", "#AD98F1FF", "#D856A7FF", "#F8B100FF", "#48439BFF")
reference_vector <- c(1:5)

# Standard Theme for all plots (Removes axes, sets font size)
common_theme <- theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(size = 35, hjust = 0.5, face = "bold")
  )

#-------------------------------------------------------------------------------
# DATA LOADING (GROUND TRUTH)
#-------------------------------------------------------------------------------
gut_df_wt_muscle <- readRDS("./gut_df_wt_muscle.rds")
load("swiss_roll_wt_muscle_boundary.RData")

# Extract and Scale Ground Truth
Y_sample <- gut_df_wt_muscle[, 1:3]
loc <- gut_df_wt_muscle[, c("x", "y")]
bnd <- boundary
z <- gut_df_wt_muscle$z

coords <- apply(loc, 2, scale)
# Create Master DataFrame for Merging
df_subset <- data.frame(loc, z)
df_subset$spots <- rownames(df_subset)

# Scale Boundary
bnd_scaled <- list(
  x = (bnd$x - mean(loc[,1]))/sd(loc[,1]),
  y = (bnd$y - mean(loc[,2]))/sd(loc[,2])
)

#-------------------------------------------------------------------------------
# HELPER FUNCTIONS
#-------------------------------------------------------------------------------

# Optimized Reorder Function (Vectorized)
reorder_based_on_reference <- function(cluster_vector, reference_vector) {
  # match() is significantly faster than a for loop
  return(match(cluster_vector, reference_vector))
}

# Generic Plot Function (No Boundary)
groups_plot_wB <- function(coords, group_assign, point_size = 2.5) {
  group_data <- as.data.frame(cbind(coords, group_assign))
  
  ggplot(group_data, aes(x = group_data[,1], y = group_data[,2], colour = as.factor(group_assign))) + 
    geom_point(size = point_size) +
    labs(x = 'Coordinate X', y = 'Coordinate Y') + 
    ggtitle('Groups of Observations') 
}

# Master Processing Function for Competitor Methods
process_and_plot_method <- function(config, ground_truth_df, ref_vec, palette, theme_obj, out_path) {
  
  message(paste("Processing:", config$name))
  
  # 1. Load Data based on extension
  if (grepl("\\.rds$", config$file, ignore.case = TRUE)) {
    raw_data <- readRDS(config$file)
  } else {
    raw_data <- read_csv(config$file, show_col_types = FALSE)
  }
  
  # 2. Extract Coordinates and Clusters based on config mapping
  # Note: logic handles cases where x and y are switched in raw files
  # If config says x_col = "y" and y_col = "x", it handles the switch.
  
  # Handle special case for SpaGCN which implies row index matching or specific columns
  # Creating a temp dataframe to standardize
  temp_df <- data.frame(
    x_raw = raw_data[[config$x_col]],
    y_raw = raw_data[[config$y_col]],
    cluster_raw = raw_data[[config$z_col]]
  )
  
  # 3. Add offset if needed (e.g. converting 0-indexed python results)
  if (!is.null(config$offset)) {
    temp_df$cluster_raw <- temp_df$cluster_raw + config$offset
  }
  
  # 4. Prepare for Merge
  # We map the raw x/y to the column names 'x' and 'y' expected by df_subset
  # IMPORTANT: If the method flipped X and Y, we map them correctly here to match ground truth
  
  merge_prep <- data.frame(
    x = temp_df$x_raw, 
    y = temp_df$y_raw,
    method_z = temp_df$cluster_raw
  )
  
  # Special merge for SpaGCN if needed (uses spots/rownames), otherwise x/y merge
  merged_df <- merge(merge_prep, ground_truth_df, by = c("x", "y"), all = TRUE)
  
  # 5. Calculate ARI
  ari <- adjustedRandIndex(merged_df$z, merged_df$method_z)
  
  # 6. Reorder Clusters for Coloring
  final_partition <- reorder_based_on_reference(merged_df$method_z, ref_vec)
  
  # 7. Plotting
  # Use the standard coordinate columns for plotting
  p <- groups_plot_wB(merged_df[, c("x", "y")], final_partition, point_size = 2.5) + # SIZE 2 ENFORCED
    scale_color_manual(values = palette) +
    theme_obj +
    ggtitle(bquote(bold(.(config$name)) * ", c = 5, " * bold("ARI = ") * bold(.(format(round(ari, 3), nsmall = 3)))))
  
  # 8. Save
  ggsave(filename = paste0(out_path, config$save_name), plot = p, width = 8, height = 8, bg = "transparent")
  
  return(ari)
}
#-------------------------------------------------------------------------------
# PART 1: PROCESS DP-RST (THE MAIN METHOD - WITH BOUNDARY)
#-------------------------------------------------------------------------------
# Keeps separate logic because it loads an RData environment and uses boundary plotting

load("./Results/Gut_DP.RST_FromNewBastPT_p3_Version2_OutputOnly.RData")
output <- Gut_DP.RST_FromNewBastPT_Version2

mode_based_partition <- partition(DP.RST_output = output, method = "mode_based", batch_size = 100)

# Calculate Accuracy
X <- table(sequence(length(mode_based_partition$groups_partition)), mode_based_partition$groups_partition)
Z <- table(sequence(length(mode_based_partition$teams_partition)), mode_based_partition$teams_partition)
obs_in_teams <- X %*% Z
obs_in_teams_vec <- obs_in_teams %*% sort(unique(mode_based_partition$teams_partition))
accur_teams <- adjustedRandIndex(df_subset$z, obs_in_teams_vec)

# Reorder
DPM_partition <- reorder_based_on_reference(obs_in_teams_vec, reference_vector)

# Plot WITH BOUNDARY (Using groups_plot from DP.RST package)
# Note: Providing point_size = 2
DPM.RST_plot <- groups_plot(coords, DPM_partition, bnd_scaled, point_size = 2.5) +
  scale_color_manual(values = my_palette1) +
  common_theme +
  ggtitle(bquote(bold("DP-RST") * ", c = 5, " * bold("ARI = ") * bold(.(format(round(accur_teams, 3), nsmall = 3)))))

ggsave(filename = paste0(OUTPUT_DIR, "DPM-RST_plot.png"), plot = DPM.RST_plot, width = 8, height = 8, bg = "transparent")

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
ggsave(filename = paste0(OUTPUT_DIR, "hist_MCMC_clusters.png"),
       plot = cluster_count_hist, width = 8, height = 6,
       units = "in", dpi = 300, bg = "transparent")

#-------------------------------------------------------------------------------
# PART 2: CONFIGURE ALL COMPETING METHODS
#-------------------------------------------------------------------------------
# Define how to read each file. 
# x_col/y_col: Which column in the RESULT file maps to the x/y in GROUND TRUTH.
# If a method swaps coords, put the 'y' column name in x_col.

methods_config <- list(
  # --- 4.2 BayesSpace ---
  list(name = "BayesSpace", file = "./Results/BayesSpace_Gut_clustering_results.rds", 
       x_col = "imagerow", y_col = "imagecol", z_col = "cluster_q5_pc3", save_name = "BS_plot.png"),
  
  # --- 4.3 SC-MEB ---
  list(name = "SC-MEB", file = "./Results/SC-MEB_Gut_clustering_results.rds", 
       x_col = "x", y_col = "y", z_col = "pc3_k5", save_name = "SC.MEB_plot.png"),
  
  # --- 4.4 DR-SC ---
  list(name = "DR-SC", file = "./Results/DR-SC_Gut_clustering_results.rds", 
       x_col = "imagerow", y_col = "imagecol", z_col = "cluster_q5_pc3", save_name = "DR.SC_plot.png"),
  
  # --- 4.5 K-Means ---
  list(name = "k-means", file = "./Results/kmeans_Gut_clustering_results.rds", 
       x_col = "x", y_col = "y", z_col = "pc3_k5", save_name = "kmeans_plot.png"),
  
  # --- 4.6 BASS ---
  list(name = "BASS", file = "./Results/BASS_Intestine_3PCs.rds", 
       x_col = "x", y_col = "y", z_col = "z", save_name = "BASS_plot.png"),
  
  # --- 4.7 SEDR (Swapped x/y) ---
  list(name = "SEDR", file = "./Results/SEDR_Gut.csv", 
       x_col = "y", y_col = "x", z_col = "label_pca3_k5", save_name = "sedr_plot.png"),
  
  # --- 4.8 GraphST (Updated path & var, Swapped x/y) ---
  list(name = "GraphST", file = "./Results/GraphST_Gut.csv", 
       x_col = "y", y_col = "x", z_col = "mclust_3PCs", save_name = "graphst_plot.png"),
  
  # --- 4.9 STAGATE (Updated path & var, Swapped x/y) ---
  list(name = "STAGATE", file = "./Results/STAGATE_Gut.csv", 
       x_col = "y", y_col = "x", z_col = "mclust", save_name = "stagate_plot.png"),
  
  # --- 4.10 SpaGCN (Without) (Updated path & var, Swapped x/y, Offset) ---
  list(name = "SpaGCN (w/o)", file = "./Results/SpaGCN_Gut.csv", 
       x_col = "y", y_col = "x", z_col = "refined_pred_3PCs_k5", 
       offset = 1, save_name = "spagcn_without_plot.png"),
  
  # --- 4.10 SpaGCN (With) (Updated path & var, Swapped x/y, Offset) ---
  list(name = "SpaGCN (w/)", file = "./Results/SpaGCN_Gut_histology.csv", 
       x_col = "y", y_col = "x", z_col = "refined_pred_3PCs_k5", 
       offset = 1, save_name = "spagcn_with_plot.png"),
  
  # --- 4.11 stLearn (With) (Updated path & var, Swapped x/y, Offset) ---
  list(name = "stLearn (w/)", file = "./Results/stLearn_Gut_with_image.csv", 
       x_col = "y", y_col = "x", z_col = "label_3PCs_k5", 
       offset = 1, save_name = "stlearn_with_plot.png"),
  
  # --- 4.11 stLearn (Without) (Updated path & var, Swapped x/y, Offset) ---
  list(name = "stLearn (w/o)", file = "./Results/stLearn_Gut_without_image.csv", 
       x_col = "y", y_col = "x", z_col = "label_3PCs_k5", 
       offset = 1, save_name = "stlearn_without_plot.png"),
  
  # --- 4.12 SpaMask (Updated path, Offset, X=X / Y=Y as per snippet code) ---
  list(name = "SpaMask", file = "./Results/SpaMask_Gut_3PCs.csv", 
       x_col = "x", y_col = "y", z_col = "kmeans", 
       offset = 1, save_name = "SpaMask_plot.png"),
  
  # --- 4.13 STAMarker (New Method, Swapped x/y, Offset) ---
  list(name = "STAMarker", file = "./Results/STAMarker_Gut_joint.csv", 
       x_col = "y", y_col = "x", z_col = "STAMarker_label", 
       offset = 1, save_name = "STAMarker_plot.png")
)
#-------------------------------------------------------------------------------
# PART 3: EXECUTE LOOP
#-------------------------------------------------------------------------------

for (config in methods_config) {
  tryCatch({
    process_and_plot_method(
      config = config,
      ground_truth_df = df_subset,
      ref_vec = reference_vector,
      palette = my_palette1,
      theme_obj = common_theme,
      out_path = OUTPUT_DIR
    )
  }, error = function(e) {
    message(paste("Error processing", config$name, ":", e$message))
  })
}

message("All plots generated successfully.")
