# ==============================================================================
# UMAP PLOTTING SCRIPT FOR MOUSE HEALING COLON RESULTS WITH 3 PCs
# ==============================================================================

##### 1. LIBRARIES & SETUP #####
library(DP.RST)
library(mclust)
library(ggplot2)
library(readr)
library(dplyr)
library(umap)       # specific for UMAP calculation
library(gridExtra)

# Set base paths based on your provided directory structure
# Assuming script is run from the project root or you adjust this path:
BASE_DATA_DIR <- "~/Desktop/DP-RST-workload/Datasets Results/Gut"
FIG_OUTPUT_DIR <- "~/Desktop/DP-RST-workload/Figures/umap"

# Create output directory if it doesn't exist
if (!dir.exists(FIG_OUTPUT_DIR)) {
  dir.create(FIG_OUTPUT_DIR, recursive = TRUE)
}

setwd(BASE_DATA_DIR)
set.seed(7303)

##### 2. LOAD & PROCESS GROUND TRUTH (COMPUTE UMAP ONCE) #####

# Load pre-processed data
gut_df <- readRDS("./gut_df_wt_muscle.rds")

# Extract PCs for UMAP calculation
# The user code used cols 1:3 for PCs
Y_sample <- gut_df[, 1:3]
loc <- gut_df[, c("x", "y")]
true_labels <- gut_df$z

# --- CALCULATION: UMAP ---
# 
# We calculate UMAP on the PCs once. All methods will be projected onto this map.
cat("Calculating UMAP on ground truth PCs...\n")
umap_config <- umap.defaults
umap_config$n_neighbors <- 30 # Adjustable based on local vs global structure preference
umap_res <- umap(Y_sample, config = umap_config)

# Create a Master Dataframe with UMAP coords and Spatial Coords
master_df <- data.frame(
  x = loc$x,
  y = loc$y,
  UMAP1 = umap_res$layout[,1],
  UMAP2 = umap_res$layout[,2],
  True_Labels = true_labels
)

##### 3. HELPER FUNCTIONS (THE OPTIMIZATION) #####

# A. Color Palette
my_palette1 <- c("#FCDD23FF", "#AD98F1FF", "#D856A7FF", "#F8B100FF", "#48439BFF")
reference_vector <- c(1:5)

# B. Reordering Function
reorder_clusters <- function(cluster_vector, ref) {
  reordered <- numeric(length(cluster_vector))
  for (i in seq_along(cluster_vector)) {
    position <- match(cluster_vector[i], ref)
    reordered[i] <- position
  }
  return(reordered)
}

# C. Unified Plotting Function
# This function takes the method's labels, calculates ARI, merges with UMAP, and saves.
generate_method_umap <- function(method_name, cluster_labels, coords_df, file_suffix) {
  
  # 1. Merge labels with master UMAP data
  # We assume 'coords_df' has columns x, y to match master_df
  # If coords_df is NULL, we assume cluster_labels are already in same order as master_df
  
  if (!is.null(coords_df)) {
    # Combine input labels with their x/y
    temp_df <- data.frame(x = coords_df$x, y = coords_df$y, raw_cluster = cluster_labels)
    
    # Merge with Master (which has UMAP coords)
    plot_df <- merge(master_df, temp_df, by = c("x", "y"))
  } else {
    plot_df <- master_df
    plot_df$raw_cluster <- cluster_labels
  }
  
  # 2. Calculate ARI
  ari_val <- adjustedRandIndex(plot_df$True_Labels, plot_df$raw_cluster)
  
  # 3. Reorder for Color consistency
  # Note: This assumes raw_cluster are integers 1-5 or similar. 
  # If they are strings, this might need adjustment.
  plot_df$final_cluster <- reorder_clusters(plot_df$raw_cluster, reference_vector)
  plot_df$final_cluster <- as.factor(plot_df$final_cluster)
  
  # 4. Plotting
  p <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = final_cluster)) +
    geom_point(size = 3, alpha = 0.8) + # Adjusted size/alpha for UMAP density
    scale_color_manual(values = my_palette1) +
    theme_minimal(base_size = 20) +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(size = 28, hjust = 0.5, face = "bold")
    ) +
    ggtitle(bquote(bold(.(method_name)) * ", ARI = " * bold(.(format(round(ari_val, 3), nsmall = 3)))))
  
  # 5. Save
  filename <- file.path(FIG_OUTPUT_DIR, paste0(file_suffix, "_UMAP.png"))
  ggsave(filename = filename, plot = p, width = 8, height = 8, bg = "transparent")
  
  cat(paste("Saved:", method_name, "| ARI:", round(ari_val, 3), "\n"))
  return(p)
}


##### 4. GENERATE PLOTS FOR ALL METHODS #####

# --- 4.1 DP-RST ---
load("./Results/Gut_DP.RST_FromNewBastPT_p3_Version2_OutputOnly.RData")
output <- Gut_DP.RST_FromNewBastPT_Version2
# Helper to get teams
partition_res <- partition(DP.RST_output = output, method = "mode_based", batch_size = 100)
X <- table(sequence(length(partition_res$groups_partition)), partition_res$groups_partition)
Z <- table(sequence(length(partition_res$teams_partition)), partition_res$teams_partition)
obs_in_teams_vec <- (X %*% Z) %*% sort(unique(partition_res$teams_partition))

# DP-RST is already in same order as gut_df, so coords_df is NULL
generate_method_umap("DP-RST", as.numeric(obs_in_teams_vec), NULL, "DP-RST")


# --- 4.2 BayesSpace ---
bs_res <- readRDS("./Results/BayesSpace_Gut_clustering_results.rds")
bs_coords <- data.frame(x = bs_res$imagerow, y = bs_res$imagecol)
generate_method_umap("BayesSpace", bs_res$cluster_q5_pc3, bs_coords, "BayesSpace")


# --- 4.3 SC-MEB ---
sc_res <- readRDS("./Results/SC-MEB_Gut_clustering_results.rds")
sc_coords <- data.frame(x = sc_res$x, y = sc_res$y)
generate_method_umap("SC-MEB", sc_res$pc3_k5, sc_coords, "SC-MEB")


# --- 4.4 DR-SC ---
dr_res <- readRDS("./Results/DR-SC_Gut_clustering_results.rds")
dr_coords <- data.frame(x = dr_res$imagerow, y = dr_res$imagecol)
generate_method_umap("DR-SC", dr_res$cluster_q5_pc3, dr_coords, "DR-SC")


# --- 4.5 K-Means ---
km_res <- readRDS("./Results/kmeans_Gut_clustering_results.rds")
km_coords <- data.frame(x = km_res$x, y = km_res$y)
generate_method_umap("k-means", km_res$pc3_k5, km_coords, "kmeans")


# --- 4.6 BASS ---
bass_res <- readRDS("./Results/BASS_Intestine_3PCs.rds")
bass_coords <- data.frame(x = bass_res$x, y = bass_res$y)
generate_method_umap("BASS", bass_res$z, bass_coords, "BASS")


# --- 4.7 SEDR ---
sedr_res <- read_csv("./Results/SEDR_Gut.csv", show_col_types = FALSE)
# Note: swapping x/y as per your original script logic
sedr_coords <- data.frame(x = sedr_res$y, y = sedr_res$x) 
generate_method_umap("SEDR", sedr_res$label_pca3_k5, sedr_coords, "SEDR")


# --- 4.8 GraphST ---
gst_res <- read_csv("./Results/GraphST_Gut.csv", show_col_types = FALSE)
gst_coords <- data.frame(x = gst_res$y, y = gst_res$x)
generate_method_umap("GraphST", gst_res$mclust_3PCs, gst_coords, "GraphST")


# --- 4.9 STAGATE ---
sta_res <- read_csv("./Results/STAGATE_Gut.csv", show_col_types = FALSE)
sta_coords <- data.frame(x = sta_res$y, y = sta_res$x)
generate_method_umap("STAGATE", sta_res$mclust, sta_coords, "STAGATE")


# --- 4.10 SpaGCN (With & Without) ---
spa_res <- read_csv("./Results/SpaGCN_Gut.csv", show_col_types = FALSE)
# SWAP x and y here:
spa_coords <- data.frame(x = spa_res$y, y = spa_res$x)
# Without Histology
generate_method_umap("SpaGCN (w/o)", (spa_res$refined_pred_3PCs_k5 + 1), spa_coords, "SpaGCN_without")

spa_res <- read_csv("./Results/SpaGCN_Gut_histology.csv", show_col_types = FALSE)
# SWAP x and y here:
spa_coords <- data.frame(x = spa_res$y, y = spa_res$x)
# With Histology (adding 1 as per your script)
generate_method_umap("SpaGCN (w/)", (spa_res$refined_pred_3PCs_k5 + 1), spa_coords, "SpaGCN_with")


# --- 4.11 stLearn (With & Without) ---
stl_with <- read_csv("./Results/stLearn_Gut_with_image.csv", show_col_types = FALSE)
stl_no <- read_csv("./Results/stLearn_Gut_without_image.csv", show_col_types = FALSE)

# Careful: stLearn inputs were separate files
stl_w_coords <- data.frame(x = stl_with$y, y = stl_with$x)
generate_method_umap("stLearn (w/)", (stl_with$label_3PCs_k5 + 1), stl_w_coords, "stLearn_with")

stl_no_coords <- data.frame(x = stl_no$y, y = stl_no$x)
generate_method_umap("stLearn (w/o)", (stl_no$label_3PCs_k5 + 1), stl_no_coords, "stLearn_without")


# --- 4.12 SpaMask ---
mask_res <- read_csv("./Results/SpaMask_Gut_3PCs.csv", show_col_types = FALSE)
# Logic from your script: Rename x=y and y=x
mask_coords <- data.frame(x = mask_res$x, y = mask_res$y)
generate_method_umap("SpaMask", (mask_res$kmeans + 1), mask_coords, "SpaMask")

cat("\nAll UMAP plots generated in:", FIG_OUTPUT_DIR, "\n")

# --- 4.13 STAMarker ---
stam_res <- read_csv("./Results/STAMarker_Gut_joint.csv", show_col_types = FALSE)
# Logic from your script: Rename x=y and y=x
stam_coords <- data.frame(x = stam_res$y, y = stam_res$x)
generate_method_umap("STAMarker", (stam_res$STAMarker_label + 1), stam_coords, "STAMarker")

cat("\nAll UMAP plots generated in:", FIG_OUTPUT_DIR, "\n")
