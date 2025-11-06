# ======================================================================
# DP-RST Sensitivity Analysis for Gut Dataset — Batch Postprocessing
# ======================================================================

library(ggplot2)
library(igraph)
library(ggpubr)
library(patchwork)
library(mclust)
library(DP.RST)
library(gridExtra)   # for grid.arrange / arrangeGrob
library(paletteer)   # optional color palettes
library(grid)

setwd("~/Desktop/DP-RST-workload/Datasets Results/Gut")

# ----------------------------------------------------------------------
# Load base data (coordinates, boundary, true labels)
# ----------------------------------------------------------------------
load("swiss_roll_wt_muscle_finaltouches1.RData")
load("swiss_roll_wt_muscle_boundary.RData")

Y_sample <- swiss_roll_wt_muscle_finaltouches1[, 1:3]
loc      <- swiss_roll_wt_muscle_finaltouches1[, 4:5]
z        <- swiss_roll_wt_muscle_finaltouches1$z
z_man    <- swiss_roll_wt_muscle_finaltouches1$z_man
n        <- nrow(Y_sample)

# Standardize
coords <- apply(loc, 2, scale)
Y_std  <- apply(Y_sample, 2, scale)

bnd <- boundary
bnd_scaled <- list(
  x = (bnd$x - mean(loc[, 1])) / sd(loc[, 1]),
  y = (bnd$y - mean(loc[, 2])) / sd(loc[, 2])
)

# ----------------------------------------------------------------------
# Sensitivity file paths (K variants)
# ----------------------------------------------------------------------

file_paths <- c(
  "Sensitivity/Sensitivity_K/PC3/Gut_DPRST_PT_k05_OutputOnly.RData",
  "Sensitivity/Sensitivity_K/PC3/Gut_DPRST_PT_k10_OutputOnly.RData",
  "Sensitivity/Sensitivity_K/PC3/Gut_DPRST_PT_k15_OutputOnly.RData",
  "Results/Gut_DP.RST_FromNewBastPT_p3_Version2_OutputOnly.RData", # K = 20
  "Sensitivity/Sensitivity_K/PC3/Gut_DPRST_PT_k25_OutputOnly.RData"
)

# file_paths <- c(
#   "Sensitivity/Sensitivity_K/PC10/Gut_DPRST_PT_k05_OutputOnly.RData",
#   "Sensitivity/Sensitivity_K/PC10/Gut_DPRST_PT_k10_OutputOnly.RData",
#   "Sensitivity/Sensitivity_K/PC10/Gut_DPRST_PT_k15_OutputOnly.RData",
#   "Results/Gut_DP.RST_FromNewBastPT_p10_Version2_OutputOnly.RData", # K = 20
#   "Sensitivity/Sensitivity_K/PC10/Gut_DPRST_PT_k25_OutputOnly.RData"
# )

# Containers for final "K Overview" figure
mode_list <- list()     # per-K mode refined partition plots
K_vals <- numeric() # actual K read from each output

#-------------------------------------------------------------------------------
##### FUNCTION TO REORDER CLUSTER LABELS BASED ON THE REFERENCE VECTOR #####

# Create a reference vector with values from 1 to 9
reference_vector <- c(1:8)

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

# Define a custom color palette for plotting
my_palette1 <- c(
  "#FCDD23FF",  # bright yellow
  "#AD98F1FF",  # light lavender
  "#D856A7FF",  # magenta pink
  "#F8B100FF",  # orange-gold
  "#48439BFF",  # deep indigo
  "#70C1B3FF",  # soft teal
  "#F28E8EFF",  # warm coral
  "#7F4F24FF"   # earthy muted brown
)

#-------------------------------------------------------------------------------

# Initialize lists to store plots for each simulation output
hist_list <- list()
mode_list <- list()

K_grid <- c(5, 10, 15, 20, 25)
counter = 0
# ----------------------------------------------------------------------
# Loop: process each run (per-K)
# ----------------------------------------------------------------------
for (fp in file_paths) {
  
  counter  = counter + 1
  # Load output object safely
  obj_name <- load(fp)     # returns the object name
  output <- get(obj_name)  # the DP.RST output list
  
  # ---------------- Best partitions (MODE only) -----------------------
  mode_based_partition <- partition(output, method = "mode_based", batch_size = 100)
  
  obs_in_teams_vec <- mode_based_partition$teams_partition[mode_based_partition$groups_partition]
  
  # ARIs for final MODE partitions
  ari_teams_mode  <- adjustedRandIndex(obs_in_teams_vec, z)
  
  # ---------------- Plots: true / groups(mode) / teams(mode) ----------
  # Reorder the partition based on a reference vector
  DPM_partition <- reorder_based_on_reference(obs_in_teams_vec, reference_vector)
  
  ##### Plot the crude spatial partition for Figure 1 #####
  mode_list[[counter]] <- groups_plot(coords, DPM_partition, bnd_scaled, point_size = 0.5) +
    scale_color_manual(values = my_palette1) +
    geom_point(size = 0.5) +
    theme_classic() +  # Remove grey background
    theme(legend.position = "none",
          plot.title = element_text(face = "italic", size = 12)) +  # Italicize title
    ggtitle(paste("k =", length(unique(output$init_val$cluster[,1])),
                  "c =", length(unique(mode_based_partition$teams_partition)),
                  ", ARI =", round(ari_teams_mode, 3)))
}

# ---------- Arrange All Plots in a Grid ----------

mode_list <- Filter(Negate(is.null), mode_list)

# Combine the three lists into one vector of plots
all_plots <- c(mode_list)

# Arrange in a grid:
# - Number of columns = number of K values (i.e. length of SR_sim_SensK)
# - 3 rows: row 1 (histograms), row 2 (mode-based plots), row 3 (frobenius plots)
grid.arrange(grobs = all_plots, ncol = 5, nrow = 1)

# Create the arranged grob
final_plot <- arrangeGrob(grobs = all_plots, ncol = 5, nrow = 1)

# Save as PDF in the "Sensitivity" folder
ggsave("Sensitivity/Sensitivity_K/PC3/Sensitivity_K_PC3.pdf", final_plot, 
       width = 15, height = 3)

