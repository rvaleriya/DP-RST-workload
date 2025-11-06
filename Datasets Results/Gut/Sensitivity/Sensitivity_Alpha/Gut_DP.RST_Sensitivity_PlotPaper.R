# ======================================================================
# DP-RST Sensitivity Analysis for Gut Dataset — Batch Postprocessing
# + Final "Alpha Overview" Figure (histograms + mode refined partitions)
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
# Sensitivity file paths (alpha variants)
# ----------------------------------------------------------------------

file_paths <- c(
  "Sensitivity/Sensitivity_Alpha/PC3/Gut_DP.RST_Sensitivity_alpha005_OutputOnly.RData",
  "Sensitivity/Sensitivity_Alpha/PC3/Gut_DP.RST_Sensitivity_alpha01_OutputOnly.RData",
  "Sensitivity/Sensitivity_Alpha/PC3/Gut_DP.RST_Sensitivity_alpha02_OutputOnly.RData",
  "Results/Gut_DP.RST_FromNewBastPT_p3_Version2_OutputOnly.RData", # alpha = 0.5
  "Sensitivity/Sensitivity_Alpha/PC3/Gut_DP.RST_Sensitivity_alpha1_OutputOnly.RData"
)

# file_paths <- c(
#   "Sensitivity/Sensitivity_Alpha/PC10/Gut_DPRST_PT_alpha0.05_OutputOnly.RData",
#   "Sensitivity/Sensitivity_Alpha/PC10/Gut_DPRST_PT_alpha0.10_OutputOnly.RData",
#   "Sensitivity/Sensitivity_Alpha/PC10/Gut_DPRST_PT_alpha0.20_OutputOnly.RData",
#   # "Sensitivity/Sensitivity_Alpha/PC10/Gut_DPRST_PT_alpha0.30_OutputOnly.RData",
#   "Results/Gut_DP.RST_FromNewBastPT_p10_Version2_OutputOnly.RData", # alpha = 0.5
#   # "Sensitivity/Sensitivity_Alpha/PC10/Gut_DPRST_PT_alpha0.70_OutputOnly.RData",
#   "Sensitivity/Sensitivity_Alpha/PC10/Gut_DPRST_PT_alpha1.00_OutputOnly.RData"
# )

# Containers for final "Alpha Overview" figure
hist_list <- list()     # per-alpha histograms of # refined clusters
mode_list <- list()     # per-alpha mode refined partition plots
alpha_vals <- numeric() # actual alpha read from each output

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

alpha_grid <- c(0.05, 0.1, 0.2, 0.5, 1)
counter = 0
# ----------------------------------------------------------------------
# Loop: process each run (per-alpha)
# ----------------------------------------------------------------------
for (fp in file_paths) {
  
  counter  = counter + 1
  # Load output object safely
  obj_name <- load(fp)     # returns the object name
  output <- get(obj_name)  # the DP.RST output list
  
  # ---------------- Diagnostics: team counts -----
  
  teams_freq <- as.data.frame(table(output$j_teams_out))
  
  hist_list[[counter]] <- ggplot(teams_freq, aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity", color = "#C497A7", fill = "#F8E6EC") +
    labs(title = bquote(alpha == .(alpha_grid[counter])),
         x = "# of Refined Clusters", y = "Frequency") +
    theme_classic()

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
    ggtitle(paste("c =", length(unique(mode_based_partition$teams_partition)),
                  ", ARI =", round(ari_teams_mode, 3)))
}

# ---------- Arrange All Plots in a Grid ----------

hist_list <- Filter(Negate(is.null), hist_list)
mode_list <- Filter(Negate(is.null), mode_list)

# Combine the three lists into one vector of plots
all_plots <- c(hist_list, mode_list)

# Arrange in a grid:
# - Number of columns = number of alpha values (i.e. length of SR_sim_SensAlpha)
# - 3 rows: row 1 (histograms), row 2 (mode-based plots), row 3 (frobenius plots)
grid.arrange(grobs = all_plots, ncol = 5, nrow = 2)

# Create the arranged grob
final_plot <- arrangeGrob(grobs = all_plots, ncol = 5, nrow = 2)

# Save as PDF in the "Sensitivity" folder
ggsave("Sensitivity/Sensitivity_Alpha/PC3/Sensetivity_Alpha_PC3.pdf", final_plot, width = 15, height = 6)
