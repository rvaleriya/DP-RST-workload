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
  "Sensitivity/Gut_DP.RST_Sensitivity_alpha005_OutputOnly.RData",
  "Sensitivity/Gut_DP.RST_Sensitivity_alpha01_OutputOnly.RData",
  "Sensitivity/Gut_DP.RST_Sensitivity_alpha02_OutputOnly.RData",
  "Results/Gut_DP.RST_FromNewBastPT_p3_Version2_OutputOnly.RData", # alpha = 0.5
  "Sensitivity/Gut_DP.RST_Sensitivity_alpha1_OutputOnly.RData"
)

# Containers for final "Alpha Overview" figure
hist_list <- list()     # per-alpha histograms of # refined clusters
mode_list <- list()     # per-alpha mode refined partition plots
alpha_vals <- numeric() # actual alpha read from each output

# ----------------------------------------------------------------------
# Loop: process each run (per-alpha)
# ----------------------------------------------------------------------
for (fp in file_paths) {
  # Load output object safely
  obj_name <- load(fp)     # returns the object name
  output <- get(obj_name)  # the DP.RST output list
  
  # ---------------- Diagnostics: marginal likelihood, team counts -----
  iters_vec <- seq_along(output$marginal_likelihood_out)
  marg_lik <- ggplot(
    data.frame(x = iters_vec,
               y = output$marginal_likelihood_out),
    aes(x, y)
  ) +
    geom_line() +
    ggtitle("Marginal Likelihood post BurnIn and Thinning") +
    xlab("Iteration") +
    ylab("Value") +
    theme_classic()
  
  teams_freq <- as.data.frame(table(output$j_teams_out))
  colnames(teams_freq) <- c("teams", "Freq")
  moves_freq_plot <- ggplot(teams_freq, aes(x = teams, y = Freq)) +
    geom_bar(stat = "identity", color = "black", fill = "grey70") +
    labs(title = "Frequency of # of refined clusters", x = "# of Teams", y = "Frequency") +
    theme_classic()
  
  # ---------------- ARI across MCMC (spatial vs refined) --------------
  iters <- seq_along(output$marginal_likelihood_out)
  accur_spatial <- numeric(length(iters))
  accur_refined <- numeric(length(iters))
  
  for (i in iters) {
    # spatial partition ARI
    Gut_df <- data.frame(res_memb = output$cluster_out[i, ], z = z)
    Gut_df <- subset(Gut_df, !is.na(z))
    accur_spatial[i] <- adjustedRandIndex(Gut_df$res_memb, Gut_df$z)
    
    # refined partition ARI (map groups->teams at iteration i)
    X <- table(seq_along(output$cluster_out[i, ]), output$cluster_out[i, ])
    Z <- table(seq_along(output$teams_out[[i]]),  output$teams_out[[i]])
    obs_in_teams <- X %*% Z
    obs_in_teams_vec <- obs_in_teams %*% sort(unique(output$teams_out[[i]]))
    
    Gut_df2 <- data.frame(res_memb = as.vector(obs_in_teams_vec), z = z)
    Gut_df2 <- subset(Gut_df2, !is.na(z))
    accur_refined[i] <- adjustedRandIndex(Gut_df2$res_memb, Gut_df2$z)
  }
  
  acc_spatial <- ggplot(data.frame(x = iters, y = accur_spatial), aes(x, y)) +
    geom_line() +
    ggtitle("Accuracy of Spatial Partition") +
    xlab("Iteration") +
    ylab("ARI") +
    theme_classic()
  
  acc_refined <- ggplot(data.frame(x = iters, y = accur_refined), aes(x, y)) +
    geom_line() +
    ggtitle("Accuracy of Refined Partition") +
    xlab("Iteration") +
    ylab("ARI") +
    theme_classic()
  
  # ---------------- Best partitions (MODE only) -----------------------
  mode_based_partition <- partition(output, method = "mode_based", batch_size = 100)
  
  # ARIs for final MODE partitions
  ari_groups_mode <- adjustedRandIndex(mode_based_partition$groups_partition, z)
  ari_teams_mode  <- adjustedRandIndex(mode_based_partition$teams_partition[mode_based_partition$groups_partition], z)
  
  # ---------------- Plots: true / groups(mode) / teams(mode) ----------
  true_clusters <- groups_plot(coords, as.factor(z), bnd_scaled) +
    theme(legend.position = "none") +
    ggtitle(paste0("True cluster memberships, k = ", length(unique(z))))
  
  groups_plotted_mode <- groups_plot(coords, mode_based_partition$groups_partition, bnd_scaled) +
    theme(legend.position = "none") +
    ggtitle(paste0("DP-RST Spatial Partition (mode), k = ",
                   length(unique(mode_based_partition$groups_partition)),
                   ", ARI = ", round(ari_groups_mode, 3)))
  
  teams_plotted_mode <- teams_plot(coords,
                                   mode_based_partition$teams_partition,
                                   mode_based_partition$groups_partition,
                                   bnd_scaled) +
    theme(legend.position = "none") +
    ggtitle(paste0("DP-RST Refined Partition (mode), k = ",
                   length(unique(mode_based_partition$teams_partition)),
                   ", ARI = ", round(ari_teams_mode, 3)))
  
  # ---------------- Title / per-alpha save ----------------------------
  sigmasq_mu <- output$hyperpar$sigmasq_mu
  nu         <- output$hyperpar$nu
  alpha      <- output$hyperpar$alpha
  M          <- output$hyperpar$M
  k_max      <- output$hyperpar$k_max
  j_max      <- output$hyperpar$j_max
  PCs        <- length(output$init_val$mu_teams[[1]])
  
  file_name <- sub("\\.RData$", "", basename(fp))
  
  title_text <- sprintf(
    "Gut, DP-RST Sensitivity: PCs = %d, alpha = %.3f\nsigmasq_mu = %.3f, nu = %d, M = %d, k_max = %d, j_max = %d",
    PCs, alpha, sigmasq_mu, nu, M, k_max, j_max
  )
  
  plots <- ggarrange(
    true_clusters,
    teams_plotted_mode,
    groups_plotted_mode,
    moves_freq_plot,
    acc_refined,
    marg_lik,
    ncol = 3, nrow = 2
  ) +
    patchwork::plot_annotation(
      title = title_text,
      theme = theme(plot.title = element_text(size = 14, color = "darkblue", face = "bold"))
    )
  
  save_path <- paste0("Sensitivity/", file_name, ".pdf")
  ggsave(filename = save_path, plot = plots, width = 8, height = 6, scale = 2, device = cairo_pdf)
  
  cat("\n✅ Finished:", file_name, "\nSaved:", save_path, "\n")
}