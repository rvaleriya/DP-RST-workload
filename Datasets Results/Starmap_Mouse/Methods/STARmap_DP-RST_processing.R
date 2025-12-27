# ==================== DP-RST post-processing: STARmap (3PCs & 10PCs) ====================

library(DP.RST)
library(ggplot2)
library(dplyr)
library(patchwork)
library(mclust)

# ---- Paths -------------------------------------------------------------------
base_dir    <- "~/Desktop/DP-RST-workload/Datasets Results/Starmap_Mouse"
results_dir <- file.path(base_dir, "Results")
setwd(base_dir)

# ---- Load coords and boundary ------------------------------------------------
pcs_full <- read.csv(file.path(base_dir, "starmap_pcs_coords_labels.csv"), check.names = FALSE)
bnd_all  <- read.csv(file.path(base_dir, "starmap_manual_boundary.csv"))

coords_scaled <- scale(as.matrix(pcs_full[, c("x", "y")]))

bnd_scaled <- list(
  x = (bnd_all$x - mean(pcs_full$x)) / sd(pcs_full$x),
  y = (bnd_all$y - mean(pcs_full$y)) / sd(pcs_full$y)
)
bnd_scaled_df <- data.frame(x = bnd_scaled$x, y = bnd_scaled$y)

z_true <- as.integer(factor(pcs_full$ground_truth))

# ---- Helpers -----------------------------------------------------------------
.safe_load_result <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (exists("result", e)) e$result else e$output
}

.to_int <- function(x) as.integer(factor(as.vector(x)))
.ari3 <- function(a, b) round(adjustedRandIndex(a, b), 3)

plot_points <- function(coords, labels, bnd_df, title) {
  df <- data.frame(x = coords[,1], y = coords[,2], lab = factor(labels))
  ggplot(df, aes(x, y, color = lab)) +
    geom_point(size = 0.7) +
    geom_path(data = bnd_df, aes(x, y), inherit.aes = FALSE, linewidth = 0.35) +
    coord_equal() +
    theme_classic(base_size = 10) +
    labs(title = title, x = NULL, y = NULL) +
    theme(legend.position = "none")
}

process_one <- function(pc) {
  cat("\n=== Processing", pc, "PCs ===\n")
  rdata_path <- file.path(results_dir, sprintf("DP-RST-files/Full_Result_DP.RST_%dPCs.RData", pc))
  output <- .safe_load_result(rdata_path)
  
  n <- nrow(coords_scaled)
  init_lab <- .to_int(output[["init_val"]][["cluster"]][,1])
  
  mode_part <- partition(output, method = "mode_based", batch_size = 100)
  groups_lab <- .to_int(mode_part$groups_partition)
  teams_lab  <- .to_int(mode_part$teams_partition)[groups_lab]
  
  ari_init  <- .ari3(z_true, init_lab)
  ari_group <- .ari3(z_true, groups_lab)
  ari_team  <- .ari3(z_true, teams_lab)
  
  k_true  <- length(unique(z_true))
  k_init  <- length(unique(init_lab))
  k_group <- length(unique(groups_lab))
  k_team  <- length(unique(teams_lab))
  
  p_true   <- plot_points(coords_scaled, z_true,      bnd_scaled_df, sprintf("True labels | K = %d", k_true))
  p_init   <- plot_points(coords_scaled, init_lab,    bnd_scaled_df, sprintf("Initialization — ARI = %.3f | K = %d", ari_init,  k_init))
  p_groups <- plot_points(coords_scaled, groups_lab,  bnd_scaled_df, sprintf("Spatial partition | Groups (mode) — ARI = %.3f | K = %d", ari_group, k_group))
  p_teams  <- plot_points(coords_scaled, teams_lab,   bnd_scaled_df, sprintf("Refined partition | Teams (mode) — ARI = %.3f | K = %d", ari_team,  k_team))
  
  panel <- (p_true | p_teams) / (p_init | p_groups) +
    plot_annotation(title = sprintf("STARmap — DP-RST Post-Processing (%d PCs)", pc),
                    theme = theme(plot.title = element_text(size = 12, face = "bold")))
  
  out_pdf <- file.path(results_dir, sprintf("STARmap_DP-RST_Postproc_%dPCs.pdf", pc))
  ggsave(out_pdf, plot = panel, device = cairo_pdf, width = 10, height = 8, dpi = 300)
  cat("Saved:", out_pdf, "\n")
}

# ---- Run for 3 PCs and 10 PCs ------------------------------------------------
for (pc in c(3, 10)) process_one(pc)
