# ==================== osmFISH: Density plots of PC1–PC10 by ground_truth ====================
# Each facet = one PC; curves colored by true label (ground_truth)

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)   # for color palettes

# ---- Paths -----------------------------------------------------------------------------
base_dir <- "~/Desktop/DP-RST-workload/Datasets Results/osmFISH"
setwd(base_dir)

# ---- Load data -------------------------------------------------------------------------
pcs_full <- read.csv(file.path(base_dir, "osmfish_pcs_coords_labels.csv"), check.names = FALSE)

# ---- Verify ground_truth ----------------------------------------------------------------
if (!"ground_truth" %in% names(pcs_full)) {
  stop("Column 'ground_truth' not found in osmfish_pcs_coords_labels.csv")
}

# ---- Select first 10 PCs ---------------------------------------------------------------
pc_cols_all <- grep("^PC\\d+$", names(pcs_full), value = TRUE)
if (length(pc_cols_all) < 10) stop("Fewer than 10 PCs found.")
pc_cols <- pc_cols_all[1:10]

# ---- Z-score PCs (column-wise) ---------------------------------------------------------
pcs_scaled <- scale(as.matrix(pcs_full[, pc_cols]))
colnames(pcs_scaled) <- pc_cols

# ---- Build long data frame: ground_truth + PCs -----------------------------------------
plot_df <- cbind.data.frame(ground_truth = pcs_full$ground_truth,
                            as.data.frame(pcs_scaled)) |>
  mutate(ground_truth = factor(ground_truth)) |>
  pivot_longer(cols = all_of(pc_cols), names_to = "PC", values_to = "value") |>
  mutate(PC = factor(PC, levels = pc_cols))  # keep PC1..PC10 order

# ---- Plot: density per label, faceted by PC (2×5) --------------------------------------
p_den <- ggplot(plot_df, aes(x = value, fill = ground_truth)) +
  geom_density(alpha = 0.45, linewidth = 0.25, position = "identity", adjust = 1) +
  facet_wrap(~ PC, ncol = 5, scales = "free_x") +
  scale_fill_viridis_d(option = "C") +
  labs(
    title = "osmFISH — Density of PC1–PC10 by Ground Truth (z-scored PCs)",
    x = "PC value (z-score)", y = "Density", fill = "Ground truth"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 8),
    strip.text = element_text(size = 10)
  )

# ---- Save PDF --------------------------------------------------------------------------
out_pdf <- file.path(base_dir, "osmFISH_PC1-10_Density_byGroundTruth.pdf")
ggsave(out_pdf, plot = p_den, device = cairo_pdf, width = 12, height = 6.5, dpi = 300)
message("Saved: ", out_pdf)
