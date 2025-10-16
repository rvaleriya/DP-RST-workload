# ==================== Load processed Visium HD data ====================
library(Seurat)
library(Matrix)
library(tidyverse)
library(sf)
library(patchwork)
library(viridis)

set.seed(42)

setwd("~/Desktop/DP-RST-workload/Datasets Results/Visium_HD_Human_Colon_Cancer")

# Paths
rds_file   <- "Nuclei_Data/colon_hd_nuclei_seurat_object_processed.rds"
rdata_file <- "Nuclei_Data/colon_hd_nuclei_data_processed.RData"
boundary_file <- "boundary_outer_filtered.csv"

# Load Seurat object
so <- readRDS(rds_file)

# Load other objects (spatial_keep, counts_keep, etc.)
load(rdata_file)

# Load boundary
bnd <- read.csv(boundary_file, check.names = FALSE)
# --- 1) Extract PCs 1–10 and join with coordinates -----------------------------
pc_mat <- Embeddings(so, "pca")
stopifnot(ncol(pc_mat) >= 10)
pc_mat <- pc_mat[, 1:10, drop = FALSE]         # columns PC_1 ... PC_10

pc_df <- as.data.frame(pc_mat)
pc_df$id <- rownames(pc_mat)

# make sure IDs are same type for joining
spatial_keep <- spatial_keep %>% mutate(id = as.character(id))
pc_df$id     <- as.character(pc_df$id)

df_pc <- spatial_keep %>%
  select(id, x, y) %>%
  inner_join(pc_df, by = "id")

# --- 2) Prepare boundary dataframe (already read as `bnd`) ---------------------
# bnd must have columns X, Y, polygon_id (as in your script)
stopifnot(all(c("X","Y","polygon_id") %in% names(bnd)))
bnd_df <- bnd %>% mutate(polygon_id = as.factor(polygon_id))

# --- 3) Helper to draw one PC --------------------------------------------------
make_pc_plot <- function(pc_colname) {
  ggplot(df_pc, aes(x = x, y = y, color = .data[[pc_colname]])) +
    geom_point(size = 0.25, alpha = 0.9) +
    geom_path(data = bnd_df, aes(x = X, y = Y, group = polygon_id),
              inherit.aes = FALSE, linewidth = 0.35) +
    coord_equal() +
    scale_color_viridis(option = "C", direction = 1, name = pc_colname) +
    labs(title = pc_colname, x = NULL, y = NULL) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.ticks = element_blank(),
      axis.text  = element_blank(),
      legend.position = "right",
      legend.key.height = unit(0.4, "cm"),
      legend.key.width  = unit(0.2, "cm"),
      legend.text = element_text(size = 6)
    )
}
# --- 4) Build PC1–PC10 plots in a 2×5 grid ------------------------------------
pc_names <- paste0("PC_", 1:10)
plots_list <- lapply(pc_names, make_pc_plot)

panel <- wrap_plots(plots_list, ncol = 5) +
  plot_annotation(
    title = "Visium HD — Spatial Maps of PC1–PC10",
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  )

# --- 5) Show in session --------------------------------------------------------
print(panel)

# --- 6) Save to PDF -------------------------------------------------
out_pdf <- file.path("Nuclei_Data", "Colon_VisiumHD_PC1-10_Spatial.pdf")
ggsave(out_pdf, plot = panel, device = cairo_pdf, width = 14, height = 6.5, dpi = 300)
message("Saved: ", out_pdf)