###############################################################################
#                     DR-SC on Lung Xenium: 10 PCs (K = 6)
###############################################################################

writeLines("R script started", stderr())
start.time_file <- Sys.time()
set.seed(42)

# ========================== Libraries & Paths ================================
my_lib_path <- "/scratch/user/varogovchenko/Rlibs"
.libPaths(my_lib_path)

library(Seurat)
library(DR.SC)
library(dplyr)
library(ggplot2)
library(mclust)
library(tibble)

work_dir   <- "/scratch/user/varogovchenko/DRSC_runs/Lung"
data_rds   <- "/scratch/user/varogovchenko/BASTION_HPRC/Lung_Xenium/Data_VUILD96MF/Lung_xenium_seurat_object_processed.rds"

dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)
out_dir <- getwd()

# ========================== Load Data ========================================
seu <- readRDS(data_rds)
meta <- seu@meta.data

meta$row <- meta$x_centroid
meta$col <- meta$y_centroid
if (!"Annotation_Type" %in% names(meta)) meta$Annotation_Type <- NA
rownames(meta) <- colnames(seu)
seu@meta.data <- meta

# ========================== Run DR-SC (10 PCs) ================================
K <- 6
q_dim <- 10

cat(sprintf(">> Running DR-SC (K=%d, q=%d) ...\n", K, q_dim))
fit <- DR.SC(
  seu,
  K        = K,
  q        = q_dim,
  platform = "Other_SRT",
  verbose  = FALSE
)

meta$label <- as.integer(fit$spatial.drsc.cluster)
ari_val <- ifelse(all(is.na(meta$Annotation_Type)), NA,
                  adjustedRandIndex(as.integer(as.factor(meta$Annotation_Type)), meta$label))

# ========================== Save & Plot ======================================
out_tbl <- tibble(
  ID    = colnames(seu),
  row   = meta$row,
  col   = meta$col,
  label = meta$label,
  truth = meta$Annotation_Type
)

rds_path <- file.path(out_dir, sprintf("DRSC_Lung_%dPCs.rds", q_dim))
saveRDS(out_tbl, rds_path)

p <- ggplot(out_tbl, aes(x = col, y = row, color = factor(label))) +
  geom_point(size = 0.4, alpha = 0.9) +
  coord_fixed() +
  labs(title = sprintf("DR-SC (K=%d, q=%d)%s",
                       K, q_dim,
                       ifelse(is.na(ari_val), "", sprintf(", ARI = %.3f", ari_val))),
       color = "Cluster") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, sprintf("DRSC_Lung_%dPCs.png", q_dim)),
       p, width = 6, height = 5, dpi = 300)

cat("\n========== SUMMARY ==========\n")
cat(sprintf("10 PCs : %s\n",
            ifelse(is.na(ari_val), "ARI = NA", sprintf("ARI = %.6f", ari_val))))
cat("Saved:\n")
cat("  - ", rds_path, "\n", sep = "")

end.time_file <- Sys.time()
print(end.time_file - start.time_file)
cat("Finished DR-SC run on Lung Xenium.\n")
###############################################################################