###############################################################################
#                     DR-SC on Colon Xenium: 10 PCs (K = 6)
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
library(FNN)

work_dir   <- "/scratch/user/varogovchenko/DRSC_runs/Colon"
data_rds   <- "/scratch/user/varogovchenko/BASTION_HPRC/Visium_HD_Human_Colon_Cancer/Nuclei_Data/colon_hd_nuclei_seurat_object_processed.rds"

dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)
out_dir <- getwd()

# ========================== Load Data ========================================
seu <- readRDS(data_rds)
meta <- seu@meta.data

meta$row <- meta$x
meta$col <- meta$y
if (!"annotation" %in% names(meta)) meta$annotation <- NA
rownames(meta) <- colnames(seu)
seu@meta.data <- meta

# ========================== Run DR-SC (10 PCs) ================================
K <- 6
q_dim <- 10

cat(sprintf(">> Running DR-SC (K=%d, q=%d) ...\n", K, q_dim))
# build adjacency once and pass it to DR-SC
coords <- as.matrix(seu@meta.data[, c("row","col")])
idx <- sample(seq_len(nrow(coords)), size = min(2000, nrow(coords)))
Csub <- coords[idx, , drop = FALSE]

# Use kNN distances to estimate a radius that hits ~8 neighbors
knn8 <- FNN::get.knn(Csub, k = 8)
# distance to the 8th neighbor per sampled spot
r8 <- median(knn8$nn.dist[, 8], na.rm = TRUE)
print(r8)

# Give getAdj_auto a search bracket around that radius
lower_target_neighbors <- 4   # median neighbors should exceed this
Adj <- DR.SC::getAdj_auto(
  pos       = coords,
  lower.med = lower_target_neighbors,
  # Let upper.med be comfortably above r8 (e.g., 3x) so the search won’t top out
  upper.med = ceiling(3 * r8)
)

# Build X (cells x genes) as a *sparse* matrix
X <- Matrix::t(SeuratObject::GetAssayData(
  seu, assay = SeuratObject::DefaultAssay(seu), layer = "data"
))   # Matrix::t keeps it sparse

# Sanity check (should print 'dgCMatrix')
message(class(X)[1])

fit <- DR.SC_fit(
  X,
  K        = 6,
  q        = 10,
  Adj_sp   = Adj,          # << forces DR-SC to use our graph
  verbose  = FALSE
)

lab <- as.integer(fit$Objdrsc[[1]]$cluster)
print(table(lab))

# Labels back to meta in the SAME order
meta <- seu@meta.data
meta$label <- lab

ari_val <- ifelse(all(is.na(meta$annotation)), NA,
                  adjustedRandIndex(as.integer(as.factor(meta$annotation)), meta$label))

# ========================== Save & Plot ======================================
out_tbl <- tibble(
  ID    = colnames(seu),
  row   = meta$row,
  col   = meta$col,
  label = meta$label,
  truth = meta$annotation
)

rds_path <- file.path(out_dir, sprintf("DRSC_Colon_%dPCs.rds", q_dim))
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

ggsave(file.path(out_dir, sprintf("DRSC_Colon_%dPCs.png", q_dim)),
       p, width = 6, height = 5, dpi = 300)

cat("\n========== SUMMARY ==========\n")
cat(sprintf("10 PCs : %s\n",
            ifelse(is.na(ari_val), "ARI = NA", sprintf("ARI = %.6f", ari_val))))
cat("Saved:\n")
cat("  - ", rds_path, "\n", sep = "")

end.time_file <- Sys.time()
print(end.time_file - start.time_file)
cat("Finished DR-SC run on Colon Xenium.\n")
###############################################################################