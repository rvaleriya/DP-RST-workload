###############################################################################
#                BayesSpace on STARmap (Radius-Based Neighbors)
#-------------------------------------------------------------------------------
#  - Handles irregular coordinates (e.g., STARmap, MERFISH)
#  - Runs for 3 PCs and 10 PCs
#  - Computes ARI, plots results, and saves labeled tibble as RDS
###############################################################################
set.seed(42)
# ========================== 1. Libraries =====================================
library(SingleCellExperiment)
library(SpatialExperiment)
library(BayesSpace)
library(Matrix)
library(FNN)
library(RANN)
library(ggplot2)
library(mclust)
library(tibble)
library(dplyr)

# ========================== 2. Paths & Load ==================================
setwd("~/Desktop/DP-RST-workload/Datasets Results/Starmap_Mouse")
sce <- readRDS("starmap_sce_with_pca.rds")

sce <- SpatialExperiment(
  assays = assays(sce),
  colData = colData(sce),
  reducedDims = reducedDims(sce),
  spatialCoords = as.matrix(colData(sce)[, c("x", "y")])
)

# ========================== 3. Helpers =======================================
Mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

plot_spatial <- function(df, color_col, title_text) {
  stopifnot(color_col %in% names(df))
  df$.grp <- factor(df[[color_col]])
  ggplot(df, aes(x = x, y = y, color = .grp)) +
    geom_point(size = 0.6, alpha = 0.9) +
    coord_fixed() +
    labs(title = title_text, color = color_col, x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"))
}

# ========================== 4. Neighbor Graph ================================
find_neighbors_radius <- function(sce, q_ref = 0.5, tau = 2, max_k = 50, mutual = TRUE) {
  coords <- as.matrix(spatialCoords(sce)[, c("x", "y")])
  d1 <- FNN::get.knn(coords, k = 1)$nn.dist[, 1]
  r  <- tau * quantile(d1, probs = q_ref, names = FALSE)
  
  nn <- RANN::nn2(coords, query = coords, k = max_k, searchtype = "radius", radius = r)
  idx <- nn$nn.idx
  n <- nrow(coords)
  
  A <- matrix(FALSE, n, n)
  for (i in seq_len(n)) {
    js <- idx[i, ]; js <- js[js > 0 & js != i]
    if (length(js)) A[i, js] <- TRUE
  }
  if (mutual) A <- A & t(A)
  diag(A) <- FALSE
  
  df_j <- lapply(seq_len(n), function(i) which(A[i, ]) - 1L)
  deg <- vapply(df_j, length, integer(1))
  msg <- paste(capture.output(summary(deg)), collapse = " | ")
  message(sprintf("Radius r=%.3f | mean deg=%.2f | mutual=%s | %s",
                  r, mean(deg), mutual, msg))
  df_j
}

# ========================== 5. Run BayesSpace ================================
run_bayes_custom <- function(sce, df_j, n_pcs, q = 7, burn.in = 40000, nrep = 50000) {
  cat(sprintf("\n>> Running BayesSpace with %d PCs...\n", n_pcs))
  Y <- reducedDim(sce, "PCA")[, 1:n_pcs, drop = FALSE]
  init <- BayesSpace:::.init_cluster(Y, q = q, init.method = "mclust")
  
  res <- BayesSpace:::cluster(
    Y, q = q, df_j = df_j, init = init,
    model = "t", precision = "equal",
    mu0 = colMeans(Y), lambda0 = diag(0.01, ncol(Y)),
    gamma = 2, alpha = 1, beta = 0.01, nrep = nrep
  )
  
  # Extract labels (mirroring spatialCluster)
  zs <- res$z[(burn.in + 1):nrep, , drop = FALSE]
  labels <- apply(zs, 2, Mode)
  as.integer(labels)
}

# ========================== 6. Execute for 3 PCs & 10 PCs ====================
df <- as.data.frame(colData(sce))

# Build radius-based neighbors once
df_j <- find_neighbors_radius(sce)

# ---- Run for 3 PCs ----
labels_3 <- run_bayes_custom(sce, df_j, n_pcs = 3)
ari_3 <- adjustedRandIndex(as.integer(as.factor(df$ground_truth)), labels_3)
cat(sprintf("ARI (3 PCs): %.4f\n", ari_3))

df_3 <- tibble(ID = rownames(df), x = df$x, y = df$y,
               label = labels_3, truth = df$ground_truth)

plot_spatial(df_3, "label", sprintf("BayesSpace (3 PCs), ARI = %.3f", ari_3))

saveRDS(df_3, "Results/BayesSpace_STARmap_3PCs.rds")

# ---- Run for 10 PCs ----
labels_10 <- run_bayes_custom(sce, df_j, n_pcs = 10)
ari_10 <- adjustedRandIndex(as.integer(as.factor(df$ground_truth)), labels_10)
cat(sprintf("ARI (10 PCs): %.4f\n", ari_10))

df_10 <- tibble(ID = rownames(df), x = df$x, y = df$y,
                label = labels_10, truth = df$ground_truth)

plot_spatial(df_10, "label", sprintf("BayesSpace (10 PCs), ARI = %.3f", ari_10))

saveRDS(df_10, "Results/BayesSpace_STARmaps_10PCs.rds")

###############################################################################
cat(" Finished both runs.\n")
###############################################################################