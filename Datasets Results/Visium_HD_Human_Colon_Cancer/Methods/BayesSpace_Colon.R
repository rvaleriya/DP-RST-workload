###############################################################################
#                    BayesSpace on Colon Visium HD (Radius Neighbors)
#   - Irregular coordinates (nuclei) handled via radius-based mutual graph
#   - Runs for 3 PCs and 10 PCs
#   - Computes ARI (if Annotation_Type present), saves RDS + PNGs
###############################################################################

writeLines("R script started", stderr())
start.time_file <- Sys.time()
set.seed(42)

# ========================== 1. Libraries & LibPath ============================
my_lib_path <- "/scratch/user/varogovchenko/Rlibs"
.libPaths(my_lib_path)

library(BayesSpace)
library(Matrix)
library(dplyr)
library(ggplot2)
library(mclust)
library(FNN)
library(RANN)
library(tibble)
library(FNN)
library(data.table)

# ========================== 2. Paths & Data ==================================
work_dir <- "/scratch/user/varogovchenko/BayesSpace_runs/Colon"
data_rdata <- "/scratch/user/varogovchenko/BASTION_HPRC/Visium_HD_Human_Colon_Cancer/Nuclei_Data/colon_hd_nuclei_data_processed.RData"

dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)

load(data_rdata)

# ========================== 3. Sanity & Align ================================

xy_mat  <- as.matrix(spatial_keep[, c("x","y")])
truth_vec <- spatial_keep$annotation

cat("Counts dim (genes x cells):", paste(dim(counts_keep), collapse=" x "), "\n")
cat("PC matrix dim (cells x PCs):", paste(dim(pca_embeddings), collapse=" x "), "\n")
cat("XY matrix dim (cells x 2):", paste(dim(xy_mat), collapse=" x "), "\n")

# ========================== 4. Helpers =======================================
Mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

plot_spatial <- function(df, label_col, title_text) {
  stopifnot(label_col %in% names(df))
  df$.lab <- factor(df[[label_col]])
  ggplot(df, aes(x = x, y = y, color = .lab)) +
    geom_point(size = 0.4, alpha = 0.9) +
    coord_equal() +
    labs(title = title_text, color = label_col, x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"))
}

# Build mutual radius neighbors -> list-of-integers df_j (0-based indices)

find_neighbors_knn <- function(xy, k = 6, mutual = TRUE) {
  stopifnot(is.matrix(xy) && ncol(xy) == 2)
  n <- nrow(xy)

  # kNN (directed i -> j)
  nn  <- FNN::get.knn(xy, k = k)
  idx <- nn$nn.index            # n × k (1-based indices)

  # Build edge list
  ii <- rep.int(seq_len(n), times = ncol(idx))
  jj <- as.vector(idx)
  E  <- data.table(i = ii, j = jj)

  if (mutual) {
    # Keep only mutual edges (undirected core), then duplicate to both directions
    E[, a := pmin(i, j)][, b := pmax(i, j)]
    Em <- E[, .N, by = .(a, b)][N >= 2L][, .(i = a, j = b)]
    E  <- rbindlist(list(Em, Em[, .(i = j, j = i)]), use.names = FALSE)
  }

  # --- Key change: split over a factor with levels 1:n (guarantees length n) ---
  E[, i := as.integer(i)]
  df_j <- split(E$j, factor(E$i, levels = seq_len(n)))  # list length n
  # Convert to 0-based integer neighbor lists (empty for isolates)
  out <- lapply(df_j, function(v) if (length(v)) as.integer(v - 1L) else integer(0L))

  # Degree summary
  deg <- vapply(out, length, integer(1))
  msg <- paste(capture.output(summary(deg)), collapse = " | ")
  message(sprintf("Neighbor graph (k=%d): mean deg=%.2f | mutual=%s | %s",
                  k, mean(deg), mutual, msg))

  out
}

# Internal BayesSpace call 
# -------------------- Add this once, right after loading data -----------------
# Align and set barcodes
stopifnot(nrow(pca_embeddings) == nrow(spatial_keep))
barcodes <- rownames(pca_embeddings)
if (is.null(barcodes)) barcodes <- spatial_keep$cell_id
stopifnot(length(barcodes) == nrow(pca_embeddings))
# -----------------------------------------------------------------------------

# ========================= Robust BayesSpace wrapper ==========================
run_bayes_custom <- function(Y, df_j, n_pcs, q = 6, burn.in = 40000, nrep = 50000) {
  if (!requireNamespace("BayesSpace", quietly = TRUE)) {
    stop("BayesSpace not installed/available on .libPaths()")
  }
  cat(sprintf("\n>> Running BayesSpace with %d PCs (q=%d, burn.in=%d, nrep=%d)...\n",
              n_pcs, q, burn.in, nrep))

  Yp <- as.matrix(Y[, seq_len(n_pcs), drop = FALSE])

  # Initialize with mclust
  init <- BayesSpace:::.init_cluster(Yp, q = q, init.method = "mclust")

  # Fit (MCMC)
  res <- BayesSpace:::cluster(
    Yp,
    q = q,
    df_j = df_j,
    init = init,
    model = "t",
    precision = "equal",
    mu0 = colMeans(Yp),
    lambda0 = diag(0.01, ncol(Yp)),
    gamma = 2,
    alpha = 1,
    beta = 0.01,
    nrep = nrep
  )

  # ---- Robust extraction of labels after burn-in ----
  z <- res$z
  n_cells <- nrow(Yp)

  # Helper: posterior mode (ties resolved by first)
  Mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

  if (is.null(z)) {
    stop("BayesSpace:::cluster returned NULL z; fitting likely failed upstream.")
  }

  # Case A: z is a vector (last state only)
  if (is.atomic(z) && is.null(dim(z))) {
    if (length(z) != n_cells) {
      stop(sprintf("Unexpected z length (%d) != n_cells (%d).", length(z), n_cells))
    }
    return(as.integer(z))
  }

  # Case B: z is a matrix; detect orientation
  if (!is.matrix(z)) {
    stop("Unexpected z type; expected matrix or vector.")
  }

  # Identify which dimension equals n_cells
  dimz <- dim(z)
  # Orientation 1: rows = iterations, cols = cells
  if (dimz[2] == n_cells) {
    n_iters <- dimz[1]
    first_keep <- min(max(burn.in + 1L, 1L), n_iters)
    if (first_keep > n_iters) {
      stop(sprintf("Burn-in (%d) exceeds number of stored iterations (%d).", burn.in, n_iters))
    }
    z_post <- z[first_keep:n_iters, , drop = FALSE]   # (iters_kept × cells)
    labels <- apply(z_post, 2, Mode)
    return(as.integer(labels))
  }

  # Orientation 2: rows = cells, cols = iterations
  if (dimz[1] == n_cells) {
    n_iters <- dimz[2]
    first_keep <- min(max(burn.in + 1L, 1L), n_iters)
    if (first_keep > n_iters) {
      stop(sprintf("Burn-in (%d) exceeds number of stored iterations (%d).", burn.in, n_iters))
    }
    z_post <- z[, first_keep:n_iters, drop = FALSE]   # (cells × iters_kept)
    labels <- apply(z_post, 1, Mode)
    return(as.integer(labels))
  }

  stop(sprintf(
    "Cannot match z dimensions %s to n_cells=%d.",
    paste(dimz, collapse = "×"), n_cells
  ))
}

# ========================== 5. Build neighbor graph ===========================
df_j <- find_neighbors_knn(xy_mat, k = 6, mutual = TRUE)

# ========================== 6. Run: 10 PCs ===================================
labels_10 <- run_bayes_custom(pca_embeddings, df_j, n_pcs = 10, q = 6, burn.in = 40000, nrep = 50000)

ari_10 <- NA_real_
if (!is.null(truth_vec)) {
  keep <- !is.na(truth_vec)
  if (sum(keep) >= 2L &&
      length(unique(labels_10[keep])) >= 2L &&
      length(unique(truth_vec[keep])) >= 2L) {
    truth_int <- as.integer(as.factor(truth_vec[keep]))
    ari_10 <- mclust::adjustedRandIndex(as.integer(labels_10[keep]), truth_int)
  }
}
cat(sprintf("ARI (10 PCs): %s\n", ifelse(is.na(ari_10), "NA", sprintf("%.4f", ari_10))))

df_10 <- tibble(
  barcode = barcodes,
  x = xy_mat[, 1],
  y = xy_mat[, 2],
  label = labels_10,
  truth = if (is.null(truth_vec)) NA else truth_vec
)
saveRDS(df_10, file = "BayesSpace_Colon_10PCs.rds")

p10 <- plot_spatial(df_10, "label",
                    sprintf("BayesSpace (Colon Xenium, 10 PCs)%s",
                            ifelse(is.na(ari_10), "", sprintf(", ARI = %.3f", ari_10))))
ggsave("BayesSpace_Colon_10PCs_clusters.png", p10, width = 6, height = 5, dpi = 300)

# ========================== 8. Summary & Timing ===============================
cat("\n========== SUMMARY (ARI excludes NA) ==========\n")
cat(sprintf("10 PCs : %s\n", ifelse(is.na(ari_10), "ARI = NA", sprintf("ARI = %.3f", ari_10))))

end.time_file <- Sys.time()
file_time <- end.time_file - start.time_file
print(file_time)
cat("Finished BayesSpace runs on Colon Xenium.\n")
###############################################################################