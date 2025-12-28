#!/usr/bin/env Rscript

# ================================
# Parallel shard preprocessing (single-node, multicore) for BAST
# ================================
t0 <- Sys.time()

# 0) Library path & RNG -------------------------
my_lib_path <- "/scratch/user/varogovchenko/Rlibs"
.libPaths(my_lib_path)

RNGkind("L'Ecuyer-CMRG")
set.seed(42)

source("/scratch/user/varogovchenko/BASTION_HPRC/Extra_fun/MST_clusters_shredding.R")
source("/scratch/user/varogovchenko/BASTION_HPRC/Extra_fun/MST_from_assignments.R")
source("/scratch/user/varogovchenko/BASTION_HPRC/Extra_fun/merge_all_MST_components.R")

# 1) Libraries ---------------------------------
library(DP.RST)
library(MASS)
library(igraph)
library(fields)
library(mclust)
library(dplyr)
library(purrr)
library(fastcluster)
library(fdaPDE)
library(BASS)
library(readr)
library(future)
library(furrr)

# 2) Paths & parameters -------------------------
base_dir <- "/scratch/user/varogovchenko/BASTION_HPRC/Lung_Xenium"
setwd(base_dir)

S <- 28
M <- 10
temp <- seq(0.1, 1, length.out = M)

out_dir <- file.path(base_dir, "BAST_Shards_3PC")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("[INFO] Working dir: %s\n", base_dir))
cat(sprintf("[INFO] Output dir : %s\n", out_dir))
cat(sprintf("[INFO] Shards (S) : %d (disjoint, no anchors)\n", S))

# 3) Load input data (FIX: use primary_file) ----
primary_file <- file.path(base_dir, "Data_VUILD96MF/Lung_xenium_data_processed.RData")
cat(sprintf("[INFO] Loading data from: %s\n", primary_file))
load(primary_file)

# Sanity checks for required objects:
req <- c("spatial_keep", "pca_embeddings", "counts_keep")  # counts_keep must exist here
missing <- req[!vapply(req, exists, logical(1))]
if (length(missing) > 0) stop("Missing objects in RData: ", paste(missing, collapse=", "))

boundary_file <- file.path(base_dir, "boundary_outer.csv")
cat(sprintf("[INFO] Loading boundary data from: %s\n", boundary_file))
bnd <- read.csv(boundary_file, check.names = FALSE)

# Extract & standardize bases
loc       <- as.matrix(spatial_keep[, c("x_centroid", "y_centroid")])
ann       <- as.matrix(spatial_keep[, c("Annotation_Type")])
Y_sample  <- as.matrix(pca_embeddings[, c("PC_1", "PC_2", "PC_3")])
n         <- nrow(Y_sample)
p         <- ncol(Y_sample)
cat(sprintf("[INFO] Observations: %d | Features (PCs): %d\n", n, p))

# 4) Disjoint shards ----------------------------
sizes <- rep(floor(n / S), S)
remainder <- n %% S
if (remainder > 0) sizes[seq_len(remainder)] <- sizes[seq_len(remainder)] + 1
indices <- sample.int(n)
shards  <- split(indices, rep(seq_len(S), times = sizes))
cat("[INFO] Constructed shard sizes:\n")
print(vapply(shards, length, integer(1)))

# 5) Helpers ------------------------------------

run_bass_once_sce <- function(counts, coords, pcs) {
  # counts: genes x cells; coords: cells x 2; pcs: cells x n_pcs
  barcodes <- colnames(counts)
  rownames(coords) <- barcodes
  
  bass <- createBASSObject(
    X  = list(counts),
    xy = list(coords),
    C  = 20,
    R  = 14,
    beta_method = "SW"
  )
  
  pcs_feat <- t(pcs)             # n_pcs x cells
  colnames(pcs_feat) <- barcodes
  bass@X_run <- pcs_feat
  
  bass <- BASS.run(bass)
  bass <- BASS.postprocess(bass, adjustLS = TRUE)
  
  zlabels <- as.integer(bass@results$z[[1]])
  
  tibble(
    ID    = barcodes,
    x     = coords[, 1],
    y     = coords[, 2],
    label = zlabels
  )
}


preprocess_shard <- function(idx, shard_id) {
  y_shard      <- Y_sample[idx, , drop = FALSE]
  loc_shard    <- loc[idx, , drop = FALSE]                 # matrix [n_s x 2]
  counts_shard <- counts_keep[, idx, drop = FALSE]         # genes x cells (align with idx)
  ann_shard    <- ann[idx]
  
  # Scaled coords for graph building
  coords_sc <- scale(loc_shard)
  Y_std     <- scale(y_shard)
  
  n_s <- nrow(y_shard)
  p_s <- ncol(y_shard)
  
  # Scale boundary using numeric columns (x,y)
  bnd_scaled <- data.frame(
    x = (bnd$x - mean(loc_shard[, 1])) / sd(loc_shard[, 1]),
    y = (bnd$y - mean(loc_shard[, 2])) / sd(loc_shard[, 2])
  )
  
  # --- Clustering with BASS (unchanged) ---------------------------------------
  res_bass <- run_bass_once_sce(
    counts = counts_shard,         # genes x cells
    coords = loc_shard,            # cells x 2
    pcs    = y_shard               # cells x p_s
  )
  res_bass$x_sc <- coords_sc[, 1]
  res_bass$y_sc <- coords_sc[, 2]
  
  # --- Build MST + spatial components; merge tiny comps ALONG MST ONLY --------
  init_partition <- compute_MST_spatial(
    data          = res_bass,
    coords_cols   = c("x_sc","y_sc"),
    cluster_col   = "label",
    bnd           = bnd_scaled,
    threshold     = 5000,
    penalty       = NULL,
    min_comp_size = 10
  )
  
  k_orig   <- length(unique(init_partition$data$spatial_cluster))
  k_merged <- length(unique(init_partition$spatial_cluster_merged))
  cat(sprintf("[Shard %03d] Spatial comps before: %d | after MST-merge: %d\n",
              shard_id, k_orig, k_merged))
  
  clust_s   <- as.integer(init_partition$spatial_cluster_merged)
  graph0    <- init_partition$graph_base
  mstgraph0 <- init_partition$mst
  
  comp_no <- igraph::components(graph0)$no
  cat(sprintf("[Shard %03d] graph0 components: %d\n", shard_id, comp_no))
  if (comp_no > 1L) {
    cat(sprintf("[Shard %03d] Connecting components in graph0...\n", shard_id))
    graph0 <- merge_all_components_to_giant(graph0, coords_sc, verbose = FALSE)
    cat(sprintf("[Shard %03d] Components after connect: %d\n",
                shard_id, igraph::components(graph0)$no))
  }
  
  if ("raw_len" %in% igraph::edge_attr_names(mstgraph0))
    mstgraph0 <- igraph::delete_edge_attr(mstgraph0, "raw_len")
  if ("pen_w" %in% igraph::edge_attr_names(mstgraph0))
    mstgraph0 <- igraph::delete_edge_attr(mstgraph0, "pen_w")
  
  # --- Build BAST initials (unchanged logic) ----------------------------------
  Y_std <- scale(y_shard)
  k_max_s <- length(unique(clust_s))
  
  cluster_means_matrix <- matrix(0, nrow = k_max_s, ncol = p_s)
  for (i in seq_len(k_max_s)) {
    idx_i <- which(clust_s == i)
    cluster_means_matrix[i, ] <- colMeans(Y_std[idx_i, , drop = FALSE])
  }
  
  cluster_mat <- matrix(rep(clust_s, times = M), nrow = n_s, ncol = M, byrow = FALSE)
  mu_list    <- vector("list", M)
  sigmasq_y  <- vector("list", M)
  mst_list   <- vector("list", M)
  covY       <- cov(Y_std)
  
  for (m in seq_len(M)) {
    mu_list[[m]]    <- cluster_means_matrix
    sigmasq_y[[m]]  <- covY
    mst_list[[m]]   <- mstgraph0
  }
  
  init_val <- list(
    trees     = mst_list,
    mu        = mu_list,
    cluster   = cluster_mat,
    sigmasq_y = sigmasq_y
  )
  
  hyperpar <- list(
    sigmasq_mu = (0.5/(2*sqrt(1)))^2,
    lambda_s   = diag(1, p_s),
    nu         = p_s,
    M          = M,
    k_max      = k_max_s,
    lambda_k   = 15
  )
  
  idx_global <- idx
  save(
    list = c("y_shard","loc_shard","graph0","init_val","hyperpar","clust_s","idx_global","ann_shard"),
    file = file.path(out_dir, sprintf("Shard_%03d_Preproc_BAST.RData", shard_id))
  )
  cat(sprintf("[OK] Shard %03d → %s\n", shard_id,
              file.path(out_dir, sprintf("Shard_%03d_Preproc_BAST.RData", shard_id))))
}

# 6) Parallel backend ----------------------------
slurm_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
num_workers <- if (!is.na(slurm_cores)) slurm_cores else future::availableCores()
plan(multicore, workers = num_workers)
options(mc.cores = num_workers)
Sys.setenv(OMP_NUM_THREADS = "1")  # avoid nested OpenMP
cat("[INFO] Using", num_workers, "workers\n")

# 7) Launch all shards in parallel -------------
# --- Use forked workers & avoid exporting big globals ---
stopifnot(.Platform$OS.type == "unix")                # ensure Linux/macOS
options(future.globals.maxSize = Inf)                 # or "8GiB"; protects rare fallbacks
plan(multicore, workers = num_workers)                # keep fork, not multisession

opts <- furrr_options(seed = TRUE,
                      scheduling = Inf,
                      stdout = TRUE,
                      globals = FALSE)                # <-- key: don't export globals

future_walk(seq_along(shards),
            function(sid) preprocess_shard(shards[[sid]], sid),
            .options = opts)

cat("✅  All shard preprocessing files written to: ", out_dir, "\n")
cat("Total runtime (min): ", round(difftime(Sys.time(), t0, units = "mins"), 2), "\n")