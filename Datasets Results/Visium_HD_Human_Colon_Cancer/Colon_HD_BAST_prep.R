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
base_dir <- "/scratch/user/varogovchenko/BASTION_HPRC/Visium_HD_Human_Colon_Cancer"
setwd(base_dir)

S <- 100
M <- 10
temp <- seq(0.1, 1, length.out = M)

out_dir <- file.path(base_dir, "BAST_Shards_3PC")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("[INFO] Working dir: %s\n", base_dir))
cat(sprintf("[INFO] Output dir : %s\n", out_dir))
cat(sprintf("[INFO] Shards (S) : %d (disjoint, no anchors)\n", S))

# 3) Load input data (FIX: use primary_file) ----
primary_file <- file.path(base_dir, "Nuclei_Data/colon_hd_nuclei_data_processed.RData")
cat(sprintf("[INFO] Loading data from: %s\n", primary_file))
load(primary_file)

# Sanity checks for required objects:
req <- c("spatial_keep", "pca_embeddings", "counts_keep")  # counts_keep must exist here
missing <- req[!vapply(req, exists, logical(1))]
if (length(missing) > 0) stop("Missing objects in RData: ", paste(missing, collapse=", "))

boundary_file <- file.path(base_dir, "boundary_outer_filtered.csv")
cat(sprintf("[INFO] Loading boundary data from: %s\n", boundary_file))
bnd <- read.csv(boundary_file, check.names = FALSE)

# Extract & standardize bases
loc       <- as.matrix(spatial_keep[, c("x", "y")])
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
    R  = 12,
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

compute_MST_spatial <- function(data,
                                coords_cols = c("x", "y"),
                                cluster_col = "label",
                                bnd = NULL,
                                threshold = 5000,
                                penalty = NULL,
                                min_comp_size = 3) {
  stopifnot(all(coords_cols %in% names(data)), cluster_col %in% names(data))
  coords <- as.matrix(data[, coords_cols])
  mesh   <- gen2dMesh(coords, bnd)
  graph0 <- constrainedDentri(n = nrow(coords), mesh = mesh, threshold = threshold)
  E(graph0)$eid <- seq_len(igraph::ecount(graph0))
  V(graph0)$vid <- seq_len(igraph::vcount(graph0))
  V(graph0)$cluster <- data[[cluster_col]]
  
  orig_w <- E(graph0)$weight
  if (is.null(penalty)) penalty <- max(orig_w, na.rm = TRUE) + 1
  
  ee <- as.data.frame(igraph::as_edgelist(graph0))
  colnames(ee) <- c("V1","V2")
  ee$V1 <- as.integer(ee$V1); ee$V2 <- as.integer(ee$V2)
  c1 <- V(graph0)$cluster[ee$V1]; c2 <- V(graph0)$cluster[ee$V2]
  inter <- (c1 != c2)
  
  E(graph0)$new_weight <- orig_w + penalty * inter
  mst_g <- igraph::mst(graph0, weights = E(graph0)$new_weight)
  
  head_idx <- as.integer(igraph::head_of(mst_g, E(mst_g)))
  tail_idx <- as.integer(igraph::tail_of(mst_g, E(mst_g)))
  drop_cross <- V(mst_g)$cluster[head_idx] != V(mst_g)$cluster[tail_idx]
  mst_same   <- igraph::delete_edges(mst_g, which(drop_cross))
  
  comp <- igraph::components(mst_same)
  spatial_cluster_orig <- comp$membership
  data$spatial_cluster <- spatial_cluster_orig
  
  if (is.null(min_comp_size) || min_comp_size <= 1) {
    return(list(
      mst = mst_g, data = data, mst_same = mst_same,
      components = comp$no, merged = FALSE,
      spatial_cluster_merged = spatial_cluster_orig
    ))
  }
  
  comp_idx_list <- split(seq_len(nrow(data)), spatial_cluster_orig)
  comp_size  <- vapply(comp_idx_list, length, integer(1))
  comp_label <- vapply(comp_idx_list, function(ix) {
    ul <- unique(data[[cluster_col]][ix]); if (length(ul) != 1L) stop("Mixed labels")
    ul
  }, FUN.VALUE = data[[cluster_col]][1])
  comp_cx <- vapply(comp_idx_list, function(ix) mean(coords[ix, 1]), numeric(1))
  comp_cy <- vapply(comp_idx_list, function(ix) mean(coords[ix, 2]), numeric(1))
  comp_df <- data.frame(
    comp_id = as.integer(names(comp_idx_list)),
    size = as.integer(comp_size), label = comp_label, cx = comp_cx, cy = comp_cy,
    stringsAsFactors = FALSE
  )
  
  tiny_mask <- comp_df$size < min_comp_size
  if (any(tiny_mask)) {
    new_membership <- spatial_cluster_orig
    choose_target <- function(row_i) {
      lbl <- comp_df$label[row_i]; cid <- comp_df$comp_id[row_i]
      cx  <- comp_df$cx[row_i];    cy  <- comp_df$cy[row_i]
      same_lbl <- which(comp_df$label == lbl & comp_df$comp_id != cid)
      if (!length(same_lbl)) return(cid)
      big_cand <- same_lbl[comp_df$size[same_lbl] >= min_comp_size]
      cand <- if (length(big_cand)) big_cand else same_lbl
      dx <- comp_df$cx[cand] - cx; dy <- comp_df$cy[cand] - cy
      cand[ which.min(dx*dx + dy*dy) ]
    }
    for (r in which(tiny_mask)) {
      src <- comp_df$comp_id[r]; tgt <- choose_target(r)
      if (tgt != src) new_membership[new_membership == src] <- tgt
    }
    uniq_ids <- sort(unique(new_membership))
    remap <- setNames(seq_along(uniq_ids), uniq_ids)
    spatial_cluster_merged <- unname(remap[as.character(new_membership)])
    data$spatial_cluster_merged <- spatial_cluster_merged
    return(list(
      mst = mst_g, data = data, mst_same = mst_same,
      components = comp$no, merged = TRUE,
      spatial_cluster_merged = spatial_cluster_merged,
      k_merged = length(uniq_ids), comp_summary_before = comp_df
    ))
  } else {
    data$spatial_cluster_merged <- spatial_cluster_orig
    return(list(
      mst = mst_g, data = data, mst_same = mst_same,
      components = comp$no, merged = FALSE,
      spatial_cluster_merged = spatial_cluster_orig,
      k_merged = comp$no, comp_summary_before = comp_df
    ))
  }
}

merge_all_components_to_giant <- function(g, coords, verbose = TRUE) {
  stopifnot(inherits(g, "igraph"))
  if (is.data.frame(coords)) coords <- as.matrix(coords)
  stopifnot(nrow(coords) == vcount(g), ncol(coords) >= 2)
  sq_dists <- function(A, B) {
    AA <- rowSums(A*A); BB <- rowSums(B*B)
    outer(AA, BB, "+") - 2*(A %*% t(B))
  }
  iter <- 0L
  repeat {
    comp <- components(g); k <- comp$no
    if (verbose) cat(sprintf("Iteration %d: components = %d\n", iter, k))
    if (k <= 1L) break
    memb <- comp$membership; sizes <- tabulate(memb, nbins = k)
    giant_id <- which.max(sizes); giant_idx <- which(memb == giant_id)
    A <- coords[giant_idx, 1:2, drop = FALSE]
    other_ids <- setdiff(seq_len(k), giant_id)
    best_pair <- c(NA_integer_, NA_integer_); best_d2 <- Inf
    for (cid in other_ids) {
      comp_idx <- which(memb == cid); B <- coords[comp_idx, 1:2, drop = FALSE]
      D2 <- sq_dists(A, B); wmin <- which.min(D2)
      ii <- ((wmin - 1L) %% nrow(D2)) + 1L
      jj <- ((wmin - 1L) %/% nrow(D2)) + 1L
      d2 <- D2[ii, jj]
      if (d2 < best_d2) { best_d2 <- d2; best_pair <- c(giant_idx[ii], comp_idx[jj]) }
    }
    u <- best_pair[1]; v <- best_pair[2]
    if (!are_adjacent(g, u, v)) {
      g <- add_edges(g, c(u, v)); E(g)$merged_edge <- FALSE; E(g)[ecount(g)]$merged_edge <- TRUE
      if (verbose) cat(sprintf("  + connected via (%d <-> %d), d=%.3f\n", u, v, sqrt(best_d2)))
    }
    iter <- iter + 1L
  }
  E(g)$eid <- seq_len(ecount(g)); g
}

preprocess_shard <- function(idx, shard_id) {
  y_shard      <- Y_sample[idx, , drop = FALSE]
  loc_shard    <- loc[idx, , drop = FALSE]                 # matrix [n_s x 2]
  counts_shard <- counts_keep[, idx, drop = FALSE]         # genes x cells (align with idx)
  
  # Scaled coords for graph building
  coords_sc <- scale(loc_shard)
  Y_std     <- scale(y_shard)
  
  n_s <- nrow(y_shard)
  p_s <- ncol(y_shard)
  
  # Scale boundary using *numeric columns* 1 and 2 (x,y)
  bnd_scaled <- data.frame(
    x = (bnd$x - mean(loc_shard[, 1])) / sd(loc_shard[, 1]),
    y = (bnd$y - mean(loc_shard[, 2])) / sd(loc_shard[, 2])
  )
  
  res_bass <- run_bass_once_sce(
    counts = counts_shard,         # genes x cells
    coords = loc_shard,            # cells x 2
    pcs    = y_shard               # cells x p_s
  )
  
  # Inject scaled coords for MST
  res_bass$x_sc <- coords_sc[, 1]
  res_bass$y_sc <- coords_sc[, 2]
  
  init_partition <- compute_MST_spatial(
    data          = res_bass,
    coords_cols   = c("x_sc","y_sc"),
    cluster_col   = "label",
    bnd           = bnd_scaled,
    threshold     = 5000,
    penalty       = NULL,
    min_comp_size = 5
  )
  
  k_max_s <- length(unique(init_partition$spatial_cluster_merged))
  cat(sprintf("[Shard %03d] Spatial comps after drop: %d\n", shard_id, k_max_s))
  
  mesh_shard <- gen2dMesh(coords_sc, bnd_scaled)
  g_shard    <- constrainedDentri(n_s, mesh_shard)
  E(g_shard)$eid <- seq_len(ecount(g_shard))
  V(g_shard)$vid <- seq_len(vcount(g_shard))
  
  g_merged <- merge_all_components_to_giant(g_shard, coords_sc, verbose = FALSE)
  el <- as_edgelist(g_merged, names = FALSE)
  distance <- sqrt(rowSums((coords_sc[el[, 1], ] - coords_sc[el[, 2], ])^2))
  E(g_merged)$weight <- distance
  
  g_mst     <- mst(g_merged)
  g_mst     <- delete_edge_attr(g_mst, "weight")
  g_merged  <- delete_edge_attr(g_merged, "weight")
  
  graph0    <- g_merged
  mstgraph0 <- g_mst
  
  clust_s <- as.integer(init_partition$spatial_cluster_merged)
  
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
    list = c("y_shard","loc_shard","graph0","init_val","hyperpar","clust_s","idx_global"),
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
opts <- furrr_options(seed = TRUE, scheduling = Inf, stdout = TRUE)
future_walk(seq_along(shards), function(sid) {
  preprocess_shard(shards[[sid]], sid)
}, .options = opts)

cat("✅  All shard preprocessing files written to: ", out_dir, "\n")
cat("Total runtime (min): ", round(difftime(Sys.time(), t0, units = "mins"), 2), "\n")