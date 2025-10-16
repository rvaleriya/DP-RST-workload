# ---------------------------------------------------------------------------
# Consensus-MCMC preprocessing for BAST (S = 100 disjoint shards, NO anchors)
# Directory: /scratch/user/varogovchenko/BASTION_HPRC/Visium_HD_Human_Colon_Cancer
# ---------------------------------------------------------------------------

# custom lib
my_lib_path <- "/scratch/user/varogovchenko/Rlibs"        
.libPaths(my_lib_path)

## 0) Libraries --------------------------------------------------------------
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

## 1) Paths & parameters -----------------------------------------------------
set.seed(42)

base_dir <- "/scratch/user/varogovchenko/BASTION_HPRC/Visium_HD_Human_Colon_Cancer"
setwd(base_dir)

# Total number of shards (disjoint)
S <- 100

# Parallel chains (for later PT runs)
M <- 10
temp <- seq(0.1, 1, length.out = M)  # retained if you want to log it

# Output directory (concrete)
out_dir <- file.path(base_dir, "BAST_Shards_3PC")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("[INFO] Working dir: %s\n", base_dir))
cat(sprintf("[INFO] Output dir : %s\n", out_dir))
cat(sprintf("[INFO] Shards (S) : %d (disjoint, no anchors)\n", S))

## 2) Load input data  ----------------
primary_file <- file.path(base_dir, "Nuclei_Data/colon_hd_nuclei_data_processed.RData")

cat(sprintf("[INFO] Loading data from: %s\n", primary_file))
load(candidate_file)  

boundary_file <- file.path(base_dir, "boundary_outer_filtered.csv")
cat(sprintf("[INFO] Loading boundary data from: %s\n", boundary_file))
bnd <- read.csv(boundary_file, check.names = FALSE) 

# Extract & standardize
loc      <- as.matrix(spatial_keep[, c("x", "y")])
Y_sample <- as.matrix(pca_embeddings[, c("PC_1", "PC_2", "PC_3")])

n <- nrow(Y_sample)
p <- ncol(Y_sample)

cat(sprintf("[INFO] Observations: %d | Features (PCs): %d\n", n, p))

## 3) Disjoint shards (balanced) ---------------------------------------------
# Make shard sizes as balanced as possible, no anchors.
sizes <- rep(floor(n / S), S)
remainder <- n %% S
if (remainder > 0) sizes[seq_len(remainder)] <- sizes[seq_len(remainder)] + 1

indices <- sample.int(n)  # shuffle for random, balanced splits
shards <- split(indices, rep(seq_len(S), times = sizes))

cat("[INFO] Constructed disjoint shards with sizes:\n")
print(vapply(shards, length, integer(1)))

## 4) Function to run BASS ---------------------------------------------

run_bass_once_sce <- function(counts, coords, pcs) {

  # Counts (genes x cells) and coordinates (cells x 2)
  barcodes  <- colnames(counts)
  rownames(coords) <- barcodes
  
  message(">> Running BASS.")
  
  # Create BASS object with counts (features x cells)
  bass <- createBASSObject(
    X  = list(counts),     # genes x cells
    xy = list(coords),        # cells x 2
    C  = 20, 
    R  = 12,
    beta_method = "SW"
  )
  
  # Inject top n_pcs PCs (cells x n_pcs) -> transpose to features x cells
  pcs_feat  <- t(pcs)                                            # n_pcs x cells
  colnames(pcs_feat) <- barcodes                                       # align explicitly
  bass@X_run <- pcs_feat
  
  # Run + postprocess
  bass <- BASS.run(bass)
  bass <- BASS.postprocess(bass, adjustLS = TRUE)
  
  # Labels and evaluation
  zlabels <- as.integer(bass@results$z[[1]])
  
  # Output table
  out_tbl <- tibble(
    ID    = barcodes,
    x     = coords[, "x"],
    y     = coords[, "y"],
    label = zlabels
  )
  
  list(labels_tbl = out_tbl, bass = bass)
}

# ------------------------- Helper: penalized MST ----------------------------
compute_MST_spatial <- function(data,
                                coords_cols = c("x", "y"),
                                cluster_col = "label",
                                bnd = NULL,
                                threshold = 5000,
                                penalty = NULL,
                                min_comp_size = 3) {
  stopifnot(all(coords_cols %in% names(data)),
            cluster_col %in% names(data))
  
  # ---- 1) Build penalized MST exactly as you had ----------------------------
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
  
  c1 <- V(graph0)$cluster[ee$V1]
  c2 <- V(graph0)$cluster[ee$V2]
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
  
  # ---- Early exit if no merging requested -----------------------------------
  if (is.null(min_comp_size) || min_comp_size <= 1) {
    return(list(
      mst = mst_g,
      data = data,
      mst_same = mst_same,
      components = comp$no,
      merged = FALSE,
      spatial_cluster_merged = spatial_cluster_orig
    ))
  }
  
  # ---- 2) Component summary (size, label, centroid) -------------------------
  comp_ids <- seq_len(comp$no)
  # Map each comp_id to indices
  # Build a list of indices per component
  comp_idx_list <- split(seq_len(nrow(data)), spatial_cluster_orig)
  
  # For each component, compute size, label, centroid
  comp_size <- vapply(comp_idx_list, length, integer(1))
  comp_label <- vapply(comp_idx_list, function(ix) {
    ul <- unique(data[[cluster_col]][ix])
    if (length(ul) != 1L) stop("Component with mixed labels detected; unexpected after cross-edge drop.")
    ul
  }, FUN.VALUE = data[[cluster_col]][1])
  comp_cx <- vapply(comp_idx_list, function(ix) mean(coords[ix, 1]), numeric(1))
  comp_cy <- vapply(comp_idx_list, function(ix) mean(coords[ix, 2]), numeric(1))
  
  comp_df <- data.frame(
    comp_id = as.integer(names(comp_idx_list)),
    size    = as.integer(comp_size),
    label   = comp_label,
    cx      = comp_cx,
    cy      = comp_cy,
    stringsAsFactors = FALSE
  )
  
  # ---- 3) Determine merge targets for tiny components -----------------------
  tiny_mask <- comp_df$size < min_comp_size
  if (any(tiny_mask)) {
    # Initialize new membership vector (mutable)
    new_membership <- spatial_cluster_orig
    
    # Precompute per-label lookups to speed distance queries
    # For each label, we’ll hold a matrix of candidate components (>= min size) and their centroids.
    labels_all <- unique(comp_df$label)
    
    # Helper to find target comp_id for a given tiny component row
    choose_target <- function(row_i) {
      lbl  <- comp_df$label[row_i]
      cid  <- comp_df$comp_id[row_i]
      cx   <- comp_df$cx[row_i]
      cy   <- comp_df$cy[row_i]
      
      same_lbl <- which(comp_df$label == lbl & comp_df$comp_id != cid)
      if (length(same_lbl) == 0L) {
        # No other component with same label; nothing to merge into
        return(cid)
      }
      # Prefer "big" components within same label
      big_cand <- same_lbl[comp_df$size[same_lbl] >= min_comp_size]
      cand <- big_cand
      if (length(cand) == 0L) {
        # If no big candidates, fall back to any same-label component
        cand <- same_lbl
      }
      # Nearest by centroid distance
      dx <- comp_df$cx[cand] - cx
      dy <- comp_df$cy[cand] - cy
      j  <- cand[ which.min(dx*dx + dy*dy) ]
      return(comp_df$comp_id[j])
    }
    
    # For each tiny component, reassign all its nodes to the chosen target
    tiny_rows <- which(tiny_mask)
    for (r in tiny_rows) {
      src_cid <- comp_df$comp_id[r]
      tgt_cid <- choose_target(r)
      if (tgt_cid != src_cid) {
        new_membership[new_membership == src_cid] <- tgt_cid
      }
    }
    
    # ---- 4) Renumber merged component ids to 1..K_merged --------------------
    # After merges, component ids might be sparse; compress them
    uniq_ids <- sort(unique(new_membership))
    remap <- setNames(seq_along(uniq_ids), uniq_ids)
    spatial_cluster_merged <- unname(remap[as.character(new_membership)])
    data$spatial_cluster_merged <- spatial_cluster_merged
    
    # Optional: provide a compact summary post-merge
    k_merged <- length(uniq_ids)
    
    return(list(
      mst = mst_g,
      data = data,
      mst_same = mst_same,
      components = comp$no,
      merged = TRUE,
      spatial_cluster_merged = spatial_cluster_merged,
      k_merged = k_merged,
      comp_summary_before = comp_df
    ))
  } else {
    # Nothing to merge
    data$spatial_cluster_merged <- spatial_cluster_orig
    return(list(
      mst = mst_g,
      data = data,
      mst_same = mst_same,
      components = comp$no,
      merged = FALSE,
      spatial_cluster_merged = spatial_cluster_orig,
      k_merged = comp$no,
      comp_summary_before = comp_df
    ))
  }
}


# ------------------------- Helper: merge separated mst components ----------------------------
merge_all_components_to_giant <- function(g, coords, verbose = TRUE) {
  stopifnot(inherits(g, "igraph"))
  if (is.data.frame(coords)) coords <- as.matrix(coords)
  stopifnot(nrow(coords) == vcount(g), ncol(coords) >= 2)
  
  # fast squared distance matrix helper
  sq_dists <- function(A, B) {
    AA <- rowSums(A*A); BB <- rowSums(B*B)
    outer(AA, BB, "+") - 2*(A %*% t(B))
  }
  
  iter <- 0L
  repeat {
    comp <- components(g)
    k <- comp$no
    if (verbose) cat(sprintf("Iteration %d: components = %d\n", iter, k))
    if (k <= 1L) break
    
    memb <- comp$membership
    sizes <- tabulate(memb, nbins = k)
    giant_id <- which.max(sizes)
    
    giant_idx <- which(memb == giant_id)
    A <- coords[giant_idx, 1:2, drop = FALSE]
    
    # for each other component, compute its min distance to the giant
    other_ids <- setdiff(seq_len(k), giant_id)
    best_comp <- NA_integer_
    best_pair <- c(NA_integer_, NA_integer_)
    best_d2 <- Inf
    
    for (cid in other_ids) {
      comp_idx <- which(memb == cid)
      B <- coords[comp_idx, 1:2, drop = FALSE]
      D2 <- sq_dists(A, B)
      wmin <- which.min(D2)
      ii <- ((wmin - 1L) %% nrow(D2)) + 1L
      jj <- ((wmin - 1L) %/% nrow(D2)) + 1L
      d2 <- D2[ii, jj]
      if (d2 < best_d2) {
        best_d2  <- d2
        best_comp <- cid
        best_pair <- c(giant_idx[ii], comp_idx[jj])  # (u in giant, v in cid)
      }
    }
    
    u <- best_pair[1]; v <- best_pair[2]
    if (!are_adjacent(g, u, v)) {
      g <- add_edges(g, c(u, v))
      E(g)$merged_edge <- FALSE
      E(g)[ecount(g)]$merged_edge <- TRUE
      if (verbose) cat(sprintf("  + connected comp %d via (%d <-> %d), d=%.3f\n",
                               best_comp, u, v, sqrt(best_d2)))
    } else if (verbose) {
      cat("  (skipped: already connected)\n")
    }
    
    iter <- iter + 1L
  }
  
  E(g)$eid <- seq_len(ecount(g))
  g
}

## 4) Per-shard preprocessing -------------------------------------------------
preprocess_shard <- function(idx, shard_id) {
  
  y_shard   <- Y_sample[idx, , drop = FALSE]
  loc_shard <- loc[idx, , drop = FALSE]
  counts_shard <- t(counts_keep)[idx, , drop = FALSE]
  
  coords <- apply(loc_shard, 2, scale)             # z-score coords
  Y_std  <- apply(y_shard, 2, scale)        # z-score expr
  
  n_s <- nrow(y_shard)
  p_s <- ncol(y_shard)
  
  bnd_scaled <- data.frame(
    x = (bnd$x - mean(loc_shard[, c("x")])) / sd(loc_shard[, c("x")]),
    y = (bnd$y - mean(loc_shard[, c("y")])) / sd(loc_shard[, c("y")])
  )
  
  # Run BASS on shard
  res_bass <- run_bass_once_sce(counts = t(counts_shard), 
                                coords = loc_shard, 
                                pcs = y_shard)$labels_tbl
  
  # Inject scaled coords for MST helper
  res_bass$x_sc <- coords[, 1]
  res_bass$y_sc <- coords[, 2]
  
  # Penalized MST and spatial components
  init_partition <- compute_MST_spatial(
    data        = res_bass,
    coords_cols = c("x_sc","y_sc"),
    cluster_col = "label",
    bnd         = bnd_scaled,
    threshold   = 5000,
    penalty     = NULL,
    min_comp_size = 5
  )
  
  k_max_s <- length(unique(init_partition$spatial_cluster_merged))
  cat(sprintf("Spatial components (after drop): %d\n", k_max_s))
  
  # Base constrained graph (weights removed)
  mesh_shard   <- gen2dMesh(coords, bnd_scaled)
  g_shard    <- constrainedDentri(n_s, mesh_shard)
  # Ensure vertex/edge ids present
  E(g_shard)$eid <- seq_len(ecount(g_shard))
  V(g_shard)$vid <- seq_len(vcount(g_shard))
  
  g_merged <- merge_all_components_to_giant(g_shard, coords, verbose = TRUE)
  cat(sprintf("Final components: %d\n", components(g_merged)$no))
  
  # Extract edgelist of merged graph
  el <- as_edgelist(g_merged, names = FALSE)
  
  # Compute edge lengths exactly like in constrainedDentri()
  distance <- sqrt(rowSums((coords[el[, 1], ] - coords[el[, 2], ])^2))
  
  # Assign back to the graph
  E(g_merged)$weight <- distance
  
  # MST (strip weight after MST)
  g_mst   <- mst(g_merged)
  g_mst   <- delete_edge_attr(g_mst, "weight")
  g_merged <- delete_edge_attr(g_merged, "weight")
  
  # Required by downstream: keep a copy named exactly graph0
  graph0 <- g_merged
  mstgraph0 <- g_mst
  
  cat(sprintf("Base graph components: %d\n", igraph::components(graph0)$no))
  cat(sprintf("Base MST components: %d\n", igraph::components(mstgraph0)$no))
  
  # Cluster seeds
  clust_s <- as.integer(init_partition$spatial_cluster_merged)
  
  # Per-component means
  cluster_means_matrix <- matrix(0, nrow = k_max_s, ncol = p)
  for (i in seq_len(k_max_s)) {
    idx_i <- which(clust_s == i)
    cluster_means_matrix[i, ] <- colMeans(Y_std[idx_i, , drop = FALSE])
  }
  
  # Init lists for M temperatures
  cluster_mat <- matrix(rep(clust_s, times = M), nrow = n, ncol = M, byrow = FALSE)
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
    lambda_s   = diag(1, p),
    nu         = p,
    M          = M,
    k_max      = k_max_s,
    lambda_k   = 15
  )
  
  ## Save bundle
  save(
    list = c("y_shard", "loc_shard", "graph0",
             "init_val", "hyperpar", "clust_s", "idx_global"),
    file = file.path(out_dir, sprintf("Shard_%03d_Preproc_BAST.RData", shard_id))
  )
  
  cat(sprintf("[OK] Shard %03d → %s\n",
              shard_id,
              file.path(out_dir, sprintf("Shard_%03d_Preproc_BAST.RData", shard_id))))
}

## 6) Run all shards ----------------------------------------------------------
for (sid in seq_along(shards)) {
  idx_global <- shards[[sid]]
  preprocess_shard(idx_global, sid)
}

cat("✅  All shard preprocessing files written to: ", out_dir, "\n")