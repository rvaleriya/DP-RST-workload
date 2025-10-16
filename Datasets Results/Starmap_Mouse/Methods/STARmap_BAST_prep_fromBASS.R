# ---------------------------------------------------------------------------
# BAST Preprocessing (FULL DATASET) for 3 PCs and 10 PCs
# ---------------------------------------------------------------------------

library(DP.RST)
library(MASS)
library(igraph)
library(fields)
library(mclust)
library(dplyr)
library(purrr)
library(SingleCellExperiment)

set.seed(42)

# ------------------------- Config -------------------------------------------
M       <- 10
base_in <- "~/Desktop/DP-RST-workload/Datasets Results/Starmap_Mouse"
out_dir <- file.path(base_in, "Results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pcs_path <- file.path(base_in, "starmap_pcs_coords_labels.csv")
bnd_path <- file.path(base_in, "starmap_manual_boundary.csv")

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
# ------------------------- One pass for n_pcs -------------------------------
run_once <- function(n_pcs) {
  cat(sprintf("\n===== Starting FULL preprocessing with %d PCs =====\n", n_pcs))
  
  pcs <- read.csv(pcs_path, check.names = FALSE)
  stopifnot(all(c("x","y") %in% colnames(pcs)))
  
  boundary_df <- read.csv(bnd_path, check.names = FALSE)
  stopifnot(all(c("x","y") %in% colnames(boundary_df)))
  
  # Pick matching BASS init file
  init_file <- if (n_pcs == 3) {
    file.path(base_in, "Results", "InitBAST_from_BASS_3PCs.csv")
  } else if (n_pcs == 10) {
    file.path(base_in, "Results", "InitBAST_from_BASS_10PCs.csv")
  } else {
    stop("n_pcs must be 3 or 10 for this script.")
  }
  init_values <- read.csv(init_file, check.names = FALSE)
  stopifnot("label" %in% colnames(init_values))
  
  # Build Y (first n_pcs columns assumed to be PCs), coords
  Y       <- as.matrix(pcs[, seq_len(n_pcs), drop = FALSE])
  coords  <- as.matrix(pcs[, c("x","y")])
  
  # Standardize on FULL data
  Y_std     <- scale(Y)
  print(head(Y_std))
  coords_sc <- scale(coords)
  
  # Scale boundary using full-data coords stats
  bnd_scaled <- data.frame(
    x = (boundary_df$x - mean(pcs$x)) / sd(pcs$x),
    y = (boundary_df$y - mean(pcs$y)) / sd(pcs$y)
  )
  
  # Align init labels to pcs rows if ID present in both
  if ("ID" %in% colnames(pcs) && "ID" %in% colnames(init_values)) {
    idx <- match(pcs$ID, init_values$ID)
    stopifnot(!any(is.na(idx)))
    init_values <- init_values[idx, , drop = FALSE]
  }
  
  # Inject scaled coords for MST helper
  init_values$x <- coords_sc[, 1]
  init_values$y <- coords_sc[, 2]
  
  # Penalized MST and spatial components
  init_partition <- compute_MST_spatial(
    data        = init_values,
    coords_cols = c("x","y"),
    cluster_col = "label",
    bnd         = bnd_scaled,
    threshold   = 5000,
    penalty     = NULL,
    min_comp_size = 4
  )
  
  k_max_s <- length(unique(init_partition$spatial_cluster_merged))
  cat(sprintf("Spatial components (after drop): %d\n", k_max_s))
  
  # Base constrained graph (weights removed)
  mesh   <- gen2dMesh(coords_sc, bnd_scaled)
  graph0 <- constrainedDentri(nrow(coords_sc), mesh)
  E(graph0)$eid <- seq_len(ecount(graph0))
  V(graph0)$vid <- seq_len(vcount(graph0))
  g_mst  <- mst(graph0)
  g_mst  <- delete_edge_attr(g_mst, "weight")
  graph0 <- delete_edge_attr(graph0, "weight")
  
  mstgraph0 <- delete_edge_attr(init_partition$mst, "weight")
  
  cat(sprintf("Base graph components: %d\n", igraph::components(graph0)$no))
  cat(sprintf("Penalized MST components kept: %d\n",
              igraph::components(init_partition$mst_same)$no))
  
  # Cluster seeds
  clust_s <- as.integer(init_partition$spatial_cluster_merged)
  n <- nrow(Y_std); p <- ncol(Y_std)
  
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
  
  save_path <- file.path(out_dir, sprintf("Full_Preproc_BAST_%dPCs.RData", n_pcs))
  save(list = c("Y_std","coords_sc","graph0","init_val","hyperpar","clust_s","mstgraph0"),
       file = save_path)
  cat(sprintf("Saved FULL preprocessing to: %s\n", save_path))
}

# ------------------------- Run both settings --------------------------------
run_once(3)
run_once(10)