compute_MST_spatial <- function(data,
                                coords_cols = c("x", "y"),
                                cluster_col = "label",
                                bnd = NULL,
                                threshold = 5000,
                                penalty = NULL,
                                min_comp_size = 3) {
  
  stopifnot(all(coords_cols %in% names(data)), cluster_col %in% names(data))
  
  # --- Base graph on scaled/mesh domain ---
  coords <- as.matrix(data[, coords_cols])
  mesh   <- gen2dMesh(coords, bnd)
  graph0 <- constrainedDentri(n = nrow(coords), mesh = mesh, threshold = threshold)
  
  E(graph0)$eid <- seq_len(igraph::ecount(graph0))
  V(graph0)$vid <- seq_len(igraph::vcount(graph0))
  V(graph0)$cluster <- data[[cluster_col]]
  
  # --- Ensure connectivity BEFORE MST (this was the missing step) ---
  comp_no <- igraph::components(graph0)$no
  if (comp_no > 1L) {
    graph0 <- merge_all_components_to_giant(graph0, coords, verbose = FALSE, max_iter = 10000L)
  }
  
  # Recompute Euclidean weights for ALL edges (includes any new bridging edges)
  el <- igraph::as_edgelist(graph0, names = FALSE)
  w  <- sqrt(rowSums((coords[el[,1], , drop=FALSE] - coords[el[,2], , drop=FALSE])^2))
  E(graph0)$weight <- w
  
  # --- Spatial components inside each original label (unchanged) ---
  unique_clusters <- unique(V(graph0)$cluster)
  new_cluster_id <- 0L
  spatial_cluster_orig <- integer(nrow(data))
  for (cl in unique_clusters) {
    cluster_verts <- which(V(graph0)$cluster == cl)
    if (length(cluster_verts) == 1L) {
      new_cluster_id <- new_cluster_id + 1L
      spatial_cluster_orig[cluster_verts] <- new_cluster_id
    } else {
      subg  <- igraph::induced_subgraph(graph0, cluster_verts)
      comps <- igraph::components(subg)
      for (comp_id in seq_len(comps$no)) {
        new_cluster_id <- new_cluster_id + 1L
        idx_sub <- which(comps$membership == comp_id)
        spatial_cluster_orig[cluster_verts[idx_sub]] <- new_cluster_id
      }
    }
  }
  V(graph0)$spatial_cluster <- spatial_cluster_orig
  data$spatial_cluster <- spatial_cluster_orig
  
  # --- Penalized weights + CONNECTED MST ---
  E(graph0)$raw_len <- E(graph0)$weight
  if (is.null(penalty)) penalty <- max(E(graph0)$raw_len, na.rm = TRUE) + 1
  
  ee <- as.data.frame(igraph::as_edgelist(graph0))
  colnames(ee) <- c("V1","V2"); ee$V1 <- as.integer(ee$V1); ee$V2 <- as.integer(ee$V2)
  c1 <- spatial_cluster_orig[ee$V1]; c2 <- spatial_cluster_orig[ee$V2]
  inter <- (c1 != c2)
  E(graph0)$pen_w <- E(graph0)$raw_len + penalty * inter
  
  mst_g <- igraph::mst(graph0, weights = E(graph0)$pen_w)
  stopifnot(igraph::components(mst_g)$no == 1L)  # assert MST is a single tree now
  
  V(mst_g)$cluster <- data[[cluster_col]]
  V(mst_g)$spatial_cluster <- spatial_cluster_orig
  # Edge attributes (raw_len, pen_w, eid, etc.) are inherited
  
  # --- Step 3: Within-cluster MST (pre-merge, for reference / analysis) -------
  h <- as.integer(igraph::head_of(mst_g, E(mst_g)))
  t <- as.integer(igraph::tail_of(mst_g, E(mst_g)))
  drop_cross <- spatial_cluster_orig[h] != spatial_cluster_orig[t]
  mst_same <- igraph::delete_edges(mst_g, which(drop_cross))
  comp_pre <- igraph::components(mst_same)
  
  # --- Step 4: Merge tiny components ALONG MST EDGES ONLY ---------------------
  if (is.null(min_comp_size) || min_comp_size <= 1L) {
    data$spatial_cluster_merged <- spatial_cluster_orig
    return(list(
      mst = mst_g,
      data = data,
      mst_same = mst_same,
      components = length(unique(spatial_cluster_orig)),
      merged = FALSE,
      spatial_cluster_merged = spatial_cluster_orig,
      graph_base = graph0
    ))
  }
  
  # Merge pass
  spatial_cluster_merged <- merge_small_by_mst(
    mst_g            = mst_g,
    init_membership  = spatial_cluster_orig,
    min_comp_size    = min_comp_size,
    weight_attr      = "raw_len",
    verbose          = FALSE
  )
  
  data$spatial_cluster_merged <- spatial_cluster_merged
  V(graph0)$spatial_cluster_merged <- spatial_cluster_merged
  V(mst_g)$spatial_cluster_merged  <- spatial_cluster_merged
  
  # Within-merged-cluster MST (each merged cluster should be a connected subtree)
  h_m <- as.integer(igraph::head_of(mst_g, E(mst_g)))
  t_m <- as.integer(igraph::tail_of(mst_g, E(mst_g)))
  drop_cross_m <- spatial_cluster_merged[h_m] != spatial_cluster_merged[t_m]
  mst_same_merged <- igraph::delete_edges(mst_g, which(drop_cross_m))
  
  # --- Optional sanity check: each merged cluster is connected in mst_same_merged
  comps_after <- igraph::components(mst_same_merged)
  # Map each vertex's component id vs its assigned cluster; they must align 1-1
  # (i.e., every merged cluster occupies exactly one connected component)
  # We won't stop() here, but this is a strong diagnostic:
  # bad <- any(tapply(seq_along(spatial_cluster_merged), spatial_cluster_merged,
  #                   function(ix) length(unique(comps_after$membership[ix])) > 1))
  # if (bad) warning("Some merged clusters are not connected in within-cluster MST.")
  
  # Pre-merge summary
  comp_idx_list <- split(seq_len(nrow(data)), spatial_cluster_orig)
  comp_df_before <- data.frame(
    comp_id = as.integer(names(comp_idx_list)),
    size    = vapply(comp_idx_list, length, integer(1)),
    stringsAsFactors = FALSE
  )
  
  list(
    mst                    = mst_g,
    data                   = data,
    mst_same               = mst_same_merged,
    components             = length(unique(spatial_cluster_merged)),
    merged                 = !identical(spatial_cluster_merged, spatial_cluster_orig),
    spatial_cluster_merged = spatial_cluster_merged,
    k_merged               = length(unique(spatial_cluster_merged)),
    comp_summary_before    = comp_df_before,
    graph_base             = graph0
  )
}