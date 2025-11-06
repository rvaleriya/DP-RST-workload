plot_mst_clusters_pdf <- function(coords, graph, clusters, out_pdf = "MST_overlay_init_clusters.pdf") {
  
  # Load required packages
  library(igraph)
  library(ggplot2)
  library(DP.RST)
  
  # --- standardize coordinates ---
  if (is.matrix(coords)) {
    coords_df <- data.frame(x = coords[, 1], y = coords[, 2])
  } else if (is.data.frame(coords)) {
    xcol <- if ("x" %in% names(coords)) "x" else names(coords)[1]
    ycol <- if ("y" %in% names(coords)) "y" else names(coords)[2]
    coords_df <- data.frame(x = coords[[xcol]], y = coords[[ycol]])
  } else {
    stop("coords must be a matrix or data.frame.")
  }
  
  # --- checks ---
  stopifnot(inherits(graph, "igraph"))
  stopifnot(vcount(graph) == nrow(coords_df))
  stopifnot(length(clusters) == nrow(coords_df))
  
  cl_fac <- as.factor(clusters)
  cluster_ids <- levels(cl_fac)
  cat(sprintf("Found %d clusters.\n", length(cluster_ids)))
  
  # --- base MST plot ---
  base_mst <- DP.RST::plotGraph(coords = coords_df, graph = graph) +
    labs(x = "x", y = "y") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank())
  
  bg_df <- coords_df
  
  # --- output PDF ---
  pdf(out_pdf, width = 7, height = 6)
  for (cid in cluster_ids) {
    sel <- which(cl_fac == cid)
    fg_df <- coords_df[sel, , drop = FALSE]
    p <- base_mst +
      geom_point(data = bg_df, aes(x = x, y = y), size = 0.15, alpha = 0.15) +
      geom_point(data = fg_df, aes(x = x, y = y), size = 0.6, alpha = 0.95) +
      ggtitle(sprintf("Cluster %s | n = %d", cid, nrow(fg_df)))
    print(p)
  }
  dev.off()
  
  cat(sprintf("Saved %d pages to %s\n", length(cluster_ids), normalizePath(out_pdf)))
}