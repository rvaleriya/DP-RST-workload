merge_all_components_to_giant <- function(g, coords, verbose = TRUE, max_iter = 10000L) {
  stopifnot(igraph::is_igraph(g))
  if (is.data.frame(coords)) coords <- as.matrix(coords)
  stopifnot(nrow(coords) == igraph::vcount(g), ncol(coords) >= 2)
  
  sq_dists <- function(A, B) {
    AA <- rowSums(A*A); BB <- rowSums(B*B)
    outer(AA, BB, "+") - 2*(A %*% t(B))
  }
  
  iter <- 0L
  repeat {
    comp <- igraph::components(g); k <- comp$no
    if (verbose) cat(sprintf("Iteration %d: components = %d\n", iter, k))
    if (k <= 1L) break
    if (iter >= max_iter) {
      warning(sprintf("[merge_all_components_to_giant] max_iter=%d reached; stopping early.", max_iter))
      break
    }
    
    memb <- comp$membership; sizes <- tabulate(memb, nbins = k)
    giant_id  <- which.max(sizes)
    giant_idx <- which(memb == giant_id)
    other_ids <- setdiff(seq_len(k), giant_id)
    
    A <- coords[giant_idx, 1:2, drop = FALSE]
    best_pair <- c(NA_integer_, NA_integer_); best_d2 <- Inf
    for (cid in other_ids) {
      comp_idx <- which(memb == cid)
      B <- coords[comp_idx, 1:2, drop = FALSE]
      D2 <- sq_dists(A, B)
      wmin <- which.min(D2)
      ii <- ((wmin - 1L) %% nrow(D2)) + 1L
      jj <- ((wmin - 1L) %/% nrow(D2)) + 1L
      d2 <- D2[ii, jj]
      if (d2 < best_d2) { best_d2 <- d2; best_pair <- c(giant_idx[ii], comp_idx[jj]) }
    }
    
    u <- best_pair[1]; v <- best_pair[2]
    if (!igraph::are_adjacent(g, u, v)) {
      g <- igraph::add_edges(g, c(u, v))
      E(g)$merged_edge <- FALSE
      E(g)[igraph::ecount(g)]$merged_edge <- TRUE
      if (verbose) cat(sprintf("  + connected via (%d <-> %d), d=%.3f\n", u, v, sqrt(best_d2)))
    }
    iter <- iter + 1L
  }
  E(g)$eid <- seq_len(igraph::ecount(g))
  g
}