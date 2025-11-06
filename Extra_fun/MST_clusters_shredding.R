# Helper: merge tiny spatial components strictly along MST edges
merge_small_by_mst <- function(mst_g,
                               init_membership,
                               min_comp_size,
                               weight_attr = "raw_len",
                               verbose = FALSE) {
  stopifnot(igraph::is_igraph(mst_g),
            length(init_membership) == igraph::vcount(mst_g),
            is.numeric(min_comp_size), min_comp_size >= 1,
            weight_attr %in% igraph::edge_attr_names(mst_g))
  
  membership <- as.integer(init_membership)
  
  build_comp_edges <- function() {
    el <- igraph::as_edgelist(mst_g, names = FALSE)
    cu <- membership[el[, 1]]
    cv <- membership[el[, 2]]
    w  <- igraph::edge_attr(mst_g, weight_attr)
    data.frame(u = cu, v = cv, w = w, stringsAsFactors = FALSE)
  }
  
  repeat {
    sizes  <- tabulate(membership, nbins = max(membership))
    small  <- which(sizes > 0L & sizes < min_comp_size)
    if (!length(small)) break
    
    ce <- build_comp_edges()
    ce <- ce[ce$u != ce$v, , drop = FALSE]
    if (!nrow(ce)) break
    
    merged_any <- FALSE
    order_small <- small[order(sizes[small])]
    for (s in order_small) {
      ce_s <- ce[(ce$u == s) | (ce$v == s), , drop = FALSE]
      if (!nrow(ce_s)) next
      tgt_vec <- ifelse(ce_s$u == s, ce_s$v, ce_s$u)
      tgt <- tgt_vec[which.min(ce_s$w)]
      if (!is.na(tgt) && tgt != s) {
        membership[membership == s] <- tgt
        merged_any <- TRUE
        if (verbose) message(sprintf("[merge_small_by_mst] %d -> %d", s, tgt))
      }
    }
    membership <- as.integer(factor(membership, levels = sort(unique(membership))))
    if (!merged_any) break
  }
  
  as.integer(factor(membership, levels = sort(unique(membership))))
}