###############################################################################
##  POSTPROC for *STARmap* — now targets the 4 "Full_*" files (3PCs & 10PCs)
##  Preproc loads: Y_std, coords_sc, graph0, init_val, hyperpar, clust_s
###############################################################################

my_lib_path <- "/scratch/user/varogovchenko/Rlibs"
.libPaths(my_lib_path)

library(Matrix)

base_dir <- "/scratch/user/varogovchenko/BASTION_HPRC/STARmap"

# ---------- helper to pick best labels (unchanged) ----------
pick_best_labels <- function(res_obj) {
  groups_assign <- res_obj$cluster_out
  k_out         <- res_obj$k_out
  trees         <- res_obj$tree_out
  M             <- length(trees[[1]])
  n_s           <- ncol(groups_assign)
  
  modeK      <- as.numeric(names(which.max(table(k_out))))
  iter_modeK <- which(k_out == modeK)
  
  W_cum <- matrix(0, n_s, n_s)
  for (it in iter_modeK) {
    X     <- table(sequence(n_s), groups_assign[it, ])
    W_cum <- W_cum + X %*% t(X)
  }
  W_bar <- W_cum / length(iter_modeK)
  
  fn <- numeric(length(iter_modeK))
  for (j in seq_along(iter_modeK)) {
    X      <- table(sequence(n_s), groups_assign[iter_modeK[j], ])
    fn[j]  <- norm(X %*% t(X) - W_bar, type = "F")
  }
  best_iter <- iter_modeK[which.min(fn)]
  
  list(cluster = groups_assign[best_iter, ],
       trees   = trees[[best_iter]],
       M       = M)
}

# ---------- process the 3PCs & 10PCs "Full_*" files ----------
for (pc in c(3, 10)) {
  message("==> Processing STARmap Full_* for ", pc, " PCs")
  
  # (1) Load result
  res_file <- file.path(base_dir, sprintf("Full_Result_BAST_%dPCs.RData", pc))
  if (!file.exists(res_file)) {
    warning("Result file missing: ", basename(res_file), " — skipping.")
    next
  }
  load(res_file)  # loads: result
  if (!exists("result")) stop("Object 'result' not found in ", res_file)
  best <- pick_best_labels(result)
  rm(result)
  
  # (2) Load preproc (Y_std, coords_sc, graph0, init_val, hyperpar, clust_s)
  pre_file <- file.path(base_dir, sprintf("Full_Preproc_BAST_%dPCs.RData", pc))
  if (!file.exists(pre_file)) stop("Preproc file missing: ", basename(pre_file))
  load(pre_file)
  
  # ---- pick the matrix actually present ----
  if (exists("Y_std")) {
    y_mat <- Y_std
  } else if (exists("y_shard")) {
    y_mat <- y_shard
  } else {
    stop("No expression matrix found (expected Y_std or y_shard) in ", pre_file)
  }
  
  # ---- labels: prefer existing clust_s; else best from result ----
  clust_vec <- if (exists("clust_s")) clust_s else best$cluster
  if (length(clust_vec) != nrow(y_mat)) {
    stop("Length of clust_s (", length(clust_vec),
         ") != nrow(y_mat) (", nrow(y_mat), ") in ", basename(pre_file))
  }
  
  k_max_s  <- length(unique(clust_vec))
  M        <- best$M
  temp_lad <- seq(0.1, 1, length.out = M)
  
  # ---- rebuild init_val slots using y_mat & clust_vec ----
  cluster_mat <- matrix(rep(clust_vec, M), ncol = M)
  
  mu_list      <- vector("list", M)
  mu_teams     <- vector("list", M)
  sigmasq_y    <- vector("list", M)
  mstgraph_lst <- vector("list", M)
  
  uniq_k <- sort(unique(clust_vec))
  for (m in seq_len(M)) {
    mu_list[[m]] <- t(vapply(
      uniq_k,
      function(k) colMeans(y_mat[clust_vec == k, , drop = FALSE]),
      numeric(ncol(y_mat))
    ))
    mu_teams[[m]]     <- colMeans(y_mat)
    sigmasq_y[[m]]    <- cov(y_mat)
    mstgraph_lst[[m]] <- best$trees
  }
  
  init_val$cluster      <- cluster_mat
  init_val$mu           <- mu_list
  init_val$mu_teams     <- mu_teams
  init_val$sigmasq_y    <- sigmasq_y
  init_val$mstgraph_lst <- mstgraph_lst
  init_val$teams        <- matrix(1, nrow = k_max_s, ncol = M)
  
  hyperpar$k_max      <- k_max_s
  hyperpar$j_max      <- k_max_s
  hyperpar$alpha      <- 0.5
  hyperpar$sigmasq_mu <- 0.3
  hyperpar$temp       <- temp_lad
  hyperpar$M          <- M
  
  # (3) Save Full PostProc
  post_file <- file.path(base_dir, sprintf("Full_Postproc_BAST_%dPCs.RData", pc))
  
  save_list <- intersect(
    c("Y_std","coords_sc","graph0","init_val","hyperpar","clust_s","idx_global","anchor_local"),
    ls(envir = .GlobalEnv)
  )
  # Ensure clust_s exists under that name
  if (!("clust_s" %in% save_list)) {
    clust_s <- clust_vec
    save_list <- c(save_list, "clust_s")
  }
  
  save(list = save_list, file = post_file)
  message("✔ Saved: ", basename(post_file), " [", paste(save_list, collapse = ", "), "]")
}

message("🏁 Done: STARmap Full_* files processed.")