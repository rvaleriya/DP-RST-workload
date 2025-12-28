###############################################################################
## (1) Load libraries and set paths
###############################################################################
my_lib_path <- "/scratch/user/varogovchenko/Rlibs"
.libPaths(my_lib_path)

library(Matrix)
library(future)
library(furrr)
library(mclust) 

res_dir  <- "/scratch/user/varogovchenko/BASTION_HPRC/Visium_HD_Human_Colon_Cancer/BAST_Results_10PC"
pre_dir  <- "/scratch/user/varogovchenko/BASTION_HPRC/Visium_HD_Human_Colon_Cancer/BAST_Shards_10PC"
post_dir <- "/scratch/user/varogovchenko/BASTION_HPRC/Visium_HD_Human_Colon_Cancer/Shards_PostBAST_10PC"
dir.create(post_dir, showWarnings = FALSE, recursive = TRUE)

###############################################################################
## (2) Helper: pick best labels from result object
###############################################################################
pick_best_labels <- function(res_obj) {
  groups_assign <- res_obj$cluster_out
  k_out         <- res_obj$k_out
  trees         <- res_obj$tree_out
  M             <- length(trees[[1]])
  
  n_s <- ncol(groups_assign)
  tbl         <- table(k_out)
  modeK       <- as.numeric(names(tbl)[tbl == max(tbl)])
  iter_modeK  <- which(k_out == modeK)
  
  W_cum <- matrix(0, n_s, n_s)
  for (it in iter_modeK) {
    X      <- table(sequence(n_s), groups_assign[it, ])
    W_cum  <- W_cum + X %*% t(X)
  }
  W_bar <- W_cum / length(iter_modeK)
  
  fn <- numeric(length(iter_modeK))
  for (j in seq_along(iter_modeK)) {
    X   <- table(sequence(n_s), groups_assign[iter_modeK[j], ])
    fn[j] <- norm(X %*% t(X) - W_bar, type = "F")
  }
  best_iter <- iter_modeK[ which.min(fn) ]
  
  list(cluster = groups_assign[best_iter, ],
       trees   = trees[[best_iter]],
       M       = M)
}

###############################################################################
## (3) List result files
###############################################################################
res_files <- list.files(
  res_dir,
  pattern = "^Shard_\\d{3}_Result_BAST\\.RData$",
  full.names = TRUE
)
cat("[INFO] Found", length(res_files), "shard result files\n")

###############################################################################
## (4) Worker function: process one shard
###############################################################################
process_shard <- function(rf) {
  sid <- sub("^Shard_(\\d{3}).*", "\\1", basename(rf))
  
  # --- Load BAST result and pick best labels ---
  load(rf)  # object: result
  best <- pick_best_labels(result)
  rm(result)
  
  # --- Load preproc file ---
  pre_file <- file.path(pre_dir, sprintf("Shard_%s_Preproc_BAST.RData", sid))
  if (!file.exists(pre_file))
    stop("Preproc file missing for shard ", sid)
  load(pre_file)  # loads y_shard, loc_shard, init_val, hyperpar, clust_s, idx_global, graph0, ann_shard
  
  # --- Update cluster info ---
  clust_s <- best$cluster
  k_max_s <- length(unique(clust_s))
  M       <- best$M
  temp_lad <- seq(0.1, 1, length.out = M)
  
  ari_val <- adjustedRandIndex(factor(ann_shard), factor(clust_s))
  cat(sprintf("[Shard %s] ARI = %.4f\n", sid, ari_val))
  
  cluster_mat <- matrix(rep(clust_s, M), ncol = M)
  mu_list   <- vector("list", M)
  mu_teams  <- vector("list", M)
  sigmasq_y <- vector("list", M)
  mstgraph_lst <- vector("list", M)
  
  for (m in seq_len(M)) {
    mu_list[[m]] <- t(vapply(sort(unique(clust_s)),
                             \(k) colMeans(y_shard[clust_s == k, , drop = FALSE]),
                             numeric(ncol(y_shard))))
    mu_teams[[m]]  <- colMeans(y_shard)
    sigmasq_y[[m]] <- cov(y_shard)
    mstgraph_lst[[m]] <- best$trees
  }
  
  init_val$cluster        <- cluster_mat
  init_val$mu             <- mu_list
  init_val$mu_teams       <- mu_teams
  init_val$sigmasq_y      <- sigmasq_y
  init_val$mstgraph_lst   <- mstgraph_lst
  init_val$teams          <- matrix(1, nrow = k_max_s, ncol = M)
  init_val[["trees"]] <- NULL
  
  hyperpar$k_max <- hyperpar$j_max <- k_max_s
  hyperpar$alpha <- 0.5
  hyperpar$sigmasq_mu <- 0.3
  hyperpar$temp  <- temp_lad
  hyperpar$M     <- M
  
  # --- Save postprocessed bundle ---
  post_file <- file.path(post_dir, sprintf("Shard_%s_PostBAST.RData", sid))
  save(list = c("y_shard", "loc_shard", "graph0",
                "init_val", "hyperpar",
                "clust_s", "idx_global", "ann_shard"),
       file = post_file)
  
  message("✔ Shard ", sid, " → ", basename(post_file))
}

###############################################################################
## (5) Parallel backend
###############################################################################
slurm_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
num_workers <- if (!is.na(slurm_cores)) slurm_cores else future::availableCores()

# Use forked workers (multicore) and avoid global export
plan(multicore, workers = num_workers)
options(future.globals.maxSize = Inf)
Sys.setenv(OMP_NUM_THREADS = "1")

opts <- furrr_options(seed = TRUE, scheduling = Inf, stdout = TRUE, globals = FALSE)
cat("[INFO] Using", num_workers, "parallel workers\n")

###############################################################################
## (6) Launch processing in parallel
###############################################################################
t0 <- Sys.time()
future_walk(res_files, process_shard, .options = opts)

cat("🏁  All shard bundles updated at", format(Sys.time()), "\n")
cat("⏱  Total runtime (min):", round(difftime(Sys.time(), t0, units = "mins"), 2), "\n")