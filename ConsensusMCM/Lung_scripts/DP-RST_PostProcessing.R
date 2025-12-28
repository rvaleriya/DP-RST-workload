#!/usr/bin/env Rscript

# ===========================================
# Parallel DP-RST post-processing for shards
# ===========================================

t0 <- Sys.time()

# ---- Library path, RNG, and packages ----
my_lib_path <- "/scratch/user/varogovchenko/Rlibs"
.libPaths(my_lib_path)

set.seed(2025)

library(DP.RST)
library(mclust)
library(ggplot2)
library(future)
library(furrr)
library(tools)

# ---- Base directories (given) ----
base_dir     <- "/scratch/user/varogovchenko/BASTION_HPRC/Lung_Xenium"
res_dir_3    <- file.path(base_dir, "DP-RST_Results_3PC")
res_dir_10   <- file.path(base_dir, "DP-RST_Results_10PC")
shards_dir_3 <- file.path(base_dir, "Shards_PostBAST_3PC")
shards_dir_10<- file.path(base_dir, "Shards_PostBAST_10PC")
out_root     <- file.path(base_dir, "DP-RST_PostProcessed_Results")

# ---- Helpers ----
safe_load_as_list <- function(rdata_path) {
  e <- new.env(parent = emptyenv())
  load(rdata_path, envir = e)
  as.list(e)
}

ensure_dirs <- function(paths) for (p in paths) if (!dir.exists(p)) dir.create(p, recursive = TRUE)

extract_shard_id <- function(filename) {
  m <- regmatches(filename, regexec("Shard[_-](\\d+)", filename))[[1]]
  if (length(m) >= 2) as.integer(m[2]) else NA_integer_
}

# ---- Parallel config (single node, multicore) ----
slurm_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
num_workers <- if (!is.na(slurm_cores)) slurm_cores else future::availableCores()
plan(multicore, workers = num_workers)
options(mc.cores = num_workers)
Sys.setenv(OMP_NUM_THREADS = "1")  # avoid nested OpenMP contention
cat(sprintf("[INFO] Using %d workers\n", num_workers))

# --------------------------------------------
# Core per-PC function: parallel post-processing
# --------------------------------------------
process_pc_parallel <- function(pc, res_dir, shards_dir, out_root_dir) {
  cat("\n==============================\n")
  cat(sprintf("Processing PC = %d\n", pc))
  cat(  "==============================\n")
  
  # Output folders
  out_pc_dir    <- file.path(out_root_dir, paste0("PC", pc))
  out_prop_dir  <- file.path(out_pc_dir, "proportions")
  out_means_dir <- file.path(out_pc_dir, "means")
  out_cov_dir <- file.path(out_pc_dir, "covariance")
  out_lists_dir <- file.path(out_pc_dir, "lists")
  out_figs_dir  <- file.path(out_pc_dir, "figs")
  ensure_dirs(c(out_pc_dir, out_prop_dir, out_means_dir, out_cov_dir, out_lists_dir, out_figs_dir))
  
  # Enumerate DP-RST result files
  res_files <- list.files(res_dir, pattern = "^Shard_\\d+_Result_DPRST\\.RData$", full.names = TRUE)
  if (length(res_files) == 0L) stop("No DP-RST result files in: ", res_dir)
  shard_ids <- vapply(basename(res_files), extract_shard_id, integer(1))
  ord <- order(shard_ids)
  res_files <- res_files[ord]
  shard_ids <- shard_ids[ord]
  
  # Parallel worker: returns small summary; writes heavy outputs to disk
  worker_fun <- function(res_path, sid) {
    # Load DP-RST 'result'
    lr <- safe_load_as_list(res_path)
    if (!("result" %in% names(lr))) stop("Object 'result' missing in: ", res_path)
    result <- lr$result
    if (!is.list(result)) stop("'result' is not a list: ", res_path)
    
    # Partition (mode-based, batched)
    best <- partition(result, method = "mode_based", batch_size = 500)
    final_clusters <- best$teams_partition[best$groups_partition]
    
    # Proportions CSV
    tab <- table(final_clusters)
    prop_df <- data.frame(
      cluster     = as.integer(names(tab)),
      count       = as.integer(tab),
      proportion  = as.numeric(tab) / length(final_clusters),
      total_n     = length(final_clusters)
    )
    write.csv(prop_df,
              file.path(out_prop_dir, sprintf("Shard_%03d_proportions.csv", sid)),
              row.names = FALSE)
    
    # Selected-iteration team means CSV
    sel_idx <- best$selected_iteration
    mu_mat  <- result$mu_teams_out[[sel_idx]]
    if (is.null(mu_mat)) stop("mu_teams_out[[selected_iteration]] is NULL for shard ", sid)
    mu_df <- as.data.frame(mu_mat)
    mu_df$team <- seq_len(nrow(mu_df))
    mu_df <- mu_df[, c("team", setdiff(names(mu_df), "team"))]
    write.csv(mu_df,
              file.path(out_means_dir, sprintf("Shard_%03d_means.csv", sid)),
              row.names = FALSE)
              
    # Selected-iteration covariance matrixes CSV
    cov_mat  <- result$sigmasq_y_out[[sel_idx]]
    if (is.null(cov_mat)) stop("sigmasq_y_out[[selected_iteration]] is NULL for shard ", sid)
    cov_df <- as.data.frame(cov_mat)
    write.csv(cov_df,
              file.path(out_cov_dir, sprintf("Shard_%03d_covariance.csv", sid)),
              row.names = FALSE)
    
    # Load PostBAST shard to get labels and meta
    post_path <- file.path(shards_dir, sprintf("Shard_%03d_PostBAST.RData", sid))
    if (!file.exists(post_path)) stop("Missing PostBAST shard: ", post_path)
    lp <- safe_load_as_list(post_path)
    
    if (!("ann_shard" %in% names(lp))) stop("ann_shard missing for shard ", sid)
    if (!("idx_global" %in% names(lp))) stop("idx_global missing for shard ", sid)
    if (!("init_val" %in% names(lp))) stop("init_val missing for shard ", sid)
    if (!all(c("loc_shard", "y_shard") %in% names(lp)))
      stop("init_val$loc_shard or init_val$y_shard missing for shard ", sid)
    
    ann_shard <- lp$ann_shard
    idx_global<- lp$idx_global
    loc_shard <- lp$loc_shard
    y_shard   <- lp$y_shard
    
    if (length(final_clusters) != length(ann_shard))
      stop(sprintf("Length mismatch (final_clusters=%d, ann_shard=%d) for shard %d",
                   length(final_clusters), length(ann_shard), sid))
    
    # ARI
    keep <- !is.na(ann_shard)
    truth_int <- as.integer(factor(ann_shard[keep]))
  pred_int  <- as.integer(final_clusters[keep])
    ari_val <- mclust::adjustedRandIndex(truth_int, pred_int)
    
    # Per-shard list (save as RDS to avoid big return values)
    shard_list <- list(
      shard_id       = sid,
      final_clusters = as.integer(final_clusters),
      loc_shard      = loc_shard,
      y_shard        = y_shard,
      ann_shard      = ann_shard,
      idx_global     = idx_global,
      mu_team_means  = mu_mat,
      selected_iter  = sel_idx
    )
    saveRDS(shard_list, file.path(out_lists_dir, sprintf("Shard_%03d_list.rds", sid)))
    
    # Small summary returned to master
    data.frame(shard = sid, ARI = ari_val)
  }
  
  # Run in parallel
  opts <- furrr_options(seed = TRUE, scheduling = Inf, stdout = TRUE, globals = FALSE)
  res_df_list <- future_map2(
    .x = res_files,
    .y = shard_ids,
    .f = worker_fun,
    .options = opts
  )
  
  # Combine ARIs
  ari_df <- do.call(rbind, res_df_list)
  ari_df <- ari_df[order(ari_df$shard), ]
  ari_csv_path <- file.path(out_pc_dir, sprintf("ARI_by_shard_PC%d.csv", pc))
  write.csv(ari_df, ari_csv_path, row.names = FALSE)
  
  # Histogram
  p <- ggplot(ari_df, aes(x = ARI)) +
    geom_histogram(bins = 20) +
    labs(title = paste0("ARI across shards (PC=", pc, ")"),
         x = "Adjusted Rand Index", y = "Count") +
    theme_minimal(base_size = 12)
  hist_path <- file.path(out_figs_dir, sprintf("ARI_hist_PC%d.png", pc))
  ggsave(hist_path, plot = p, width = 6.5, height = 4.0, dpi = 300)
  
  # Build combined list-of-lists from per-shard RDS files
  rds_files <- list.files(out_lists_dir, pattern = "^Shard_\\d+_list\\.rds$", full.names = TRUE)
  if (length(rds_files) != length(shard_ids)) {
    warning("Number of shard RDS files (", length(rds_files),
            ") != number of shards (", length(shard_ids), ").")
  }
  # Keep order by shard id
  rids <- vapply(basename(rds_files), extract_shard_id, integer(1))
  ord2 <- order(rids)
  rds_files <- rds_files[ord2]
  rids      <- rids[ord2]
  
  shard_lists <- vector("list", length(rds_files))
  names(shard_lists) <- paste0("Shard_", rids)
  for (i in seq_along(rds_files)) shard_lists[[i]] <- readRDS(rds_files[i])
  
  combined_path <- file.path(out_lists_dir, sprintf("Combined_Shards_PC%d.RData", pc))
  save(shard_lists, file = combined_path)
  
  # Console summary
  cat(sprintf("[PC=%d] Shards processed: %d\n", pc, nrow(ari_df)))
  cat(sprintf("[PC=%d] ARI mean=%.6f, median=%.6f, min=%.6f, max=%.6f\n",
              pc, mean(ari_df$ARI), median(ari_df$ARI),
              min(ari_df$ARI), max(ari_df$ARI)))
  cat("[OUT] Proportions CSVs : ", out_prop_dir, "\n")
  cat("[OUT] Means CSVs       : ", out_means_dir, "\n")
  cat("[OUT] ARI CSV          : ", ari_csv_path, "\n")
  cat("[OUT] ARI histogram    : ", hist_path, "\n")
  cat("[OUT] Combined .RData  : ", combined_path, "\n")
}

# ---- Execute for PC=3 and PC=10 ----
ensure_dirs(out_root)
process_pc_parallel(pc = 3,
                    res_dir = res_dir_3,
                    shards_dir = shards_dir_3,
                    out_root_dir = out_root)

process_pc_parallel(pc = 10,
                    res_dir = res_dir_10,
                    shards_dir = shards_dir_10,
                    out_root_dir = out_root)

cat("\n✅ All post-processing completed in ",
    round(difftime(Sys.time(), t0, units = "mins"), 2), " minutes.\n", sep = "")