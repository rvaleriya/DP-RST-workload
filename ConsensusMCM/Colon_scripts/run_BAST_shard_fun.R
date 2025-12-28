run_shard <- function(file_path) {
  # Helper to extract numeric shard ID from filename
  get_id <- function(path)
    as.integer(sub("^.*Shard_([0-9]+)_Preproc_BAST\\.RData$", "\\1", basename(path)))

  shard_id <- get_id(file_path)
  t0 <- Sys.time()
  message(sprintf("\n[Shard %03d | %s] START", shard_id, format(t0)))

  # --- Load data ---
  message(sprintf("[Shard %03d] Loading objects from %s", shard_id, file_path))
  load(file_path)
  t_load <- Sys.time()
  message(sprintf("[Shard %03d] Loaded in %.2f sec", shard_id,
                  as.numeric(difftime(t_load, t0, units = "secs"))))
  flush.console()

  # --- Run fitBAST ---
  message(sprintf("[Shard %03d] Running fitBAST...", shard_id))
  t_bast_start <- Sys.time()
  result <- fitBAST(
    y_shard, graph0, init_val, hyperpar,
    MCMC = 30000, BURNIN = 25000, THIN = 5,
    temp = seq(0.1, 1, 0.1), PT = TRUE,
    seed = 100 + shard_id,
    backup_d = sprintf("BAST_Shards_10PC_Backup/Shard_%03d_backup.RData", shard_id)
  )
  t_bast_end <- Sys.time()
  message(sprintf("[Shard %03d] fitBAST finished in %.2f min",
                  shard_id,
                  as.numeric(difftime(t_bast_end, t_bast_start, units = "mins"))))
  flush.console()

  # --- Save result ---
  results_dir <- "/scratch/user/varogovchenko/BASTION_HPRC/Visium_HD_Human_Colon_Cancer/BAST_Results_10PC"
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(results_dir, sprintf("Shard_%03d_Result_BAST.RData", shard_id))
  save(result, file = out_file)

  t_end <- Sys.time()
  message(sprintf("[Shard %03d] Saved → %s", shard_id, out_file))
  message(sprintf("[Shard %03d | TOTAL RUNTIME] %.2f min\n",
                  shard_id,
                  as.numeric(difftime(t_end, t0, units = "mins"))))
  flush.console()
}