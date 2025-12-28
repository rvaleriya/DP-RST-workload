run_shard <- function(file_path) {
  # Helper to extract numeric shard ID from filename
  get_id <- function(path)
    as.integer(sub("^.*Shard_([0-9]+)_PostBAST\\.RData$", "\\1", basename(path)))

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
  
  message(sprintf("[Shard %03d] Running DP-RST...", shard_id))
  t_dprst_start <- Sys.time()
  result <- DP.RST(
    y_shard, graph0, init_val, hyperpar,
    MCMC = 25000, BURNIN = 20000, THIN = 5,
    PT = TRUE, PT_diff = 0.1,
    seed = 100 + shard_id
  )
  t_dprst_end <- Sys.time()
  message(sprintf("[Shard %03d] DP-RST finished in %.2f min",
                  shard_id,
                  as.numeric(difftime(t_dprst_end, t_dprst_start, units = "mins"))))
  flush.console()

  results_dir <- "/scratch/user/varogovchenko/BASTION_HPRC/Visium_HD_Human_Colon_Cancer/DP-RST_Results_10PC"
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(results_dir, sprintf("Shard_%s_Result_DPRST.RData", shard_id))
  save(result, file = out_file)
  
  t_end <- Sys.time()
  message(sprintf("[Shard %03d] Saved → %s", shard_id, out_file))
  message(sprintf("[Shard %03d | TOTAL RUNTIME] %.2f min\n",
                  shard_id,
                  as.numeric(difftime(t_end, t0, units = "mins"))))
  flush.console()
}
