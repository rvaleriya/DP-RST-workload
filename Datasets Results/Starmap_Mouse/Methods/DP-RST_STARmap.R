#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Run DP.RST on STARmap *Full* postproc files (3PCs & 10PCs) in parallel
# ---------------------------------------------------------------------------
t0 <- Sys.time()
cat("========== DP.RST STARmap FULL inference started at",
    format(Sys.time()), "==========\n")

# ---------------------------------------------------------------------------
# 0. Library setup -----------------------------------------------------------
my_lib_path <- "/scratch/user/varogovchenko/Rlibs"
.libPaths(my_lib_path)

library(DP.RST)
library(future)
library(furrr)

# ---------------------------------------------------------------------------
# 1. User parameters ---------------------------------------------------------
base_dir <- "/scratch/user/varogovchenko/BASTION_HPRC/STARmap"

MCMC   <- 25000          # sampler iterations
BURNIN <- 20000
THIN   <- 5
PTdiff <- 0.1            # parallel-tempering diff

# ---------------------------------------------------------------------------
# 2. Parallel backend --------------------------------------------------------
slurm_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
num_workers <- if (!is.na(slurm_cores)) slurm_cores else future::availableCores()
plan(multicore, workers = num_workers)
opts <- furrr_options(seed = TRUE, scheduling = Inf)
cat("[INFO] Using", num_workers, "parallel workers\n")

# ---------------------------------------------------------------------------
# 3. Locate STARmap FULL postproc files --------------------------------------
# Expect exactly these two (created by your STARmap postproc step):
#   Full_Postproc_BAST_3PCs.RData
#   Full_Postproc_BAST_10PCs.RData
full_postproc_files <- file.path(
  base_dir,
  c("Full_Postproc_BAST_3PCs.RData", "Full_Postproc_BAST_10PCs.RData")
)
full_postproc_files <- full_postproc_files[file.exists(full_postproc_files)]
stopifnot(length(full_postproc_files) > 0)
cat("[INFO] Found", length(full_postproc_files), "FULL postproc file(s)\n")

# Helper: extract PCs count from filename
get_pcs <- function(path) {
  as.integer(sub("^.*_([0-9]+)PCs\\.RData$", "\\1", basename(path)))
}

# ---------------------------------------------------------------------------
# 4. Worker: load FULL postproc + run DP.RST ---------------------------------
run_full <- function(file_path) {
  pcs <- get_pcs(file_path)
  cat(sprintf("[FULL %dPCs] loading objects …\n", pcs))
  
  # Loads (per your STARmap postproc): Y_std, coords, graph0, init_val, hyperpar, clust_s (may vary)
  load(file_path)
  
  # --- Backward compatibility aliases
  if (!exists("y_full")) {
    if (!exists("Y_std")) stop("[FULL ", pcs, "PCs] Y_std missing in ", file_path)
    y_full <- Y_std
  }
  if (!exists("loc_full")) {
    if (!exists("coords_sc")) stop("[FULL ", pcs, "PCs] coords missing in ", file_path)
    loc_full <- coords_sc
  }
  if (!exists("graph0"))   stop("[FULL ", pcs, "PCs] graph0 missing in ", file_path)
  if (!exists("init_val")) stop("[FULL ", pcs, "PCs] init_val missing in ", file_path)
  if (!exists("hyperpar")) stop("[FULL ", pcs, "PCs] hyperpar missing in ", file_path)
  
  cat(sprintf("[FULL %dPCs] starting DP.RST …\n", pcs))
  result <- DP.RST(
    y_full, graph0,
    init_val, hyperpar,
    MCMC = MCMC, BURNIN = BURNIN, THIN = THIN,
    PT = TRUE, PT_diff = PTdiff,
    seed = 100 + pcs
  )
  
  out_file <- file.path(
    dirname(file_path),
    sprintf("Full_Result_DP.RST_%dPCs.RData", pcs)
  )
  save(result, file = out_file)
  cat(sprintf("[FULL %dPCs] result saved → %s\n", pcs, out_file))
}

# ---------------------------------------------------------------------------
# 5. Launch both PCs in parallel --------------------------------------------
future_walk(full_postproc_files, run_full, .options = opts)

cat("========== All FULL runs completed at", format(Sys.time()), "==========\n")
cat("Total runtime:",
    round(difftime(Sys.time(), t0, units = "mins"), 2), "minutes\n")