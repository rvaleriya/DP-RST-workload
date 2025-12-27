# ---------------------------------------------------------------------------
# Run fitBAST on osmFISH preproc files (3PCs & 10PCs) in parallel
# ---------------------------------------------------------------------------
t0 <- Sys.time()
cat("========== BAST osmFISH inference (3PCs & 10PCs, parallel) ==========\n")

# Library setup ------------------------------------------------------------
my_lib_path <- "/scratch/user/varogovchenko/Rlibs"
.libPaths(my_lib_path)

library(future)
library(furrr)

# BAST code
source("/scratch/user/varogovchenko/BASTION_HPRC/BASTFun_2.R")
source("/scratch/user/varogovchenko/BASTION_HPRC/ComplexDomainFun.R")

# Paths & MCMC controls ----------------------------------------------------
base_dir <- "/scratch/user/varogovchenko/BASTION_HPRC/osmFISH/Outputs"

MCMC   <- 30000
BURNIN <- 25000
THIN   <- 5

# Parallel backend ---------------------------------------------------------
slurm_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
num_workers <- if (!is.na(slurm_cores)) slurm_cores else future::availableCores()
plan(multicore, workers = min(2, num_workers))  # two jobs: 3PCs, 10PCs
opts <- furrr_options(seed = TRUE, scheduling = Inf)
cat("[INFO] Using", min(2, num_workers), "workers\n")

# Find preproc files (supports both naming styles) ------------------------
candidates <- c(
  file.path(base_dir, "Preproc_BAST_3PCs.RData"),
  file.path(base_dir, "Preproc_BAST_10PCs.RData")
)
preproc_files <- candidates[file.exists(candidates)]

get_pcs <- function(p) as.integer(sub("^.*_([0-9]+)PCs\\.RData$", "\\1", basename(p)))

# Worker: load preproc & run fitBAST --------------------------------------
run_full <- function(pre_file) {
  pcs <- get_pcs(pre_file)
  cat(sprintf("[BAST %2dPCs] Loading %s\n", pcs, pre_file))
  load(pre_file)
  # Expect: Y_std, graph0, init_val, hyperpar, coords, clust_s
  
  # PT ladder aligned to number of chains
  temp <- seq(0.1, 1.0, length.out = hyperpar$M)
  
  y_shard <- Y_std
  set.seed(100 + pcs)
  
  out_file <- file.path(base_dir, sprintf("Full_Result_BAST_%dPCs.RData", pcs))
  backup_d <- file.path(base_dir, sprintf("Full_Result_BAST_%dPCs_backup.RData", pcs))
  
  cat(sprintf("[BAST %2dPCs] n=%d p=%d k_max=%d M=%d | MCMC=%d, BURNIN=%d, THIN=%d\n",
              pcs, nrow(y_shard), ncol(y_shard), hyperpar$k_max, hyperpar$M,
              MCMC, BURNIN, THIN))
  
  result <- fitBAST(
    y_shard, graph0,
    init_val, hyperpar,
    MCMC    = MCMC,
    BURNIN  = BURNIN,
    THIN    = THIN,
    temp    = temp,
    PT      = TRUE,
    seed    = 100 + pcs,
    backup_d = backup_d
  )
  
  save(result, file = out_file)
  cat(sprintf("[BAST %2dPCs] Saved → %s (backup: %s)\n", pcs, out_file, backup_d))
}

# 5) Launch in parallel -------------------------------------------------------
future_walk(preproc_files, run_full, .options = opts)

cat("========== All BAST FULL runs completed at", format(Sys.time()), "==========\n")
cat("Total runtime:", round(difftime(Sys.time(), t0, units = "mins"), 2), "minutes\n")