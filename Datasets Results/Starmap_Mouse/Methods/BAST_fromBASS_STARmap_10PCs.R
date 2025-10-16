# ---------------------------------------------------------------------------
# BAST: Single FULL dataset (sequential; no parallel)
# ---------------------------------------------------------------------------
t0 <- Sys.time()
cat("========== BAST FULL inference started at",
    format(Sys.time()), "==========\n")

# ---------------------------------------------------------------------------
# 0) Library setup -----------------------------------------------------------
my_lib_path <- "/scratch/user/varogovchenko/Rlibs"
.libPaths(my_lib_path)

library(rlang)
library(MASS)
library(igraph)
library(fields)
library(mclust)
library(dplyr)

# Your BAST code (sampler, helpers)
source("/scratch/user/varogovchenko/BASTION_HPRC/BASTFun_2.R")
source("/scratch/user/varogovchenko/BASTION_HPRC/ComplexDomainFun.R")

# ---------------------------------------------------------------------------
# 1) User parameters ---------------------------------------------------------
# Folder containing the preprocessed FULL data .RData you created earlier
full_dir      <- "/scratch/user/varogovchenko/BASTION_HPRC/STARmap"
preproc_file  <- file.path(full_dir, "Full_Preproc_BAST_10PCs.RData")  # <- 10PCs file

# Sampler controls
MCMC   <- 30000
BURNIN <- 25000
THIN   <- 5

# Parallel tempering ladder
temp <- seq(0.1, 1, by = 0.1)

# Output/backup paths
out_file   <- file.path(full_dir, "Full_Result_BAST_10PCs.RData")
backup_d   <- file.path(full_dir, "Full_Result_BAST_10PCs_backup.RData")

# ---------------------------------------------------------------------------
# 2) Load FULL preprocessed objects -----------------------------------------
stopifnot(file.exists(preproc_file))
cat("[INFO] Loading:", preproc_file, "\n")
load(preproc_file)
# Expecting: Y_std, coords_sc, graph0, init_val, hyperpar, clust_s, mstgraph0

# Alias for legacy fitBAST signature naming (if your code assumes y_shard)
y_shard <- Y_std

# Basic sanity prints
cat(sprintf("[INFO] n = %d, p = %d, k_max = %d, M = %d\n",
            nrow(Y_std), ncol(Y_std), hyperpar$k_max, hyperpar$M))

# ---------------------------------------------------------------------------
# 3) Run BAST on FULL data ---------------------------------------------------
cat("[INFO] Starting BAST on FULL data …\n")
set.seed(2025)
result <- fitBAST(
  y_shard, graph0,
  init_val, hyperpar,
  MCMC   = MCMC,
  BURNIN = BURNIN,
  THIN   = THIN,
  temp   = temp,
  PT     = TRUE,
  seed   = 7394,
  backup_d = backup_d
)

# ---------------------------------------------------------------------------
# 4) Save results ------------------------------------------------------------
save(result, file = out_file)
cat(sprintf("[INFO] Result saved → %s\n", out_file))
cat(sprintf("[INFO] Backup path  → %s\n", backup_d))

cat("========== BAST FULL inference completed at", format(Sys.time()), "==========\n")
cat("Total runtime:",
    round(difftime(Sys.time(), t0, units = "mins"), 2), "minutes\n")