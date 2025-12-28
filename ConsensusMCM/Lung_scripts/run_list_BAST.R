args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
list_file <- args[1]

## ---------- LOAD LIBS & YOUR FUNCTIONS -----------
my_lib_path <- "/scratch/user/varogovchenko/Rlibs"
.libPaths(my_lib_path)

library(rlang)
library(MASS)
library(igraph)
library(fields)
library(mclust)
library(dplyr)
library(future)
library(furrr)

source("/scratch/user/varogovchenko/BASTION_HPRC/BASTFun_2.R")
source("/scratch/user/varogovchenko/BASTION_HPRC/ComplexDomainFun.R")

## ---------- now plan() so workers inherit everything -----------
plan(multicore, workers = 28)

## ---------- read list & run shards ------------------------------
shards <- scan(list_file, what = "", quiet = TRUE)
source("/scratch/user/varogovchenko/BASTION_HPRC/Lung_Xenium/run_BAST_shard_fun.R")

start_time <- Sys.time()
cat(sprintf("\n[START] Processing %d shards at %s\n", length(shards), format(start_time)))
flush.console()

# Show what’s coming
cat(sprintf("[INFO] Workers: %d | Host: %s\n", nbrOfWorkers(), Sys.info()[["nodename"]]))
cat(sprintf("[INFO] List file: %s\n", list_file))
flush.console()

future_walk(
  shards,
  run_shard,
  .options = furrr_options(seed = TRUE, scheduling = Inf)
)

end_time <- Sys.time()
cat(sprintf("\n[END] Finished processing %d shards at %s\n", length(shards), format(end_time)))
cat(sprintf("[RUNTIME TOTAL] %s\n", format(difftime(end_time, start_time, units = "mins"))))
flush.console()

# ---------- OPTIONAL SUMMARY ----------
w <- warnings()
if (length(w)) {
  cat("\n---------- WARNING SUMMARY ----------\n")
  print(w)             # will print first 50, same as interactive
  cat("---------- END WARNING SUMMARY ------\n")
}