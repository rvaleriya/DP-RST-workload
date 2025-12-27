# ---------------------------------------------------------------------------
# Prepare BAST init bundles for osmFISH (3PCs and 10PCs) from BASS init CSVs
# ---------------------------------------------------------------------------

# Library path 
my_lib_path <- "/scratch/user/varogovchenko/Rlibs"
.libPaths(my_lib_path)

library(DP.RST)
library(MASS)
library(igraph)
library(dplyr)
library(mclust)

set.seed(42)

# -------------------- Config --------------------
M        <- 10
base_in  <- "/scratch/user/varogovchenko/BASTION_HPRC/osmFISH"
save_f   <- file.path(base_in, "Outputs")

pc_source_csv <- file.path(base_in, "osmfish_pcs_coords_labels.csv")
pcs <- read.csv(pc_source_csv, check.names = FALSE)

# boundary (manual) → use rows with type == "original"
bnd_path <- file.path(base_in, "osmfish_manual_boundary.csv")
boundary_df <- read.csv(bnd_path, check.names = FALSE)

# helper code 
source("/scratch/user/varogovchenko/BASTION_HPRC/Extra_fun/MST_clusters_shredding.R")
source("/scratch/user/varogovchenko/BASTION_HPRC/Extra_fun/MST_from_assignments.R")
source("/scratch/user/varogovchenko/BASTION_HPRC/Extra_fun/merge_all_MST_components.R")

# -------------------- Helpers --------------------

# Build init bundle (labels + MST) for a given PC count using Full_Run CSVs
build_init_for_pcs <- function(pc_count) {
  
  # Input CSV from BASS init (must contain: x, y, label, PC_1..PC_pc_count)
  init_csv <- file.path(save_f, sprintf("InitBAST_from_BASS_%dPCs.csv", pc_count))
  init_df <- read.csv(init_csv, check.names = FALSE)
  
  # Data matrices
  coords_raw <- as.matrix(init_df[, c("x","y")])        # unscaled for boundary scaling
  coords_sc  <- scale(coords_raw)                        # scaled for MST geometry
  
  coords_jittered <- coords_sc
  coords_jittered[,1] <- jitter(coords_sc[,1], amount = 0.0001)
  coords_jittered[,2] <- jitter(coords_sc[,2], amount = 0.0001)
  
  Y          <- as.matrix(pcs[, sprintf("PC%d", 1:pc_count), drop = FALSE])
  Y_std      <- scale(Y)
  
  n <- nrow(Y_std); p <- ncol(Y_std)
  
  # Boundary scaled 
  bnd_scaled <- data.frame(
    x = (boundary_df$x - mean(coords_raw[,1])) / sd(coords_raw[,1]),
    y = (boundary_df$y - mean(coords_raw[,2])) / sd(coords_raw[,2])
  )
  
  # Prepare the frame for compute_MST_spatial()
  work_df <- init_df
  work_df$x_sc <- coords_jittered[,1]
  work_df$y_sc <- coords_jittered[,2]
  
  # --- Build MST & merge tiny spatial components along MST-only  ---
  init_partition <- compute_MST_spatial(
    data          = work_df,
    coords_cols   = c("x_sc", "y_sc"),
    cluster_col   = "label",
    bnd           = bnd_scaled,
    threshold     = 5000,
    penalty       = NULL,
    min_comp_size = 10
  )
  
  k_orig   <- length(unique(init_partition$data$spatial_cluster))
  k_merged <- length(unique(init_partition$spatial_cluster_merged))
  cat(sprintf("[Full_Run %2dPCs] Spatial comps before: %d | after MST-merge: %d\n",
              pc_count, k_orig, k_merged))
  
  clust_s   <- as.integer(init_partition$spatial_cluster_merged)
  graph0    <- init_partition$graph_base
  mstgraph0 <- init_partition$mst
  
  # Ensure base graph is connected (some meshes may still be multi-component)
  comp_no <- igraph::components(graph0)$no
  cat(sprintf("[Full_Run %2dPCs] graph0 components: %d\n", pc_count, comp_no))
  if (comp_no > 1L) {
    graph0 <- merge_all_components_to_giant(graph0, coords_jittered, verbose = FALSE)
    cat(sprintf("[Full_Run %2dPCs] Components after connect: %d\n",
                pc_count, igraph::components(graph0)$no))
  }
  
  # Clean MST edge attrs that are specific to the penalized build
  if ("raw_len" %in% igraph::edge_attr_names(mstgraph0))
    mstgraph0 <- igraph::delete_edge_attr(mstgraph0, "raw_len")
  if ("pen_w" %in% igraph::edge_attr_names(mstgraph0))
    mstgraph0 <- igraph::delete_edge_attr(mstgraph0, "pen_w")
  
  # --- Build BAST init lists (DP-RST expects mstgraph_lst; no 'trees') ---
  k_max_s <- length(unique(clust_s))
  
  # cluster means in standardized PC space
  cluster_means_matrix <- matrix(0, nrow = k_max_s, ncol = p)
  for (i in 1:k_max_s) {
    cluster_means_matrix[i, ] <- colMeans(Y_std[clust_s == i, , drop = FALSE])
  }
  
  cluster_mat <- matrix(rep(clust_s, M), nrow = n, ncol = M, byrow = FALSE)
  mu_list     <- vector("list", M)
  sigmasq_y   <- vector("list", M)
  mst_list    <- vector("list", M)
  covY        <- cov(Y_std)
  
  for (m in seq_len(M)) {
    mu_list[[m]]    <- cluster_means_matrix
    sigmasq_y[[m]]  <- covY
    mst_list[[m]]   <- mstgraph0
  }
  
  init_val <- list(
    trees = mst_list,     
    mu           = mu_list,
    cluster      = cluster_mat,
    sigmasq_y    = sigmasq_y
  )
  
  hyperpar <- list(
    sigmasq_mu = (0.5/(2*sqrt(1)))^2,
    lambda_s   = diag(1, p),
    nu         = p,
    M          = M,
    k_max      = k_max_s,
    lambda_k   = k_max_s 
  )
  
  list(
    Y_std   = Y_std,
    coords  = coords_jittered,
    graph0  = graph0,
    init_val = init_val,
    hyperpar = hyperpar,
    clust_s = clust_s
  )
}

# -------------------- Build & Save (3PCs and 10PCs) --------------------

for (pc_count in c(3L, 10L)) {
  out <- build_init_for_pcs(pc_count)
  
  out_file <- file.path(
    save_f,
    sprintf("Preproc_BAST_%dPCs.RData", pc_count)
  )
  
  # Save exactly the objects DP-RST steps will need
  Y_std   <- out$Y_std
  coords  <- out$coords
  graph0  <- out$graph0
  init_val <- out$init_val
  hyperpar <- out$hyperpar
  clust_s  <- out$clust_s
  
  save(list = c("Y_std","coords","graph0","init_val","hyperpar","clust_s"),
       file = out_file)
  cat(sprintf("✔ Saved: %s\n", out_file))
}
cat("Done preparing osmFISH for BAST")