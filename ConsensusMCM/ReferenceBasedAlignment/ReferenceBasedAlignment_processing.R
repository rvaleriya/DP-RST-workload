#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(ggplot2)
library(mclust)

# Configuration
# Usage: Rscript ReferenceBasedAlignment_processing.R <dataset> <pc_version>
# Example: Rscript ReferenceBasedAlignment_processing.R Lung PC10
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  dataset <- "Colon"
  pc_version <- "PC3"
} else {
  dataset <- args[1]
  pc_version <- args[2]
}

base_dir <- "~/Desktop/DP-RST-workload/Shards_results"
ReferenceBasedAlignment_dir <- file.path(base_dir, "ReferenceBasedAlignment")

# Map dataset names to their ReferenceBasedAlignment result folders
dataset_map <- list(
  Lung = "results_Lung_Xenium",
  Colon = "results_Visium_HD_Human_Colon_Cancer"
)

ReferenceBasedAlignment_results <- file.path(ReferenceBasedAlignment_dir, paste0(dataset_map[[dataset]], "_", pc_version))
data_dir <- file.path(base_dir, dataset, pc_version)

cat("Processing:", dataset, pc_version, "\n")

# Load ReferenceBasedAlignment assignments
assignments_path <- file.path(ReferenceBasedAlignment_results, "local_to_global_hard_assignments.csv")
ReferenceBasedAlignment_assignments <- read_csv(assignments_path, show_col_types = FALSE)
cat("Loaded", nrow(ReferenceBasedAlignment_assignments), "ReferenceBasedAlignment assignments\n")

# Build local cluster to metacluster map
prop_dir <- file.path(data_dir, "proportions")
prop_files <- list.files(prop_dir, pattern = "^Shard_\\d+_proportions\\.csv$", full.names = TRUE)

cluster_map <- prop_files %>%
  sort() %>%
  map_dfr(function(f) {
    shard_id <- as.integer(str_extract(f, "(?<=Shard_)\\d+"))
    df <- read_csv(f, show_col_types = FALSE)
    df %>% transmute(shard_id = shard_id, local_cluster_id = cluster)
  }) %>%
  mutate(local_cluster_index = row_number() - 1L)

cat("Built map for", nrow(cluster_map), "clusters\n")

# Join with ReferenceBasedAlignment results
full_map <- inner_join(cluster_map, ReferenceBasedAlignment_assignments, by = "local_cluster_index")

# Load shard data
lists_dir <- file.path(data_dir, "lists")
rds_files <- list.files(lists_dir, pattern = "^Shard_\\d+_list\\.rds$", full.names = TRUE)

all_shards_df <- rds_files %>%
  sort() %>%
  map_dfr(function(f) {
    shard_id <- as.integer(str_extract(f, "(?<=Shard_)\\d+"))
    shard <- readRDS(f)
    
    loc <- as.matrix(shard$loc_shard)
    ann <- if ("ann_shard" %in% names(shard)) shard$ann_shard else NA
    
    tibble(
      shard_id = shard_id,
      x = loc[, 1],
      y = loc[, 2],
      shard_cluster = as.integer(shard$final_clusters),
      ann_shard = ann,
      idx_global = shard$idx_global
    )
  })

cat("Loaded", length(rds_files), "shards,", nrow(all_shards_df), "cells\n")

# Merge metacluster assignments
final_df <- all_shards_df %>%
  left_join(
    full_map %>% transmute(
      shard_id,
      shard_cluster = as.integer(local_cluster_id),
      global_metacluster = global_cluster_assignment
    ),
    by = c("shard_id", "shard_cluster")
  )

# Save results
output_csv <- file.path(ReferenceBasedAlignment_results, "final_assignments.csv")
write_csv(final_df, output_csv)
cat("Saved to:", output_csv, "\n")

# Compute ARI if annotations available
if (any(!is.na(final_df$ann_shard))) {
  valid <- final_df %>% filter(!is.na(ann_shard), !is.na(global_metacluster))
  if (nrow(valid) > 0) {
    ari <- adjustedRandIndex(valid$ann_shard, valid$global_metacluster)
    ari_subtext <- sprintf("ARI = %.4f", ari)
    cat("ARI:", round(ari, 4), "\n")
  }
}

# Generate spatial plot
p <- ggplot(final_df, aes(x = x, y = y, color = factor(global_metacluster))) +
  geom_point(size = 0.2, alpha = 0.8) +
  coord_fixed() +
  labs(title = paste(dataset, pc_version, "Metaclusters"),
       subtitle = ari_subtext,
       x = "X", y = "Y", color = "Metacluster") +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 3)))

plot_filename <- sprintf("spatial_plot_%s_%s.png", dataset, pc_version)
plot_path <- file.path(ReferenceBasedAlignment_results, plot_filename)
ggsave(plot_path, p, width = 10, height = 8, dpi = 300)
cat("Plot saved to:", plot_path, "\n")
