#!/usr/bin/env Rscript

# Required packages
library(dplyr)
library(readr)
library(mclust)
library(tibble)
library(stringr)

# --- Configuration ---

# Directory containing the SEDR result files
sedr_dir <- "/DP.RST/Competing_Methods/SEDR_runs/res_realdata_csv"

# Pattern to match SEDR files (e.g., "sedr_gut_reduced_pca10_results.csv")
sedr_file_pattern <- "^sedr_.*_pca\\d+_results\\.csv$"

# True-label configuration for each dataset
dataset_config <- list(
  brain = list(
    true_labels_file = "/DP.RST/Brain/brain_true_labels.csv",
    true_label_col   = "type"
  ),
  breast = list(
    true_labels_file = "/DP.RST/Breast/breast_true_labels.csv",
    true_label_col   = "z"
  ),
  gut = list(
    true_labels_file = "/DP.RST/Gut/gut_true_labels.csv",
    true_label_col   = "z"
  ),
  `gut_reduced` = list(
    true_labels_file = "/DP.RST/Gut/gut_true_labels.csv",
    true_label_col   = "z"
  ),
  prostate = list(
    true_labels_file = "/DP.RST/Prostate/prostate_true_labels.csv",
    true_label_col   = "z"
  )
)

# --- Function to compute ARI for one file ---
compute_sedr_ari <- function(csv_file, ds_config) {
  # Read the SEDR results
  results_df <- read_csv(csv_file, show_col_types = FALSE)

  # Read the true labels
  true_labels_df <- read_csv(ds_config$true_labels_file, show_col_types = FALSE)

  # Merge using x and y coordinates
  merged_data <- merge(results_df, true_labels_df, by = c("x", "y"))

  # If the merge doesn't have the expected number of rows, try swapping coordinates
  if (nrow(merged_data) != nrow(true_labels_df)) {
    true_labels_df <- true_labels_df[, c("y", "x", ds_config$true_label_col)]
    names(true_labels_df) <- c("x", "y", ds_config$true_label_col)
    merged_data <- merge(results_df, true_labels_df, by = c("x", "y"))
  }

  # Filter out rows with missing true labels
  lbl_col <- ds_config$true_label_col
  merged_data <- merged_data[!is.na(merged_data[[lbl_col]]), ]

  # Identify clustering result columns.
  # This regex matches column names like "mclust_pca3_3", "mclust_pca3_5", etc.
  cluster_cols <- grep("^(newmethod_|mclust_|kmeans_k=).+\\d+$", names(merged_data), value = TRUE)

  ari_list <- list()
  true_labels <- merged_data[[lbl_col]]

  for (col in cluster_cols) {
    cluster_labels <- merged_data[[col]]
    ari_val <- tryCatch({
      adjustedRandIndex(as.factor(true_labels), as.factor(cluster_labels))
    }, error = function(e) {
      warning(paste("ARI error in", csv_file, "for", col, ":", e$message))
      NA_real_
    })
    ari_list[[col]] <- ari_val
  }

  return(ari_list)
}

# --- Main Execution ---

# List all matching SEDR result files
new_method_files <- list.files(sedr_dir, pattern = sedr_file_pattern, full.names = TRUE)
if (length(new_method_files) == 0) {
  stop("No files found matching the pattern in: ", sedr_dir)
}

all_results <- list()

# Regex to extract dataset name and PCA number.
# This pattern captures dataset names (even with underscores, e.g., "gut_reduced")
# and the PCA number.
regex_pattern <- "^sedr_([^_]+(?:_[^_]+)*)_pca(\\d+)_results\\.csv$"

for (f in new_method_files) {
  bn <- basename(f)
  match <- regexec(regex_pattern, bn, ignore.case = TRUE)
  captures <- regmatches(bn, match)

  if (length(captures[[1]]) < 3) {
    warning("Filename does not match expected pattern: ", bn)
    next
  }

  ds_name <- tolower(captures[[1]][2])
  pc_number <- as.numeric(captures[[1]][3])

  if (!ds_name %in% names(dataset_config)) {
    warning("No dataset config found for: ", ds_name, " (file: ", bn, ")")
    next
  }

  ari_list <- compute_sedr_ari(f, dataset_config[[ds_name]])

  # Store results without the file name
  all_results[[bn]] <- list(
    dataset = ds_name,
    PCs = pc_number,
    ari_values = ari_list
  )
}

# Convert the list of results into a tidy data frame.
# Extract the cluster number from the clustering column name (trailing digits)
results_df <- do.call(rbind, lapply(all_results, function(entry) {
  if (length(entry$ari_values) == 0) {
    return(NULL)
  }
  data.frame(
    dataset = entry$dataset,
    PCs = entry$PCs,
    cluster = as.numeric(str_extract(names(entry$ari_values), "\\d+$")),
    ARI = unlist(entry$ari_values),
    stringsAsFactors = FALSE
  )
}))

if (is.null(results_df)) {
  stop("No ARI results were computed.")
}

# Remove any row names so file names are not printed
rownames(results_df) <- NULL

results_df <- results_df %>% arrange(dataset, PCs, cluster)

# Create two tables: one for PC = 3 and one for PC = 10
results_df_pc3 <- results_df %>% filter(PCs == 3) %>% select(dataset, PCs, cluster, ARI)
results_df_pc10 <- results_df %>% filter(PCs == 10) %>% select(dataset, PCs, cluster, ARI)

# --- Output Results to a Text File in the Desired Folder ---
output_txt <- file.path("/DP.RST/Competing_Methods", "SEDR_ARI_results.txt")
sink(output_txt)
cat("\n--- Adjusted Rand Index (ARI) Results for PC = 3 ---\n\n")
results_df_pc3
cat("\n\n--- Adjusted Rand Index (ARI) Results for PC = 10 ---\n\n")
results_df_pc10
cat("\nSaved combined ARI results to:\n", output_txt, "\n")
sink()  # Reset output to console

cat("Script finished. Results saved to:", output_txt, "\n")
