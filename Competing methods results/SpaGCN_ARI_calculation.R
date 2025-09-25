#!/usr/bin/env Rscript

# Required packages
library(dplyr)
library(readr)
library(mclust)
library(stringr)

# --- Configuration ---

# Directory containing SpaGCN clustering results files
results_dir <- "DP.RST/Competing_Methods/SpaGCN_runs/res_realdata_csv"

# Pattern to match clustering results files
# e.g., "all_clustering_results_Brain.csv", "all_clustering_results_Breast.csv", etc.
file_pattern <- "^all_clustering_results_(Brain|Breast|Gut|Gut_reduced|Prostate)\\.csv$"

# True label configuration for each dataset.
# (Assuming you already created files with barcodes, e.g., "brain_true_labels_with_barcodes.csv")
true_label_config <- list(
  brain = list(
    true_labels_file = "/DP.RST/Brain/brain_true_labels_with_barcodes.csv",
    true_label_col   = "type"
  ),
  breast = list(
    true_labels_file = "/DP.RST/Breast/breast_true_labels_with_barcodes.csv",
    true_label_col   = "z"
  ),
  gut = list(
    true_labels_file = "/DP.RST/Gut/gut_true_labels_with_barcodes.csv",
    true_label_col   = "z"
  ),
  gut_reduced = list(
    true_labels_file = "/DP.RST/Gut/gut_true_labels_with_barcodes.csv",
    true_label_col   = "z"
  ),
  prostate = list(
    true_labels_file = "/DP.RST/Prostate/prostate_true_labels_with_barcodes.csv",
    true_label_col   = "z"
  )
)

# --- Function to process one SpaGCN results file ---
compute_ari_file <- function(file, config) {
  # Read the clustering results file.
  df <- read_csv(file, show_col_types = FALSE)

  # The first column is the barcode but its header might be missing.
  if (!"barcode" %in% names(df)) {
    names(df)[1] <- "barcode"
  }

  # Read the corresponding true labels file.
  true_df <- read_csv(config$true_labels_file, show_col_types = FALSE)

  # Merge on the barcode.
  merged_df <- merge(df, true_df, by = "barcode")

  # Identify the prediction columns; these start with "refined_pred_"
  pred_cols <- grep("^refined_pred_", names(merged_df), value = TRUE)

  results <- list()

  # For each prediction column, compute ARI and extract parameters.
  # Expected column name format:
  # refined_pred_<clusters>clusters_<PCs>pcs_(with_histology|without_histology)
  pattern <- "^refined_pred_(\\d+)clusters_(\\d+)pcs_(with_histology|without_histology)$"

  for (col in pred_cols) {
    m <- str_match(col, pattern)
    if (is.na(m[1,1])) next  # Skip if the column name doesn't match the expected format.

    cluster_num <- as.numeric(m[1,2])
    PCs <- as.numeric(m[1,3])
    modality <- m[1,4]

    # Compute ARI comparing the predicted labels to the true labels.
    ari_val <- tryCatch({
      adjustedRandIndex(as.factor(merged_df[[col]]),
                        as.factor(merged_df[[config$true_label_col]]))
    }, error = function(e) {
      warning(paste("Error computing ARI in file", file, "for column", col, ":", e$message))
      NA_real_
    })

    results[[length(results) + 1]] <- data.frame(
      dataset = tolower(gsub("all_clustering_results_(.*)\\.csv", "\\1", basename(file))),
      PCs = PCs,
      clusters = cluster_num,
      modality = modality,
      ARI = ari_val,
      stringsAsFactors = FALSE
    )
  }

  if (length(results) > 0) {
    return(do.call(rbind, results))
  } else {
    return(NULL)
  }
}

# --- Main Execution ---

# List all SpaGCN clustering results files
files <- list.files(results_dir, pattern = file_pattern, full.names = TRUE, ignore.case = TRUE)
if (length(files) == 0) {
  stop("No clustering results files found in: ", results_dir)
}

all_results <- list()
for (f in files) {
  # Extract dataset name from filename using the pattern.
  m <- str_match(basename(f), file_pattern)
  if (is.na(m[1,1])) next
  dataset <- tolower(m[1,2])

  if (!(dataset %in% names(true_label_config))) {
    warning("No true label configuration for dataset: ", dataset)
    next
  }

  res <- compute_ari_file(f, true_label_config[[dataset]])
  if (!is.null(res)) {
    all_results[[length(all_results) + 1]] <- res
  }
}

# Combine all results into a single data frame.
final_df <- do.call(rbind, all_results)
final_df <- final_df %>% arrange(dataset, PCs, modality, clusters)

# Split the results into four tables:
# PC = 3 with_histology, PC = 3 without_histology, PC = 10 with_histology, PC = 10 without_histology.
results_pc3_with <- final_df %>% filter(PCs == 3, modality == "with_histology")
results_pc3_without <- final_df %>% filter(PCs == 3, modality == "without_histology")
results_pc10_with <- final_df %>% filter(PCs == 10, modality == "with_histology")
results_pc10_without <- final_df %>% filter(PCs == 10, modality == "without_histology")

# --- Output Results to a Text File ---
output_txt <- file.path("/DP.RST/Competing_Methods", "SpaGCN_ARI_results.txt")
sink(output_txt)
cat("\n--- Adjusted Rand Index (ARI) Results for PC = 3, with_histology ---\n\n")
results_pc3_with
cat("\n\n--- Adjusted Rand Index (ARI) Results for PC = 3, without_histology ---\n\n")
results_pc3_without
cat("\n\n--- Adjusted Rand Index (ARI) Results for PC = 10, with_histology ---\n\n")
results_pc10_with
cat("\n\n--- Adjusted Rand Index (ARI) Results for PC = 10, without_histology ---\n\n")
results_pc10_without
cat("\nSaved combined ARI results to:\n", output_txt, "\n")
sink()  # Reset output to console

cat("Script finished. Results saved to:", output_txt, "\n")
