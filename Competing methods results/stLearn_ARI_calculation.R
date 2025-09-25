#!/usr/bin/env Rscript

# Required packages
library(dplyr)
library(readr)
library(mclust)
library(tibble)
library(stringr)

# --- Configuration ---

# Directory containing stLearn clustering result files (both with_image and no_image)
stlearn_dir <- "DP.RST/Competing_Methods/stLearn_runs/res_realdata_csv"

# File pattern: e.g., "Brain_clustering_10PCs_5clusters_no_image.csv" or "Brain_clustering_10PCs_5clusters_with_image.csv"
file_pattern <- "^(Brain|Breast|Gut|Gut_reduced|Prostate)_clustering_(\\d+)PCs_(\\d+)clusters_(with|no)_image\\.csv$"

# True-label configuration for each dataset.
# Adjust the true_labels_file and true_label_col as needed.
true_label_config <- list(
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
  gut_reduced = list(
    true_labels_file = "/DP.RST/Gut/gut_true_labels.csv",
    true_label_col   = "z"
  ),
  prostate = list(
    true_labels_file = "/DP.RST/Prostate/prostate_true_labels.csv",
    true_label_col   = "z"
  )
)

# --- Function to compute ARI for one stLearn file ---
compute_stlearn_ari <- function(file, config) {
  # Read the stLearn clustering results file.
  st_df <- read_csv(file, show_col_types = FALSE)

  # Expect the file to have columns: spot, x, y, and a predicted cluster column (e.g., "cluster_5")
  pred_col <- grep("^cluster_\\d+$", names(st_df), value = TRUE)
  if (length(pred_col) != 1) {
    stop("File ", file, " must have exactly one prediction column matching '^cluster_\\d+$'.")
  }

  # Read the true labels.
  true_df <- read_csv(config$true_labels_file, show_col_types = FALSE)

  # Merge stLearn file with true labels.
  # In stLearn file, coordinates are (x, y) but in the true labels the coordinates are swapped.
  merged_df <- merge(st_df, true_df, by.x = c("x", "y"), by.y = c("y", "x"))

  # Filter out rows with missing true labels.
  true_label_col <- config$true_label_col
  merged_df <- merged_df[!is.na(merged_df[[true_label_col]]), ]

  # Compute ARI between predicted labels and true labels.
  ari_val <- adjustedRandIndex(as.factor(merged_df[[pred_col]]),
                               as.factor(merged_df[[true_label_col]]))
  return(ari_val)
}

# --- Main Execution ---

# List all stLearn files in the directory that match the pattern.
all_files <- list.files(stlearn_dir, pattern = file_pattern, full.names = TRUE, ignore.case = TRUE)
if (length(all_files) == 0) {
  stop("No stLearn files found matching the pattern in: ", stlearn_dir)
}

results_list <- list()

# Process each file.
for (f in all_files) {
  bn <- basename(f)
  # Use regex to extract:
  # Group 1: dataset name, Group 2: PCs, Group 3: clusters, Group 4: image indicator ("with" or "no")
  m <- str_match(bn, file_pattern)
  if (is.na(m[1,1])) next
  dataset <- tolower(m[1,2])
  PCs <- as.numeric(m[1,3])
  clusters <- as.numeric(m[1,4])
  image_indicator <- m[1,5]
  image_status <- ifelse(tolower(image_indicator) == "with", "with_image", "no_image")

  # Check that we have a true label config for this dataset.
  if (!(dataset %in% names(true_label_config))) {
    warning("No true label configuration for dataset: ", dataset, " (file: ", bn, ")")
    next
  }

  # Compute ARI.
  ari_val <- tryCatch({
    compute_stlearn_ari(f, true_label_config[[dataset]])
  }, error = function(e) {
    warning("Error processing file ", f, ": ", e$message)
    NA_real_
  })

  # Store the result.
  results_list[[length(results_list) + 1]] <- data.frame(
    dataset = dataset,
    PCs = PCs,
    clusters = clusters,
    image_status = image_status,
    ARI = ari_val,
    stringsAsFactors = FALSE
  )
}

# Combine all results into one data frame.
final_df <- do.call(rbind, results_list)
final_df <- final_df %>% arrange(dataset, PCs, image_status, clusters)

# Split into four tables:
results_pc3_with <- final_df %>% filter(PCs == 3, image_status == "with_image")
results_pc3_no <- final_df %>% filter(PCs == 3, image_status == "no_image")
results_pc10_with <- final_df %>% filter(PCs == 10, image_status == "with_image")
results_pc10_no <- final_df %>% filter(PCs == 10, image_status == "no_image")

# --- Output Results to a Text File ---
output_txt <- file.path("/DP.RST/Competing_Methods", "stLearn_ARI_results.txt")
sink(output_txt)
cat("\n--- Adjusted Rand Index (ARI) Results for PC = 3, with image ---\n\n")
results_pc3_with
cat("\n\n--- Adjusted Rand Index (ARI) Results for PC = 3, no image ---\n\n")
results_pc3_no
cat("\n\n--- Adjusted Rand Index (ARI) Results for PC = 10, with image ---\n\n")
results_pc10_with
cat("\n\n--- Adjusted Rand Index (ARI) Results for PC = 10, no image ---\n\n")
results_pc10_no
cat("\nSaved combined ARI results to:\n", output_txt, "\n")
sink()  # Reset output to console

cat("Script finished. Results saved to:", output_txt, "\n")
