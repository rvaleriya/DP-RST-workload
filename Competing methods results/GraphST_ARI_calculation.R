# Required packages
library(dplyr)
library(readr)
library(mclust)
library(tibble)
library(stringr)

# --- Configuration ---

# Base paths - Adjust if necessary
results_dir <- "/DP.RST/Competing_Methods/GraphST_runs/res_realdata_csv"
true_labels_base_dir <- "/DP.RST"

datasets <- list(
  Brain = list(
    results_file = file.path(results_dir, "graphst_brain_results.csv"),
    true_labels_file = file.path(true_labels_base_dir, "Brain", "brain_true_labels.csv"),
    true_label_col = "type"
  ),
  Breast = list(
    results_file = file.path(results_dir, "graphst_breast_results.csv"),
    true_labels_file = file.path(true_labels_base_dir, "Breast", "breast_true_labels.csv"),
    true_label_col = "z"
  ),
  Gut = list(
    results_file = file.path(results_dir, "graphst_gut_reduced_results.csv"),
    true_labels_file = file.path(true_labels_base_dir, "Gut", "gut_true_labels.csv"),
    true_label_col = "z"
  ),
  Prostate = list(
    results_file = file.path(results_dir, "graphst_prostate_results.csv"),
    true_labels_file = file.path(true_labels_base_dir, "Prostate", "prostate_true_labels.csv"),
    true_label_col = "z"
  )
)

# Pattern to identify clustering result columns
cluster_col_pattern <- "^(mclust_|kmeans_k=)\\d+$"

# --- Main Processing Loop ---

all_results_simple <- list()

for (dataset_name in names(datasets)) {
  message(paste("\nProcessing:", dataset_name))
  config <- datasets[[dataset_name]]
  true_label_col_name <- config$true_label_col

  # Read data
  results_df <- read_csv(config$results_file, show_col_types = FALSE, lazy = FALSE)
  true_labels_df <- read_csv(config$true_labels_file, show_col_types = FALSE, lazy = FALSE)

  # Merge using x and y columns
  merged_data <- merge(results_df, true_labels_df, by = c("x", "y"))
  # If merge doesn't match the expected number of rows, swap coordinates in true_labels_df
  if(nrow(merged_data) != nrow(true_labels_df)) {
    true_labels_df <- true_labels_df[, c("y", "x", true_label_col_name)]
    names(true_labels_df) <- c("x", "y", true_label_col_name)
    merged_data <- merge(results_df, true_labels_df, by = c("x", "y"))
  }

  # Filter out rows with NA true labels
  merged_data_subset <- merged_data[ !is.na(merged_data[[true_label_col_name]]), ]

  # Identify clustering result columns
  cluster_cols <- grep(cluster_col_pattern, names(merged_data_subset), value = TRUE)

  # Calculate ARI for each clustering method column
  ari_results <- list()
  true_labels <- merged_data_subset[[true_label_col_name]]

  for (col in cluster_cols) {
    cluster_labels <- merged_data_subset[[col]]
    ari_value <- tryCatch({
      adjustedRandIndex(as.factor(true_labels), as.factor(cluster_labels))
    }, error = function(e) {
      warning(paste("ARI calculation error for", col, "in", dataset_name, ":", e$message))
      NA_real_
    })

    ari_results[[col]] <- ari_value
  }

  all_results_simple[[dataset_name]] <- ari_results
}

# --- Output Results ---

# Save results to the specified folder
output_file <- "/DP.RST/Competing_Methods/GraphST_ARI_results.txt"
sink(output_file)
cat("\n\n--- Adjusted Rand Index (ARI) Results ---\n")

for (dataset_name in names(all_results_simple)) {
  cat("\n-- Dataset:", dataset_name, "--\n")
  result_list <- all_results_simple[[dataset_name]]

  if (is.null(result_list) || length(result_list) == 0) {
    cat("No results calculated.\n")
    next
  }

  # Create a simple tibble for output
  result_table <- tibble(
    Cluster_Method = names(result_list),
    ARI = unlist(result_list)
  )

  # Sort numerically based on trailing numbers in the cluster method name
  k_numeric <- suppressWarnings(as.numeric(str_extract(result_table$Cluster_Method, "\\d+$")))
  if (!all(is.na(k_numeric))) {
    result_table <- result_table %>% arrange(k_numeric)
  }

  print(result_table, n = Inf)
}

cat("\n\nScript finished.\n")
sink()  # Reset output to the console
