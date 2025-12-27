# ==============================================================================
# SPATIAL TRANSCRIPTOMICS - NMI CALCULATION FOR DP-RST
# ==============================================================================

library(aricode)  # For NMI
library(dplyr)
library(DP.RST)  # For partition function

# Set Root Directory
ROOT_DIR <- "~/Desktop/DP-RST-workload/Datasets Results"
setwd(ROOT_DIR)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

get_ground_truth <- function(dataset, folder_name) {
  gt_path <- ""
  df <- NULL

  # For Gut/Brain/Breast/Prostate, try RData files first (they have correct coordinates for methods)
  if (dataset %in% c("Gut", "Brain", "Breast", "Prostate")) {
    rdata_files <- c(
      "Gut" = "Gut/swiss_roll_wt_muscle_finaltouches1.RData",
      "Brain" = "Brain/DFPLC_151510_data_for_iIMPACT.RData",
      "Breast" = "Breast/10x_breast_cancer_iIMPACT_data.RData",
      "Prostate" = "Prostate/10x_prostate_cancer_iIMPACT_data.RData"
    )

    rdata_path <- file.path(ROOT_DIR, rdata_files[[dataset]])
    if (file.exists(rdata_path)) {
      tryCatch({
        env <- new.env()
        load(rdata_path, envir = env)

        # Handle Gut's special case
        if (dataset == "Gut" && exists("swiss_roll_wt_muscle_finaltouches1", envir = env)) {
          df_raw <- env$swiss_roll_wt_muscle_finaltouches1
          if ("z" %in% colnames(df_raw) && "x" %in% colnames(df_raw) && "y" %in% colnames(df_raw)) {
            df <- data.frame(x = df_raw$x, y = df_raw$y, true_label = df_raw$z)
            df <- df[!is.na(df$true_label), ]
            message(paste("  Ground truth loaded from RData:", nrow(df), "observations"))
            return(df)
          }
        }

        # Handle Brain/Breast/Prostate with ground_truth_df
        if (exists("ground_truth_df", envir = env)) {
          df <- env$ground_truth_df
          # Ensure columns are x, y, type
          if ("type" %in% colnames(df) && "x" %in% colnames(df) && "y" %in% colnames(df)) {
            # Keep ALL rows including NAs - row indices correspond to DP-RST output
            df <- data.frame(x = df$x, y = df$y, true_label = df$type)
            message(paste("  Ground truth loaded from RData:", nrow(df), "observations total,",
                          sum(!is.na(df$true_label)), "with labels"))
            return(df)
          }
        }

        # Handle Breast/Prostate with z and loc variables (like Gut but keep NAs)
        if (dataset %in% c("Breast", "Prostate") && exists("z", envir = env) && exists("loc", envir = env)) {
          z_vec <- env$z
          loc_mat <- env$loc
          if (length(z_vec) == nrow(loc_mat)) {
            df <- data.frame(x = loc_mat[,1], y = loc_mat[,2], true_label = z_vec)
            message(paste("  Ground truth loaded from RData:", nrow(df), "observations total,",
                          sum(!is.na(df$true_label)), "with labels"))
            return(df)
          }
        }
      }, error = function(e) {
        message(paste("  Could not load RData, trying CSV:", e$message))
      })
    }
  }

  # Fallback to CSV files
  if (dataset == "Gut") {
    gt_path <- file.path(ROOT_DIR, "Gut", "gut_true_labels_with_barcodes.csv")
  } else if (dataset == "Brain") {
    gt_path <- file.path(ROOT_DIR, "Brain", "brain_true_labels_with_barcodes.csv")
  } else if (dataset == "Breast") {
    gt_path <- file.path(ROOT_DIR, "Breast", "breast_true_labels_with_barcodes.csv")
  } else if (dataset == "Prostate") {
    gt_path <- file.path(ROOT_DIR, "Prostate", "prostate_true_labels_with_barcodes.csv")
  } else if (dataset == "STARmap") {
    gt_path <- file.path(ROOT_DIR, "Starmap_Mouse", "starmap_pcs_coords_labels.csv")
  } else if (dataset == "osmFISH") {
    gt_path <- file.path(ROOT_DIR, "osmFISH", "osmfish_pcs_coords_labels.csv")
  } else if (dataset == "Lung") {
    gt_path <- file.path(ROOT_DIR, "Lung_Xenium", "Data_VUILD96MF", "cells_annotated.csv")
  } else if (dataset == "Colon") {
    seurat_path <- file.path(ROOT_DIR, "Visium_HD_Human_Colon_Cancer", "Nuclei_Data", "colon_hd_nuclei_seurat_object_processed.rds")
    if (file.exists(seurat_path)) {
      tryCatch({
        seurat_obj <- readRDS(seurat_path)
        meta_data <- seurat_obj@meta.data
        if ("annotation" %in% colnames(meta_data) && "x" %in% colnames(meta_data) && "y" %in% colnames(meta_data)) {
          df <- data.frame(
            x = meta_data$x,
            y = meta_data$y,
            true_label = meta_data$annotation
          )
          df <- df[!is.na(df$true_label), ]
          message(paste("  Ground truth loaded:", nrow(df), "observations with labels"))
          return(df)
        }
      }, error = function(e) {
        warning(paste("Error loading Colon:", e$message))
        return(NULL)
      })
    }
    return(NULL)
  }

  if (gt_path == "" || !file.exists(gt_path)) {
    warning(paste("Ground truth not found for", dataset))
    return(NULL)
  }

  raw <- read.csv(gt_path, stringsAsFactors = FALSE)

  # Find x, y, label columns
  x_col <- NULL
  y_col <- NULL
  label_col <- NULL

  # Check common column names
  for (col in c("x", "imagerow", "row", "image_row", "x_centroid")) {
    if (col %in% colnames(raw)) {
      x_col <- col
      break
    }
  }

  for (col in c("y", "imagecol", "col", "image_col", "y_centroid")) {
    if (col %in% colnames(raw)) {
      y_col <- col
      break
    }
  }

  for (col in c("label", "z", "true_label", "ground_truth", "type", "Annotation_Type", "annotation")) {
    if (col %in% colnames(raw)) {
      label_col <- col
      break
    }
  }

  if (is.null(x_col) || is.null(y_col) || is.null(label_col)) {
    warning(paste("Failed to find columns for", dataset))
    return(NULL)
  }

  df <- data.frame(
    x = raw[[x_col]],
    y = raw[[y_col]],
    true_label = raw[[label_col]]
  )
  # Keep all rows including NAs - row indices correspond to DP-RST output
  message(paste("  Ground truth loaded:", nrow(df), "observations total,",
                sum(!is.na(df$true_label)), "with labels"))
  return(df)
}

# ==============================================================================
# MAIN SCRIPT
# ==============================================================================

datasets_map <- list(
  "Gut" = "Gut",
  "Brain" = "Brain",
  "Breast" = "Breast",
  "Prostate" = "Prostate",
  "STARmap" = "Starmap_Mouse",
  "Lung" = "Lung_Xenium",
  "Colon" = "Visium_HD_Human_Colon_Cancer",
  "osmFISH" = "osmFISH"
)

results_pc3 <- matrix(NA, nrow = length(datasets_map), ncol = 1)
results_pc10 <- matrix(NA, nrow = length(datasets_map), ncol = 1)
rownames(results_pc3) <- names(datasets_map)
colnames(results_pc3) <- "DP-RST"
rownames(results_pc10) <- names(datasets_map)
colnames(results_pc10) <- "DP-RST"

for (d_name in names(datasets_map)) {
  folder_name <- datasets_map[[d_name]]
  message(paste("\n===== Processing Dataset:", d_name, "====="))

  tryCatch({
    gt <- get_ground_truth(d_name, folder_name)
    if (is.null(gt)) {
      message("  -> Skipping (No Ground Truth)")
      next
    }

    res_path <- file.path(ROOT_DIR, folder_name, "Results")

    for (pc in c(3, 10)) {
      target_mat <- if (pc == 3) results_pc3 else results_pc10
      message(paste("\n  --- Processing PC =", pc, "---"))

      # Special handling for Lung and Colon - use ConsensusMCM results
      if (d_name %in% c("Lung", "Colon")) {
        consensus_base <- "~/Desktop/DP-RST-workload/ConsensusMCM/ReferenceBasedAlignment"
        dataset_folder_map <- c(
          "Lung" = "results_Lung_Xenium",
          "Colon" = "results_Visium_HD_Human_Colon_Cancer"
        )
        consensus_path <- file.path(consensus_base, paste0(dataset_folder_map[[d_name]], "_PC", pc))
        consensus_file <- file.path(consensus_path, "final_assignments.csv")

        if (!file.exists(consensus_file)) {
          # Try to unzip
          zip_file <- paste0(consensus_path, ".zip")
          if (file.exists(zip_file)) {
            message(paste("    Unzipping:", basename(zip_file)))
            system(paste("unzip -q -o", shQuote(zip_file), "-d", shQuote(consensus_base)))
          }
        }

        if (file.exists(consensus_file)) {
          message(paste("    Loading DP-RST from ConsensusMCM:", basename(consensus_file)))
          consensus_df <- read.csv(consensus_file, stringsAsFactors = FALSE)

          if ("global_metacluster" %in% colnames(consensus_df) && "ann_shard" %in% colnames(consensus_df)) {
            # Filter to valid rows (non-NA in both ground truth and clustering)
            valid_rows <- !is.na(consensus_df$ann_shard) & !is.na(consensus_df$global_metacluster)
            if (sum(valid_rows) > 0) {
              nmi_val <- NMI(consensus_df$ann_shard[valid_rows], consensus_df$global_metacluster[valid_rows])
              target_mat[d_name, "DP-RST"] <- nmi_val
              message(paste("    DP-RST NMI:", round(nmi_val, 4), "- Clusters:", length(unique(consensus_df$global_metacluster[valid_rows])),
                            "(", sum(valid_rows), "of", nrow(consensus_df), "obs with labels)"))
            } else {
              message("    DP-RST: No valid observations with both labels and clustering")
            }
          } else {
            message("    DP-RST: Missing required columns in ConsensusMCM file")
          }
        } else {
          message("    DP-RST: ConsensusMCM file not found")
        }

        if (pc == 3) results_pc3 <- target_mat else results_pc10 <- target_mat
        next
      }

      # Standard DP-RST file patterns for other datasets
      dprst_files <- c()

      # Dataset-specific file paths
      if (d_name == "Brain") {
        dprst_file <- file.path(res_path, paste0("Brain_DP.RST_FromNewBastPT_newBS_p", pc, "_Version2_OutputOnly.rds"))
        if (file.exists(dprst_file)) dprst_files <- c(dprst_file)
      } else if (d_name == "Breast") {
        dprst_file <- file.path(res_path, paste0("Breast_DP.RST_newBS_p", pc, "_Version2_OutputOnly.RData"))
        if (file.exists(dprst_file)) dprst_files <- c(dprst_file)
      } else if (d_name == "Prostate") {
        dprst_file <- file.path(res_path, paste0("Prostate_DP.RST_newBS_p", pc, "_Version2_OutputOnly.rds"))
        if (file.exists(dprst_file)) dprst_files <- c(dprst_file)
      } else if (d_name == "STARmap") {
        dprst_file <- file.path(res_path, "DP-RST-files", paste0("Full_Result_DP.RST_", pc, "PCs.RData"))
        if (file.exists(dprst_file)) dprst_files <- c(dprst_file)
      } else if (d_name == "osmFISH") {
        dprst_file <- file.path(res_path, paste0("Full_Result_DP.RST_", pc, "PCs.RData"))
        if (file.exists(dprst_file)) dprst_files <- c(dprst_file)
      } else {
        # Generic search for Gut and others
        dprst_pattern <- paste0("*DP.RST.*p", pc, ".*OutputOnly\\.RData")
        dprst_files <- list.files(res_path, pattern = dprst_pattern, full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
      }

      if (length(dprst_files) == 0) {
        message("    DP-RST: File not found")
        next
      }

      message(paste("    Loading DP-RST from:", basename(dprst_files[1])))

      tryCatch({
        # Load the DP-RST output - handle both .rds and .RData files
        output_obj <- NULL

        if (grepl("\\.rds$", dprst_files[1], ignore.case = TRUE)) {
          # Load .rds file directly
          output_obj <- readRDS(dprst_files[1])
          message(paste("      Loaded .rds file directly"))
        } else {
          # Load .RData file
          env <- new.env()
          load(dprst_files[1], envir = env)

          # Find the output object - get the first object that is a list
          for (obj_name in ls(envir = env)) {
            obj <- get(obj_name, envir = env)
            if (is.list(obj) && ("cluster_out" %in% names(obj) || "teams_out" %in% names(obj) || "mcmc.draws.alpha_kl" %in% names(obj))) {
              output_obj <- obj
              message(paste("      Found output object:", obj_name))
              break
            }
          }
        }

        if (is.null(output_obj)) {
          message("    DP-RST: Could not find valid output object")
          next
        }

        # Extract mode-based partition as specified by user
        mode_based_partition <- partition(DP.RST_output = output_obj, method = "mode_based", batch_size = 100)

        # Extract teams partition
        X <- table(sequence(length(mode_based_partition$groups_partition)),
                   mode_based_partition$groups_partition)
        Z <- table(sequence(length(mode_based_partition$teams_partition)),
                   mode_based_partition$teams_partition)
        obs_in_teams <- X %*% Z
        obs_in_teams_vec_mode_based <- obs_in_teams %*% sort(unique(mode_based_partition$teams_partition))

        # Convert to cluster labels (vector of team assignments)
        cluster_labels <- as.vector(obs_in_teams_vec_mode_based)

        # Match with ground truth - filter out NAs
        if (length(cluster_labels) != nrow(gt)) {
          message(paste("    DP-RST: Length mismatch - clusters:", length(cluster_labels), "vs gt:", nrow(gt)))
          next
        }

        # Filter to only observations with non-NA ground truth labels
        valid_idx <- !is.na(gt$true_label)
        if (sum(valid_idx) == 0) {
          message("    DP-RST: No valid ground truth labels")
          next
        }

        cluster_labels_valid <- cluster_labels[valid_idx]
        gt_labels_valid <- gt$true_label[valid_idx]

        # Calculate NMI
        nmi_val <- NMI(gt_labels_valid, cluster_labels_valid)
        target_mat[d_name, "DP-RST"] <- nmi_val
        message(paste("    DP-RST NMI:", round(nmi_val, 4), "- Clusters:", length(unique(cluster_labels_valid)),
                      "(", sum(valid_idx), "of", length(cluster_labels), "obs with labels)"))

      }, error = function(e) {
        message(paste("    DP-RST error:", e$message))
        traceback()
      })

      if (pc == 3) results_pc3 <- target_mat else results_pc10 <- target_mat
    }
  }, error = function(e) {
    message(paste("  !!! ERROR processing", d_name, ":", e$message))
  })
}

# ==============================================================================
# PRINT AND SAVE RESULTS
# ==============================================================================

message("\n\n===== RESULTS SUMMARY =====\n")
message("PC = 3 Results:")
print(round(results_pc3, 4))
message("\n\nPC = 10 Results:")
print(round(results_pc10, 4))

write.csv(results_pc3, file = file.path(ROOT_DIR, "NMI_DP-RST_PC3.csv"), row.names = TRUE)
write.csv(results_pc10, file = file.path(ROOT_DIR, "NMI_DP-RST_PC10.csv"), row.names = TRUE)

message("\n\n===== RESULTS SAVED =====")
message("Results saved to:")
message("  - NMI_DP-RST_PC3.csv")
message("  - NMI_DP-RST_PC10.csv")
