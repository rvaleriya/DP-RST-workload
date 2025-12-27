# ==============================================================================
# SPATIAL TRANSCRIPTOMICS - NMI CALCULATION FOR COMPETING METHODS (FIXED)
# (Excluding DP-RST)
# ==============================================================================

library(aricode)  # For NMI
library(dplyr)
library(readr)
library(stringr)

# Set Root Directory
ROOT_DIR <- "~/Desktop/DP-RST-workload/Datasets Results"
setwd(ROOT_DIR)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

find_col_name <- function(df, candidates) {
  cols <- colnames(df)
  match <- intersect(candidates, cols)
  if (length(match) > 0) return(match[1])

  lower_cols <- tolower(cols)
  for (cand in candidates) {
    idx <- which(lower_cols == tolower(cand))
    if (length(idx) > 0) return(cols[idx[1]])
  }

  for (cand in candidates) {
    idx <- grep(cand, cols, ignore.case = TRUE)
    if (length(idx) > 0) return(cols[idx[1]])
  }

  return(NULL)
}

# Helper function to merge by barcode/ID instead of coordinates
merge_by_barcode <- function(gt, result, gt_barcode_col = "ID", res_barcode_col = "barcode") {
  if (gt_barcode_col %in% colnames(gt) && res_barcode_col %in% colnames(result)) {
    m <- merge(gt, result, by.x = gt_barcode_col, by.y = res_barcode_col)
    return(m)
  }
  return(NULL)
}

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
            df <- data.frame(x = df$x, y = df$y, true_label = df$type)
            df <- df[!is.na(df$true_label), ]
            message(paste("  Ground truth loaded from RData:", nrow(df), "observations"))
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

  raw <- if (grepl("\\.rds$", gt_path, ignore.case = T)) readRDS(gt_path) else read_csv(gt_path, show_col_types = FALSE)

  x_col <- find_col_name(raw, c("x", "imagerow", "row", "image_row", "x_centroid"))
  y_col <- find_col_name(raw, c("y", "imagecol", "col", "image_col", "y_centroid"))
  label_col <- find_col_name(raw, c("label", "z", "true_label", "ground_truth", "type", "Annotation_Type", "annotation"))
  id_col <- find_col_name(raw, c("ID", "barcode", "cell_id", "CellID"))

  if (is.null(x_col) || is.null(y_col) || is.null(label_col)) {
    warning(paste("Failed to find columns for", dataset))
    return(NULL)
  }

  # Include ID column if available (for osmFISH and others)
  if (!is.null(id_col)) {
    df <- data.frame(
      ID = raw[[id_col]],
      x = raw[[x_col]],
      y = raw[[y_col]],
      true_label = raw[[label_col]]
    )
  } else {
    df <- data.frame(
      x = raw[[x_col]],
      y = raw[[y_col]],
      true_label = raw[[label_col]]
    )
  }
  df <- df[!is.na(df$true_label), ]
  message(paste("  Ground truth loaded:", nrow(df), "observations with labels"))
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

method_names <- c("BayesSpace", "SC-MEB", "DR-SC", "SEDR", "GraphST",
                  "STAGATE", "SpaGCN (w/)", "SpaGCN (w/o)", "stLearn (w/)",
                  "stLearn (w/o)", "k-means", "BASS", "STAMarker", "SpaMask")

results_pc3 <- matrix(NA, nrow = length(datasets_map), ncol = length(method_names))
results_pc10 <- matrix(NA, nrow = length(datasets_map), ncol = length(method_names))
rownames(results_pc3) <- names(datasets_map); colnames(results_pc3) <- method_names
rownames(results_pc10) <- names(datasets_map); colnames(results_pc10) <- method_names

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

      # --- BayesSpace ---
      message("    Loading BayesSpace...")
      bs_file <- file.path(res_path, paste0("BayesSpace_", gsub("_Mouse", "", folder_name), "_", pc, "PCs.rds"))
      if (!file.exists(bs_file)) {
        bs_file <- file.path(res_path, paste0("BayesSpace_", gsub("Starmap_Mouse", "STARmap", gsub("Visium_HD_Human_Colon_Cancer", "Colon", gsub("Lung_Xenium", "Lung", folder_name))), "_", pc, "PCs.rds"))
      }
      if (!file.exists(bs_file)) {
        # Try generic name (for Gut/Brain/Breast/Prostate)
        bs_file <- file.path(res_path, paste0("BayesSpace_", folder_name, "_clustering_results.rds"))
      }

      if (file.exists(bs_file)) {
        tryCatch({
          res_df <- readRDS(bs_file)

          # Check if file already has truth column (Lung/Colon format)
          if ("truth" %in% colnames(res_df) && "label" %in% colnames(res_df)) {
            valid_rows <- !is.na(res_df$truth)
            if (sum(valid_rows) > 0) {
              nmi_val <- NMI(res_df$truth[valid_rows], res_df$label[valid_rows])
              target_mat[d_name, "BayesSpace"] <- nmi_val
              message(paste("      BayesSpace NMI:", round(nmi_val, 4), "using label column"))
            }
          } else {
            # Find the right cluster column
            k_cols <- grep(paste0("cluster.*pc", pc, "|cluster.*pc", pc), colnames(res_df), ignore.case=TRUE, value=TRUE)
            if (length(k_cols) == 0) k_cols <- grep("cluster|label", colnames(res_df), ignore.case=TRUE, value=TRUE)

            if (length(k_cols) > 0) {
              # Use the column with highest k value or first one
              cluster_col <- k_cols[length(k_cols)]

              # Get x,y columns
              if ("imagerow" %in% colnames(res_df)) {
                result <- data.frame(x = res_df$imagerow, y = res_df$imagecol, pred_label = res_df[[cluster_col]])
              } else {
                result <- data.frame(x = res_df$x, y = res_df$y, pred_label = res_df[[cluster_col]])
              }

              m <- merge(gt, result, by=c("x","y"))
              if (nrow(m) > 0) {
                nmi_val <- NMI(m$true_label, m$pred_label)
                target_mat[d_name, "BayesSpace"] <- nmi_val
                message(paste("      BayesSpace NMI:", round(nmi_val, 4), "using", cluster_col))
              } else {
                message("      BayesSpace: No matching observations")
              }
            }
          }
        }, error = function(e) message(paste("      BayesSpace error:", e$message)))
      } else {
        message("      BayesSpace: File not found")
      }

      # --- SC-MEB ---
      message("    Loading SC-MEB...")
      sc_file <- file.path(res_path, paste0("SC-MEB_", gsub("_Mouse", "", gsub("Starmap_Mouse", "STARmap", gsub("Lung_Xenium", "Lung", gsub("Visium_HD_Human_Colon_Cancer", "Colon", folder_name)))), "_", pc, "PCs.rds"))
      if (!file.exists(sc_file)) {
        sc_file <- file.path(res_path, paste0("SC-MEB_", folder_name, "_clustering_results.rds"))
      }

      if (file.exists(sc_file)) {
        tryCatch({
          res_df <- readRDS(sc_file)

          # Check if file already has truth column (Lung/Colon format)
          if ("truth" %in% colnames(res_df) && "label" %in% colnames(res_df)) {
            valid_rows <- !is.na(res_df$truth)
            if (sum(valid_rows) > 0) {
              nmi_val <- NMI(res_df$truth[valid_rows], res_df$label[valid_rows])
              target_mat[d_name, "SC-MEB"] <- nmi_val
              message(paste("      SC-MEB NMI:", round(nmi_val, 4), "using label column"))
            }
          } else {
            k_cols <- grep(paste0("pc", pc, "_k"), colnames(res_df), ignore.case=TRUE, value=TRUE)
            if (length(k_cols) == 0) k_cols <- grep("label|cluster", colnames(res_df), ignore.case=TRUE, value=TRUE)

            if (length(k_cols) > 0) {
              cluster_col <- k_cols[length(k_cols)]
              result <- data.frame(x = res_df$x, y = res_df$y, pred_label = res_df[[cluster_col]])

              m <- merge(gt, result, by=c("x","y"))
              if (nrow(m) > 0) {
                nmi_val <- NMI(m$true_label, m$pred_label)
                target_mat[d_name, "SC-MEB"] <- nmi_val
                message(paste("      SC-MEB NMI:", round(nmi_val, 4), "using", cluster_col))
              } else {
                message("      SC-MEB: No matching observations")
              }
            }
          }
        }, error = function(e) message(paste("      SC-MEB error:", e$message)))
      } else {
        message("      SC-MEB: File not found")
      }

      # --- DR-SC ---
      message("    Loading DR-SC...")
      dr_file <- file.path(res_path, paste0("DR-SC_", gsub("_Mouse", "", gsub("Starmap_Mouse", "STARmap", folder_name)), "_", pc, "PCs.rds"))
      if (!file.exists(dr_file)) {
        dr_file <- file.path(res_path, paste0("DRSC_", gsub("Visium_HD_Human_Colon_Cancer", "Colon", gsub("Lung_Xenium", "Lung", folder_name)), "_", pc, "PCs.rds"))
      }
      if (!file.exists(dr_file)) {
        dr_file <- file.path(res_path, paste0("DR-SC_", folder_name, "_clustering_results.rds"))
      }

      if (file.exists(dr_file)) {
        tryCatch({
          res_df <- readRDS(dr_file)

          # Check if file already has truth column (STARmap/Lung/Colon/osmFISH format with row/col)
          if ("truth" %in% colnames(res_df) && "label" %in% colnames(res_df) && "row" %in% colnames(res_df)) {
            valid_rows <- !is.na(res_df$truth)
            if (sum(valid_rows) > 0) {
              nmi_val <- NMI(res_df$truth[valid_rows], res_df$label[valid_rows])
              target_mat[d_name, "DR-SC"] <- nmi_val
              message(paste("      DR-SC NMI:", round(nmi_val, 4), "using label column"))
            }
          } else {
            k_cols <- grep(paste0("cluster.*pc", pc), colnames(res_df), ignore.case=TRUE, value=TRUE)
            if (length(k_cols) == 0) k_cols <- grep("cluster|label", colnames(res_df), ignore.case=TRUE, value=TRUE)

            if (length(k_cols) > 0) {
              cluster_col <- k_cols[length(k_cols)]

              if ("imagerow" %in% colnames(res_df)) {
                result <- data.frame(x = res_df$imagerow, y = res_df$imagecol, pred_label = res_df[[cluster_col]])
              } else {
                result <- data.frame(x = res_df$x, y = res_df$y, pred_label = res_df[[cluster_col]])
              }

              m <- merge(gt, result, by=c("x","y"))
              if (nrow(m) > 0) {
                nmi_val <- NMI(m$true_label, m$pred_label)
                target_mat[d_name, "DR-SC"] <- nmi_val
                message(paste("      DR-SC NMI:", round(nmi_val, 4), "using", cluster_col))
              } else {
                message("      DR-SC: No matching observations")
              }
            }
          }
        }, error = function(e) message(paste("      DR-SC error:", e$message)))
      } else {
        message("      DR-SC: File not found")
      }

      # --- k-means ---
      message("    Loading k-means...")
      km_file <- file.path(res_path, paste0("KMeans_", gsub("_Mouse", "", gsub("Starmap_Mouse", "STARmap", gsub("Lung_Xenium", "Lung", gsub("Visium_HD_Human_Colon_Cancer", "Colon", folder_name)))), "_", pc, "PCs.rds"))
      if (!file.exists(km_file)) {
        km_file <- file.path(res_path, paste0("kmeans_", folder_name, "_clustering_results.rds"))
      }

      if (file.exists(km_file)) {
        tryCatch({
          res_df <- readRDS(km_file)

          # Check if file already has truth column (Lung/Colon format)
          if ("truth" %in% colnames(res_df) && "label" %in% colnames(res_df)) {
            valid_rows <- !is.na(res_df$truth)
            if (sum(valid_rows) > 0) {
              nmi_val <- NMI(res_df$truth[valid_rows], res_df$label[valid_rows])
              target_mat[d_name, "k-means"] <- nmi_val
              message(paste("      k-means NMI:", round(nmi_val, 4), "using label column"))
            }
          } else {
            k_cols <- grep(paste0("pc", pc, "_k"), colnames(res_df), ignore.case=TRUE, value=TRUE)
            if (length(k_cols) == 0) k_cols <- grep("cluster|label", colnames(res_df), ignore.case=TRUE, value=TRUE)

            if (length(k_cols) > 0) {
              cluster_col <- k_cols[length(k_cols)]
              result <- data.frame(x = res_df$x, y = res_df$y, pred_label = res_df[[cluster_col]])

              m <- merge(gt, result, by=c("x","y"))
              if (nrow(m) > 0) {
                nmi_val <- NMI(m$true_label, m$pred_label)
                target_mat[d_name, "k-means"] <- nmi_val
                message(paste("      k-means NMI:", round(nmi_val, 4), "using", cluster_col))
              } else {
                message("      k-means: No matching observations")
              }
            }
          }
        }, error = function(e) message(paste("      k-means error:", e$message)))
      } else {
        message("      k-means: File not found")
      }

      # --- SEDR ---
      message("    Loading SEDR...")
      sedr_files <- list.files(res_path, pattern = "SEDR.*\\.csv", full.names = T)
      if (length(sedr_files) > 0) {
        tryCatch({
          res_df <- read.csv(sedr_files[1])
          k_cols <- grep(paste0("label.*pca", pc, "_k|label.*pc", pc), colnames(res_df), ignore.case=TRUE, value=TRUE)
          if (length(k_cols) == 0) k_cols <- grep("label", colnames(res_df), value=TRUE)

          if (length(k_cols) > 0) {
            cluster_col <- k_cols[length(k_cols)]

            # Check if true_label or ground_truth already exists in results file
            if ("true_label" %in% colnames(res_df)) {
              valid_rows <- !is.na(res_df$true_label)
              if (sum(valid_rows) > 0) {
                nmi_val <- NMI(res_df$true_label[valid_rows], res_df[[cluster_col]][valid_rows])
                target_mat[d_name, "SEDR"] <- nmi_val
                message(paste("      SEDR NMI:", round(nmi_val, 4), "using", cluster_col))
              } else {
                message("      SEDR: No valid true_label")
              }
            } else if ("ground_truth" %in% colnames(res_df)) {
              valid_rows <- !is.na(res_df$ground_truth)
              if (sum(valid_rows) > 0) {
                nmi_val <- NMI(res_df$ground_truth[valid_rows], res_df[[cluster_col]][valid_rows])
                target_mat[d_name, "SEDR"] <- nmi_val
                message(paste("      SEDR NMI:", round(nmi_val, 4), "using", cluster_col))
              } else {
                message("      SEDR: No valid ground_truth")
              }
            } else {
              # Need to merge - swap x/y for scanpy-based methods when using RData ground truth
              if (d_name %in% c("Gut", "Brain", "Breast", "Prostate")) {
                result <- data.frame(x = res_df$y, y = res_df$x, pred_label = res_df[[cluster_col]])
              } else {
                result <- data.frame(x = res_df$x, y = res_df$y, pred_label = res_df[[cluster_col]])
              }

              m <- merge(gt, result, by=c("x","y"))
              if (nrow(m) > 0) {
                nmi_val <- NMI(m$true_label, m$pred_label)
                target_mat[d_name, "SEDR"] <- nmi_val
                message(paste("      SEDR NMI:", round(nmi_val, 4), "using", cluster_col))
              } else {
                message("      SEDR: No matching observations")
              }
            }
          }
        }, error = function(e) message(paste("      SEDR error:", e$message)))
      } else {
        message("      SEDR: File not found")
      }

      # --- GraphST ---
      message("    Loading GraphST...")
      gst_files <- list.files(res_path, pattern = "GraphST.*\\.csv", full.names = T)
      if (length(gst_files) > 0) {
        tryCatch({
          res_df <- read.csv(gst_files[1])
          # For osmFISH, use PC-specific column
          cluster_col <- paste0("mclust_", pc, "PCs")
          if (!(cluster_col %in% colnames(res_df))) {
            # Fallback for other datasets
            k_cols <- grep(paste0("mclust.*", pc, "PC|mclust"), colnames(res_df), ignore.case=TRUE, value=TRUE)
            if (length(k_cols) == 0) k_cols <- grep("cluster|label", colnames(res_df), value=TRUE)
            if (length(k_cols) > 0) cluster_col <- k_cols[1] else cluster_col <- NULL
          }

          if (!is.null(cluster_col) && cluster_col %in% colnames(res_df)) {
            # For osmFISH, use barcode-based matching
            if (d_name == "osmFISH") {
              result <- data.frame(barcode = res_df$barcode, pred_label = res_df[[cluster_col]])
              m <- merge_by_barcode(gt, result, gt_barcode_col = "ID", res_barcode_col = "barcode")
            } else {
              # Swap x/y for scanpy methods with RData ground truth
              if (d_name %in% c("Gut", "Brain", "Breast", "Prostate")) {
                result <- data.frame(x = res_df$y, y = res_df$x, pred_label = res_df[[cluster_col]])
              } else {
                result <- data.frame(x = res_df$x, y = res_df$y, pred_label = res_df[[cluster_col]])
              }
              m <- merge(gt, result, by=c("x","y"))
            }

            if (!is.null(m) && nrow(m) > 0) {
              nmi_val <- NMI(m$true_label, m$pred_label)
              target_mat[d_name, "GraphST"] <- nmi_val
              message(paste("      GraphST NMI:", round(nmi_val, 4), "using", cluster_col))
            } else {
              message("      GraphST: No matching observations")
            }
          }
        }, error = function(e) message(paste("      GraphST error:", e$message)))
      } else {
        message("      GraphST: File not found")
      }

      # --- STAGATE ---
      message("    Loading STAGATE...")
      sta_files <- list.files(res_path, pattern = "STAGATE.*\\.csv", full.names = T)
      if (length(sta_files) > 0) {
        tryCatch({
          res_df <- read.csv(sta_files[1])
          cluster_col <- if ("mclust" %in% colnames(res_df)) "mclust" else "cluster"

          if (cluster_col %in% colnames(res_df)) {
            # For osmFISH, use barcode-based matching
            if (d_name == "osmFISH") {
              result <- data.frame(barcode = res_df$barcode, pred_label = res_df[[cluster_col]])
              m <- merge_by_barcode(gt, result, gt_barcode_col = "ID", res_barcode_col = "barcode")
            } else {
              # Swap x/y for scanpy methods with RData ground truth
              if (d_name %in% c("Gut", "Brain", "Breast", "Prostate")) {
                result <- data.frame(x = res_df$y, y = res_df$x, pred_label = res_df[[cluster_col]])
              } else {
                result <- data.frame(x = res_df$x, y = res_df$y, pred_label = res_df[[cluster_col]])
              }
              m <- merge(gt, result, by=c("x","y"))
            }

            if (!is.null(m) && nrow(m) > 0) {
              nmi_val <- NMI(m$true_label, m$pred_label)
              target_mat[d_name, "STAGATE"] <- nmi_val
              message(paste("      STAGATE NMI:", round(nmi_val, 4)))
            } else {
              message("      STAGATE: No matching observations")
            }
          }
        }, error = function(e) message(paste("      STAGATE error:", e$message)))
      } else {
        message("      STAGATE: File not found")
      }

      # --- SpaGCN (with histology) ---
      message("    Loading SpaGCN (w/)...")
      spa_w_files <- list.files(res_path, pattern = "SpaGCN.*histology\\.csv", full.names = T)
      if (length(spa_w_files) > 0) {
        tryCatch({
          res_df <- read.csv(spa_w_files[1])
          k_cols <- grep("refined_pred", colnames(res_df), value=TRUE)

          if (length(k_cols) > 0) {
            cluster_col <- k_cols[1]
            # Swap x/y for scanpy methods with RData ground truth
            if (d_name %in% c("Gut", "Brain", "Breast", "Prostate")) {
              result <- data.frame(x = res_df$y, y = res_df$x, pred_label = res_df[[cluster_col]] + 1)
            } else {
              result <- data.frame(x = res_df$x, y = res_df$y, pred_label = res_df[[cluster_col]] + 1)
            }

            m <- merge(gt, result, by=c("x","y"))
            if (nrow(m) > 0) {
              nmi_val <- NMI(m$true_label, m$pred_label)
              target_mat[d_name, "SpaGCN (w/)"] <- nmi_val
              message(paste("      SpaGCN (w/) NMI:", round(nmi_val, 4)))
            } else {
              message("      SpaGCN (w/): No matching observations")
            }
          }
        }, error = function(e) message(paste("      SpaGCN (w/) error:", e$message)))
      } else {
        message("      SpaGCN (w/): File not found")
      }

      # --- SpaGCN (without histology) ---
      message("    Loading SpaGCN (w/o)...")
      spa_files <- list.files(res_path, pattern = "SpaGCN.*\\.csv", full.names = T)
      spa_wo_files <- spa_files[!grepl("histology", spa_files)]
      if (length(spa_wo_files) > 0) {
        tryCatch({
          res_df <- read.csv(spa_wo_files[1])
          # For osmFISH, use PC-specific column
          cluster_col <- paste0("refined_pred_", pc, "PCs_k11")
          if (!(cluster_col %in% colnames(res_df))) {
            # Fallback for other datasets
            k_cols <- grep("refined_pred", colnames(res_df), value=TRUE)
            if (length(k_cols) > 0) cluster_col <- k_cols[1] else cluster_col <- NULL
          }

          if (!is.null(cluster_col) && cluster_col %in% colnames(res_df)) {
            # For osmFISH, use barcode-based matching
            if (d_name == "osmFISH") {
              result <- data.frame(barcode = res_df$barcode, pred_label = res_df[[cluster_col]] + 1)
              m <- merge_by_barcode(gt, result, gt_barcode_col = "ID", res_barcode_col = "barcode")
            } else {
              # Swap x/y for scanpy methods with RData ground truth
              if (d_name %in% c("Gut", "Brain", "Breast", "Prostate")) {
                result <- data.frame(x = res_df$y, y = res_df$x, pred_label = res_df[[cluster_col]] + 1)
              } else {
                result <- data.frame(x = res_df$x, y = res_df$y, pred_label = res_df[[cluster_col]] + 1)
              }
              m <- merge(gt, result, by=c("x","y"))
            }

            if (!is.null(m) && nrow(m) > 0) {
              nmi_val <- NMI(m$true_label, m$pred_label)
              target_mat[d_name, "SpaGCN (w/o)"] <- nmi_val
              message(paste("      SpaGCN (w/o) NMI:", round(nmi_val, 4)))
            } else {
              message("      SpaGCN (w/o): No matching observations")
            }
          }
        }, error = function(e) message(paste("      SpaGCN (w/o) error:", e$message)))
      } else {
        message("      SpaGCN (w/o): File not found")
      }

      # --- stLearn (with image) ---
      message("    Loading stLearn (w/)...")
      stl_w_files <- list.files(res_path, pattern = "stLearn.*with.*image\\.csv", full.names = T)
      if (length(stl_w_files) > 0) {
        tryCatch({
          res_df <- read.csv(stl_w_files[1])
          k_cols <- grep(paste0("label.*", pc, "PC"), colnames(res_df), ignore.case=TRUE, value=TRUE)
          if (length(k_cols) == 0) k_cols <- grep("label", colnames(res_df), value=TRUE)

          if (length(k_cols) > 0) {
            cluster_col <- k_cols[length(k_cols)]
            # Swap x/y for scanpy methods with RData ground truth
            if (d_name %in% c("Gut", "Brain", "Breast", "Prostate")) {
              result <- data.frame(x = res_df$y, y = res_df$x, pred_label = res_df[[cluster_col]] + 1)
            } else {
              result <- data.frame(x = res_df$x, y = res_df$y, pred_label = res_df[[cluster_col]] + 1)
            }

            m <- merge(gt, result, by=c("x","y"))
            if (nrow(m) > 0) {
              nmi_val <- NMI(m$true_label, m$pred_label)
              target_mat[d_name, "stLearn (w/)"] <- nmi_val
              message(paste("      stLearn (w/) NMI:", round(nmi_val, 4), "using", cluster_col))
            } else {
              message("      stLearn (w/): No matching observations")
            }
          }
        }, error = function(e) message(paste("      stLearn (w/) error:", e$message)))
      } else {
        message("      stLearn (w/): File not found")
      }

      # --- stLearn (without image) ---
      message("    Loading stLearn (w/o)...")
      stl_wo_files <- list.files(res_path, pattern = "stLearn.*without.*image\\.csv", full.names = T)
      if (length(stl_wo_files) == 0) {
        # Try generic stLearn file
        stl_files <- list.files(res_path, pattern = "stLearn.*\\.csv", full.names = T)
        stl_wo_files <- stl_files[!grepl("with", stl_files)]
      }

      if (length(stl_wo_files) > 0) {
        tryCatch({
          res_df <- read.csv(stl_wo_files[1])
          # For osmFISH, use PC-specific column
          cluster_col <- paste0("label_", pc, "PCs_k11")
          if (!(cluster_col %in% colnames(res_df))) {
            # Fallback for other datasets
            k_cols <- grep(paste0("label.*", pc, "PC"), colnames(res_df), ignore.case=TRUE, value=TRUE)
            if (length(k_cols) == 0) k_cols <- grep("label", colnames(res_df), value=TRUE)
            if (length(k_cols) > 0) cluster_col <- k_cols[length(k_cols)] else cluster_col <- NULL
          }

          if (!is.null(cluster_col) && cluster_col %in% colnames(res_df)) {
            # For osmFISH, use barcode-based matching
            if (d_name == "osmFISH") {
              result <- data.frame(barcode = res_df$ID, pred_label = res_df[[cluster_col]] + 1)
              m <- merge_by_barcode(gt, result, gt_barcode_col = "ID", res_barcode_col = "barcode")
            } else {
              # Swap x/y for scanpy methods with RData ground truth
              if (d_name %in% c("Gut", "Brain", "Breast", "Prostate")) {
                result <- data.frame(x = res_df$y, y = res_df$x, pred_label = res_df[[cluster_col]] + 1)
              } else {
                result <- data.frame(x = res_df$x, y = res_df$y, pred_label = res_df[[cluster_col]] + 1)
              }
              m <- merge(gt, result, by=c("x","y"))
            }

            if (!is.null(m) && nrow(m) > 0) {
              nmi_val <- NMI(m$true_label, m$pred_label)
              if (is.na(target_mat[d_name, "stLearn (w/)"])) {
                target_mat[d_name, "stLearn (w/)"] <- nmi_val
                message(paste("      stLearn NMI:", round(nmi_val, 4), "using", cluster_col, "(generic)"))
              } else {
                target_mat[d_name, "stLearn (w/o)"] <- nmi_val
                message(paste("      stLearn (w/o) NMI:", round(nmi_val, 4), "using", cluster_col))
              }
            } else {
              message("      stLearn (w/o): No matching observations")
            }
          }
        }, error = function(e) message(paste("      stLearn (w/o) error:", e$message)))
      } else {
        message("      stLearn (w/o): File not found")
      }

      # --- BASS ---
      message("    Loading BASS...")
      bass_files <- list.files(res_path, pattern = paste0("BASS.*", pc, "PCs?\\.rds"), full.names = T, ignore.case = T)
      if (length(bass_files) > 0) {
        tryCatch({
          res_df <- readRDS(bass_files[1])
          # Check for cluster column - try domain, z, label, cluster
          cluster_col <- NULL
          for (col in c("domain", "z", "label", "cluster")) {
            if (col %in% colnames(res_df)) {
              cluster_col <- col
              break
            }
          }

          if (!is.null(cluster_col) && cluster_col %in% colnames(res_df)) {
            result <- data.frame(x = res_df$x, y = res_df$y, pred_label = res_df[[cluster_col]])

            m <- merge(gt, result, by=c("x","y"))
            if (nrow(m) > 0) {
              nmi_val <- NMI(m$true_label, m$pred_label)
              target_mat[d_name, "BASS"] <- nmi_val
              message(paste("      BASS NMI:", round(nmi_val, 4), "using", cluster_col))
            } else {
              message("      BASS: No matching observations")
            }
          } else {
            message("      BASS: No cluster column found")
          }
        }, error = function(e) message(paste("      BASS error:", e$message)))
      } else {
        message("      BASS: File not found")
      }

      # --- STAMarker ---
      message("    Loading STAMarker...")
      stam_files <- list.files(res_path, pattern = "STAMarker.*joint\\.csv", full.names = T)
      if (length(stam_files) > 0) {
        tryCatch({
          res_df <- read.csv(stam_files[1])
          cluster_col <- if ("STAMarker_label" %in% colnames(res_df)) "STAMarker_label" else "cluster"

          if (cluster_col %in% colnames(res_df)) {
            # Check if file already has truth column (osmFISH format)
            if ("true_label" %in% colnames(res_df)) {
              valid_rows <- !is.na(res_df$true_label)
              if (sum(valid_rows) > 0) {
                nmi_val <- NMI(res_df$true_label[valid_rows], (res_df[[cluster_col]][valid_rows] + 1))
                target_mat[d_name, "STAMarker"] <- nmi_val
                message(paste("      STAMarker NMI:", round(nmi_val, 4), "using pre-merged truth"))
              }
            } else {
              # Swap x/y for scanpy methods with RData ground truth
              if (d_name %in% c("Gut", "Brain", "Breast", "Prostate")) {
                result <- data.frame(x = res_df$y, y = res_df$x, pred_label = res_df[[cluster_col]] + 1)
              } else {
                result <- data.frame(x = res_df$x, y = res_df$y, pred_label = res_df[[cluster_col]] + 1)
              }

              m <- merge(gt, result, by=c("x","y"))
              if (nrow(m) > 0) {
                nmi_val <- NMI(m$true_label, m$pred_label)
                target_mat[d_name, "STAMarker"] <- nmi_val
                message(paste("      STAMarker NMI:", round(nmi_val, 4)))
              } else {
                message("      STAMarker: No matching observations")
              }
            }
          }
        }, error = function(e) message(paste("      STAMarker error:", e$message)))
      } else {
        message("      STAMarker: File not found")
      }

      # --- SpaMask ---
      message("    Loading SpaMask...")
      mask_files <- list.files(res_path, pattern = paste0("SpaMask.*", pc, "PCs?\\.csv"), full.names = T, ignore.case = T)
      if (length(mask_files) > 0) {
        tryCatch({
          res_df <- read.csv(mask_files[1])
          cluster_col <- if ("kmeans" %in% colnames(res_df)) "kmeans" else "cluster"

          if (cluster_col %in% colnames(res_df)) {
            result <- data.frame(x = res_df$x, y = res_df$y, pred_label = res_df[[cluster_col]] + 1)

            m <- merge(gt, result, by=c("x","y"))
            if (nrow(m) > 0) {
              nmi_val <- NMI(m$true_label, m$pred_label)
              target_mat[d_name, "SpaMask"] <- nmi_val
              message(paste("      SpaMask NMI:", round(nmi_val, 4)))
            } else {
              message("      SpaMask: No matching observations")
            }
          }
        }, error = function(e) message(paste("      SpaMask error:", e$message)))
      } else {
        message("      SpaMask: File not found")
      }

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
print(round(results_pc3, 3))
message("\n\nPC = 10 Results:")
print(round(results_pc10, 3))

write.csv(results_pc3, file = file.path(ROOT_DIR, "NMI_CompetingMethods_PC3.csv"), row.names = TRUE)
write.csv(results_pc10, file = file.path(ROOT_DIR, "NMI_CompetingMethods_PC10.csv"), row.names = TRUE)

message("\n\n===== RESULTS SAVED =====")
message("Results saved to:")
message("  - NMI_CompetingMethods_PC3.csv")
message("  - NMI_CompetingMethods_PC10.csv")
