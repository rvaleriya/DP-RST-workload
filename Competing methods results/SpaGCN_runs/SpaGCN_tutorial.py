import sys
import os
import csv
import re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random
import torch
import warnings
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import cv2
from pathlib import Path

# Define log filename (will be created in the current directory)
SINGLE_LOG_FILENAME = "SpaGCN_analysis_all_datasets.log" 

# Setup logging globally
logfile = open(SINGLE_LOG_FILENAME, "w")
sys.stdout = logfile
sys.stderr = logfile
print(f"--- Script Start: Logging to {SINGLE_LOG_FILENAME} ---")


def setup_environment(r_seed=42, t_seed=42, n_seed=42):
    """Sets up device and random seeds."""
    # print(f"SpaGCN Version: {spg.__version__}") # Logged elsewhere

    # Setup device
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    # print(f"Using device: {device}") # Logged elsewhere

    # Set seeds
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    # print(f"Seeds set: random={r_seed}, torch={t_seed}, numpy={n_seed}") # Logged elsewhere

    warnings.filterwarnings("ignore")
    return device, r_seed, t_seed, n_seed

def load_data(data_path, count_file, image_path):
    """Loads Visium data and histology image."""
    print(f"Loading data from: {data_path}")
    adata = sc.read_visium(data_path, count_file=count_file)
    print("Initial AnnData object:")
    print(adata)

    # Convert sparse matrix to dense matrix if necessary
    if issparse(adata.X):
        print("Converting sparse data matrix to dense.")
        adata.X = adata.X.toarray()

    print(f"Loading image from: {image_path}")
    img = cv2.imread(image_path)
    if img is None:
        print(f"Warning: Could not load image from {image_path}. Proceeding without histology image.")

    adata.var_names_make_unique()
    return adata, img

def preprocess_data(adata, img, use_histology=True, p=0.5, s=1, b=49):
    """Extracts coordinates, calculates adjacency matrix, normalizes data, and finds l."""
    # Extract coordinates
    print("Extracting coordinates.")
    x_array = adata.obs["array_col"].tolist()
    y_array = adata.obs["array_row"].tolist()
    x_pixel = adata.obsm['spatial'][:, 0].tolist()
    y_pixel = adata.obsm['spatial'][:, 1].tolist()

    # Calculate adjacency matrix
    print("Calculating adjacency matrix.")
    if use_histology and img is not None:
        print("Using histology image for adjacency calculation.")
        adj = spg.calculate_adj_matrix(x=x_pixel, y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, image=img, beta=b, alpha=s, histology=True)
    else:
        print("Calculating adjacency matrix based on coordinates only.")
        adj = spg.calculate_adj_matrix(x=x_pixel, y=y_pixel, histology=False)

    # Normalize and log transform
    print("Normalizing and log-transforming data.")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable].copy()
    print("A subset of data with 2000 highly variable genes:")
    print(adata)

    # Find the l value
    print(f"Searching for l parameter with p={p}.")
    l = spg.search_l(p, adj)
    print(f"Found l = {l}")

    return adata, adj, l, x_array, y_array

def run_spagcn(adata, adj, l, n_clusters, num_pcs, r_seed, t_seed, n_seed):
    """Searches for resolution, trains SpaGCN, predicts, and refines clusters."""
    # Search for suitable resolution
    print(f"Searching for resolution for {n_clusters} clusters.")
    res = spg.search_res(adata, adj, l, n_clusters, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)
    print(f"Found resolution = {res}")

    # Train SpaGCN
    print("Training SpaGCN model.")
    clf = spg.SpaGCN()
    clf.set_l(l)
    clf.train(adata, adj, num_pcs=num_pcs, init_spa=True, init="kmeans", n_clusters=n_clusters, res=res)
    y_pred, prob = clf.predict()
    adata.obs["pred"] = y_pred
    adata.obs["pred"] = adata.obs["pred"].astype('category')
    print("Initial prediction stored in adata.obs['pred']")

    # Refine clusters
    print("Refining clusters.")
    refined_pred = spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj, shape="hexagon")
    adata.obs["refined_pred"] = refined_pred
    adata.obs["refined_pred"] = adata.obs["refined_pred"].astype('category')
    print("Refined prediction stored in adata.obs['refined_pred']")

    return adata

def visualize_results(adata, output_dir="./", filename="spatial_refined_prediction.png", spot_size=250):
    """Plots the refined spatial domains and saves the figure."""
    print(f"Generating visualization: {filename}")
    output_path = Path(output_dir) / filename
    n_clusters_found = len(adata.obs["refined_pred"].cat.categories)

    fig, ax = plt.subplots(figsize=(8, 5))
    sc.pl.spatial(adata, color=["refined_pred"], spot_size=spot_size,
                  title=[f'Refined Prediction (K={n_clusters_found})'], show=False, ax=ax)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Plot saved to {output_path}")

def main():
    """Main execution function."""
    # --- Configuration ---
    BASE_DATA_DIR = Path("/ST_Datasets") # Base directory for all datasets
    COUNT_FILE_NAME = "filtered_feature_bc_matrix.h5" # Standard count file name
    IMAGE_FILE_NAME = "spatial/tissue_tiff_image.tif" # Standard image file name relative to dataset dir
    PLOT_OUTPUT_DIR = Path("./res_realdata_plots") # Directory for saving plots
    CSV_OUTPUT_DIR = Path("./res_realdata_csv") # Directory for saving CSV results
    NUM_PCS_LIST = [3, 10]  # List of PC numbers to try

    # Define datasets and cluster numbers to try per dataset
    dataset_config = {
        "Brain": [5, 7, 9],
        "Breast": [3, 5, 7],
        "Gut": [3, 5, 7, 9, 12],
        "Gut_reduced": [3, 5, 7, 9, 12],
        "Prostate": [3, 5, 7]  
    }
    # ---------------------

    # Create output directories (excluding logs)
    PLOT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    CSV_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    # LOG_FILE_DIR.mkdir(parents=True, exist_ok=True) # Removed

    # Log initial info (will go to the globally redirected stdout/stderr)
    print(f"SpaGCN Version: {spg.__version__}")

    # 1. Setup device and seeds (function no longer handles logging)
    device, r_seed, t_seed, n_seed = setup_environment(r_seed=42, t_seed=42, n_seed=42)
    print(f"Using device: {device}")
    print(f"Seeds set: random={r_seed}, torch={t_seed}, numpy={n_seed}")
    
    # --- Dataset Loop ---
    for dataset_name, n_clusters_list in dataset_config.items():
        print(f"\n{'#'*70}\n# Processing Dataset: {dataset_name}\n{'#'*70}\n")

        # Define dataset-specific paths
        current_data_dir = BASE_DATA_DIR / dataset_name
        current_count_file = COUNT_FILE_NAME
        current_image_path = str(current_data_dir / IMAGE_FILE_NAME) # Keep as string for cv2.imread
        # log_file_path = LOG_FILE_DIR / f"SpaGCN_analysis_{dataset_name}.log" # Removed: Logging setup moved outside loop

        # Check if data directory and count file exist
        if not current_data_dir.is_dir():
            print(f"Error: Data directory not found for dataset {dataset_name} at {current_data_dir}. Skipping dataset.")
            continue
        if not (current_data_dir / current_count_file).is_file():
             print(f"Error: Count file not found for dataset {dataset_name} at {current_data_dir / current_count_file}. Skipping dataset.")
             continue

        # 2. Load Data
        try:
            # Use str() for path compatibility with sc.read_visium if needed, Path works too
            adata, img = load_data(str(current_data_dir), current_count_file, current_image_path)
            img_loaded = img is not None
        except Exception as e:
            print(f"Error loading data for {dataset_name}: {e}. Skipping dataset.")
            continue # Move to the next dataset

        # --- Add Gut_reduced specific subsetting --- Start ---
        if dataset_name == "Gut_reduced":
            print("Subsetting Gut_reduced to selected spots...")
            subset_csv_path = current_data_dir / "gut_df_wt_muscle_rownames.csv" # Use Path object
            if subset_csv_path.is_file():
                try:
                    spot_ids = pd.read_csv(subset_csv_path, header=None).squeeze().astype(str)
                    # Ensure spot IDs from CSV exist in adata
                    valid_spots = spot_ids[spot_ids.isin(adata.obs_names)]
                    if len(valid_spots) < len(spot_ids):
                         print(f"Warning: {len(spot_ids) - len(valid_spots)} spot IDs from CSV not found in adata for Gut_reduced.")
                    if len(valid_spots) > 0:
                        adata = adata[valid_spots].copy()
                        print(f"✓ Subsetted to {adata.n_obs} spots")
                        print(adata)
                    else:
                        print("Error: No valid spots found after subsetting Gut_reduced. Skipping dataset.")
                        continue # Skip this dataset
                except Exception as e:
                     print(f"Error reading or applying subset CSV {subset_csv_path}: {e}. Skipping Gut_reduced subsetting.")
                     # Decide whether to continue with the full dataset or skip
                     print("Proceeding with the full Gut_reduced dataset instead.")
            else:
                print(f"Warning: Subset CSV not found at {subset_csv_path}. Proceeding with full Gut_reduced dataset.")
        # --- Add Gut_reduced specific subsetting --- End ---

        # Initialize DataFrame to store all results for this dataset
        print("Initializing results DataFrame.")
        results_df = pd.DataFrame(index=adata.obs.index)
        results_df['array_row'] = adata.obs['array_row']
        results_df['array_col'] = adata.obs['array_col']

        # Run analysis for both with and without histology
        for use_histology in [True, False]:
            
            # Skip 'with histology' if image wasn't loaded
            if use_histology and not img_loaded:
                print(f"Skipping 'with histology' analysis for {dataset_name} as image was not loaded.")
                continue

            print(f"\n{'='*50}")
            print(f"Running analysis {'with' if use_histology else 'without'} histology image")
            print(f"{'='*50}\n")

            # Create a copy of the AnnData object for this run
            adata_histology_copy = adata.copy() # Use a distinct name

            # 3. Preprocessing
            try:
                adata_histology_copy, adj, l, x_array, y_array = preprocess_data(adata_histology_copy, img, use_histology=use_histology)
            except Exception as e:
                 print(f"Error during preprocessing for {dataset_name} ({'with' if use_histology else 'without'} histology): {e}. Skipping this run.")
                 continue # Skip to next histology setting or dataset

            # 4. Run SpaGCN for each combination
            for num_pcs in NUM_PCS_LIST:
                for n_clusters in n_clusters_list: # Use dataset-specific cluster list
                    print(f"\nRunning analysis with {num_pcs} PCs and {n_clusters} clusters")
                    
                    # Create a fresh copy for this specific combination run
                    adata_combination = adata_histology_copy.copy()
                    
                    try:
                        # Run SpaGCN
                        adata_combination = run_spagcn(adata_combination, adj, l, n_clusters, num_pcs, r_seed, t_seed, n_seed)
                        
                        # Define result name for this combination
                        result_name = f"refined_pred_{n_clusters}clusters_{num_pcs}pcs_{'with' if use_histology else 'without'}_histology"
                        
                        # Add results to the main DataFrame for this dataset
                        results_df[result_name] = adata_combination.obs["refined_pred"]
                        print(f"Added column '{result_name}' to results DataFrame.")

                        # Also store in adata_histology_copy for plotting this round
                        adata_histology_copy.obs[result_name] = adata_combination.obs["refined_pred"]
                        adata_histology_copy.obs[result_name] = adata_histology_copy.obs[result_name].astype('category')
                    except Exception as e:
                         print(f"Error during SpaGCN run for {dataset_name}, PC={num_pcs}, K={n_clusters} ({'with' if use_histology else 'without'} histology): {e}. Skipping this combination.")
                         continue # Skip to next cluster/PC

            # 5. Create comparison plot for this histology setting
            print("\nCreating comparison plot")
            n_cols = len(n_clusters_list) # Dynamic number of columns based on dataset config
            n_rows = len(NUM_PCS_LIST)
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows), squeeze=False) # Ensure axes is always 2D
            histology_label = "with" if use_histology else "without"
            fig.suptitle(f'{dataset_name}: SpaGCN Results ({histology_label} histology)', fontsize=16)

            plot_successful = False
            for i, num_pcs in enumerate(NUM_PCS_LIST):
                for j, n_clusters in enumerate(n_clusters_list):
                    ax = axes[i, j]
                    plot_result_name = f"refined_pred_{n_clusters}clusters_{num_pcs}pcs_{histology_label}_histology"
                    if plot_result_name in adata_histology_copy.obs:
                        try:
                            n_clusters_found = len(adata_histology_copy.obs[plot_result_name].cat.categories)
                            sc.pl.spatial(adata_histology_copy, color=[plot_result_name], spot_size=250,
                                        title=[f'PC={num_pcs}, K={n_clusters_found}'], 
                                        show=False, ax=ax)
                            plot_successful = True
                        except Exception as e:
                             print(f"Error plotting {plot_result_name} for {dataset_name}: {e}")
                             ax.text(0.5, 0.5, f"Error plotting:\n{plot_result_name}", ha='center', va='center', transform=ax.transAxes)
                             ax.set_title(f'PC={num_pcs}, K={n_clusters} (Error)')
                             ax.set_xticks([])
                             ax.set_yticks([])
                    else: # Result missing (likely due to earlier error)
                         ax.text(0.5, 0.5, f"Result not found:\n{plot_result_name}", 
                                        horizontalalignment='center', verticalalignment='center', 
                                        transform=ax.transAxes)
                         ax.set_title(f'PC={num_pcs}, K={n_clusters} (Missing)')
                         ax.set_xticks([])
                         ax.set_yticks([])

            if plot_successful:
                plt.tight_layout(rect=[0, 0.03, 1, 0.95])
                output_filename_png = PLOT_OUTPUT_DIR / f'spatial_comparison_{histology_label}_histology_{dataset_name}.png' # Add dataset name
                try:
                    plt.savefig(output_filename_png, dpi=300, bbox_inches='tight')
                    print(f"Saved comparison plot to {output_filename_png}")
                except Exception as e:
                    print(f"Error saving plot {output_filename_png}: {e}")
            else:
                 print(f"Skipping saving plot for {dataset_name} ({histology_label} histology) as no results were plotted.")
            plt.close(fig) # Close figure

        # 6. Save the combined results DataFrame for this dataset to CSV
        if not results_df.empty and len(results_df.columns) > 2: # Check if results were added
            results_csv_path = CSV_OUTPUT_DIR / f"all_clustering_results_{dataset_name}.csv" # Add dataset name
            print(f"\nSaving all clustering results for {dataset_name} to: {results_csv_path}")
            try:
                results_df.to_csv(results_csv_path)
            except Exception as e:
                print(f"Error saving CSV {results_csv_path}: {e}")
        else:
            print(f"\nNo clustering results generated for {dataset_name}. Skipping CSV save.")

    print("\nAnalysis complete for all datasets.")

if __name__ == "__main__":
    main()
