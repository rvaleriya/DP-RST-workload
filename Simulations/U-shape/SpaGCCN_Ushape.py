import sys
import os
import pandas as pd
import numpy as np
import scanpy as sc
import math
from pathlib import Path  # Ensure Path is available
# Ensure we use the local modified SpaGCN package
local_spg_path = Path(__file__).resolve().parent / "SpaGCN" / "SpaGCN_package"
sys.path.insert(0, str(local_spg_path))
import SpaGCN as spg
from scipy.sparse import issparse
import random
import torch
import warnings
import matplotlib.pyplot as plt
import anndata as ad # Import AnnData explicitly
import inspect
from sklearn.metrics import adjusted_rand_score # Import ARI

# Define log filename (will be created in the current directory)
SIM_LOG_FILENAME = "SpaGCN_analysis_ushape_multi_seed.log"

# Setup logging globally
logfile = open(SIM_LOG_FILENAME, "w")
sys.stdout = logfile
sys.stderr = logfile
print(f"--- Script Start: Logging to {SIM_LOG_FILENAME} ---")

def setup_environment(r_seed=42, t_seed=42, n_seed=42):
    """Sets up device and random seeds."""
    # Setup device
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    # Print device only once at the start
    # print(f"Using device: {device}")

    # Set seeds
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    print(f"Seeds set: random={r_seed}, torch={t_seed}, numpy={n_seed}")

    warnings.filterwarnings("ignore")
    return device

def load_simulated_data(csv_path):
    """Loads simulated data from a CSV file and creates an AnnData object."""
    print(f"Loading simulated data from: {csv_path}")
    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        print(f"Error: CSV file not found at {csv_path}")
        return None
    except Exception as e:
        print(f"Error reading CSV file {csv_path}: {e}")
        return None

    print(f"Loaded dataframe with shape: {df.shape}")
    print("Dataframe head:\n", df.head())

    # Check for required columns
    required_cols = ['x', 'y', 'Y1', 'Y2', 'Y3', 'cluster'] # Added 'cluster'
    if not all(col in df.columns for col in required_cols):
        print(f"Error: CSV must contain columns: {required_cols}")
        return None

    # Create AnnData object
    print("Creating AnnData object.")
    # Use Y1, Y2, Y3 as the expression matrix X
    adata = ad.AnnData(X=df[['Y1', 'Y2', 'Y3']].values)
    # Store coordinates in obsm using lowercase x, y
    adata.obsm['spatial'] = df[['x', 'y']].values
    # Use dataframe index as observation names
    adata.obs_names = df.index.astype(str)
    # Optionally store original coordinates in obs if needed later (using lowercase x, y)
    adata.obs['X_coord'] = df['x'].values
    adata.obs['Y_coord'] = df['y'].values
    # Store ground truth clusters
    adata.obs['ground_truth_cluster'] = pd.Categorical(df['cluster'].values)

    print("Created AnnData object:")
    print(adata)
    return adata

def preprocess_simulated_data(adata, p=0.5):
    """Extracts coordinates, calculates adjacency matrix, and finds l for simulated data."""
    # Extract coordinates from obsm
    print("Extracting coordinates from adata.obsm['spatial'].")
    spatial_coords = adata.obsm['spatial']
    x_pixel = spatial_coords[:, 0].tolist()
    y_pixel = spatial_coords[:, 1].tolist()

    # Calculate adjacency matrix based on coordinates only
    print("Calculating adjacency matrix (histology=False).")
    adj = spg.calculate_adj_matrix(x=x_pixel, y=y_pixel, histology=False)

    # Find the l value
    # Note: search_l is deterministic, so l will be the same for all runs unless adj changes
    # If adj calculation had randomness, l could vary. Assuming it's deterministic here.
    print(f"Searching for l parameter with p={p}.")
    l = spg.search_l(p, adj)
    print(f"Found l = {l}")

    # No normalization/feature selection needed as we use pre-processed Y1, Y2, Y3
    print("Skipping normalization and feature selection for simulated data.")

    return adj, l # Return only adj and l, adata is modified inplace or copied in main loop

def run_spagcn_simulated(adata, adj, l, n_clusters, num_pcs, r_seed, t_seed, n_seed):
    """Searches for resolution, trains SpaGCN, predicts, and refines clusters for simulated data.
       Accepts AnnData object and returns it with predictions.
    """
    # Store the input features (Y1, Y2, Y3) in adata.obsm['X_pca']
    # This signals SpaGCN functions (like search_res and train) to use these pre-computed components.
    print("Storing input features (Y1, Y2, Y3) into adata.obsm['X_pca']")
    adata.obsm['X_pca'] = adata.X

    # Search for suitable resolution
    print(f"Searching for resolution for {n_clusters} clusters.")
    # search_res uses adata.obsm['X_pca'] if available
    # Pass seeds to search_res as it calls train internally
    res = spg.search_res(adata, adj, l, n_clusters, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed, num_pcs=num_pcs)
    print(f"Found resolution = {res}")

    # Train SpaGCN
    print("Training SpaGCN model.")
    clf = spg.SpaGCN()
    clf.set_l(l)
    # Train using the features stored in adata.obsm['X_pca']
    # Change init method to louvain to see effect on results
    clf.train(adata, adj, num_pcs=num_pcs, init_spa=True, init="louvain", n_clusters=n_clusters, res=res)
    y_pred, prob = clf.predict()
    adata.obs["pred"] = y_pred
    adata.obs["pred"] = adata.obs["pred"].astype('category')
    print("Initial prediction stored in adata.obs['pred']")

    # Refine clusters
    print("Refining clusters.")
    # Ensure sample_id uses the correct index if it wasn't set properly before
    if adata.obs.index.empty:
         print("Warning: AnnData observation index is empty. Using default range.")
         sample_ids = [str(i) for i in range(adata.n_obs)]
    else:
         sample_ids = adata.obs.index.tolist()

    refined_pred = spg.refine(sample_id=sample_ids, pred=adata.obs["pred"].tolist(), dis=adj, shape="square") # Adjust shape if needed
    adata.obs["refined_pred"] = refined_pred
    adata.obs["refined_pred"] = adata.obs["refined_pred"].astype('category')
    print("Refined prediction stored in adata.obs['refined_pred']")

    return adata # Return the AnnData object with predictions

def main():
    """Main execution function for simulated data multi-seed run."""
    # --- Configuration ---
    SIM_DATA_CSV = "/scratch/user/varogovchenko/BASTION_HPRC/Simulations/U_sim_6k/sim_data_ushape.csv"
    N_CLUSTERS = 6 # Number of clusters requested (used for kmeans init)
    NUM_PCS = 3 # Since Y1, Y2, Y3 are the features used
    OUTPUT_DIR = Path("./res_ushape_data")
    P_PARAM = 0.5 # p parameter for search_l

    seeds_to_run = [
        141, 549, 75, 492, 676, 179, 587, 592, 601, 916, 518, 339, 921, 423, 330,
        388, 273, 286, 61, 807, 283, 127, 165, 952, 311, 597, 473, 594, 605, 101
    ]
    # ---------------------

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Log initial info
    print(f"SpaGCN Version: {spg.__version__}")
    # Setup device once
    device = setup_environment() # Use default seeds initially, just to get device
    print(f"Using device: {device}")

    # 1. Load Simulated Data (once)
    adata_master = load_simulated_data(SIM_DATA_CSV)
    if adata_master is None:
        print("Exiting due to data loading error.")
        sys.exit(1) # Exit if data loading failed

    # Extract ground truth and coordinates for the results dataframe
    ground_truth_clusters = adata_master.obs['ground_truth_cluster'].astype(str).tolist() # Ensure comparable type
    all_results_df = pd.DataFrame({
        'X_coord': adata_master.obs['X_coord'],
        'Y_coord': adata_master.obs['Y_coord'],
        'ground_truth_cluster': adata_master.obs['ground_truth_cluster']
    }, index=adata_master.obs_names)

    # 2. Preprocessing (once, assuming deterministic adj and l)
    # If adj/l calculation depends on seeds, move this inside the loop
    try:
        print("\n--- Performing Initial Preprocessing ---")
        adj, l = preprocess_simulated_data(adata_master, p=P_PARAM)
    except Exception as e:
        print(f"Error during initial preprocessing: {e}. Exiting.")
        if logfile and not logfile.closed:
             print(f"--- Script End (Error): Log file {SIM_LOG_FILENAME} closed ---")
             logfile.close()
        sys.exit(1)

    # --- Loop through seeds ---
    ari_scores_initial = [] # Store ARI before refinement
    ari_scores_refined = [] # Store ARI after refinement
    print(f"\n--- Starting Multi-Seed Analysis ({len(seeds_to_run)} seeds) ---")
    for i, seed in enumerate(seeds_to_run):
        print(f"\n{'='*20} Running Seed {i+1}/{len(seeds_to_run)} (Seed Value: {seed}) {'='*20}")

        # Set seeds for this run
        _ = setup_environment(r_seed=seed, t_seed=seed, n_seed=seed)

        # Create a copy for this run to avoid interference
        adata_run = adata_master.copy()

        # 3. Run SpaGCN (uses pre-calculated adj, l)
        try:
            adata_run = run_spagcn_simulated(adata_run, adj, l, N_CLUSTERS, NUM_PCS, r_seed=seed, t_seed=seed, n_seed=seed)

            # 4. Calculate ARI for this run (Initial and Refined)
            initial_pred = adata_run.obs['pred'].astype(str).tolist()
            refined_pred = adata_run.obs['refined_pred'].astype(str).tolist()

            ari_initial = float('nan') # Default to NaN
            ari_refined = float('nan') # Default to NaN

            if len(ground_truth_clusters) == len(initial_pred):
                try:
                    ari_initial = adjusted_rand_score(ground_truth_clusters, initial_pred)
                    print(f"ARI Initial (seed {seed}): {ari_initial:.4f}")
                    ari_scores_initial.append(ari_initial)
                except ValueError as e:
                    print(f"Error calculating Initial ARI for seed {seed}: {e}")
            else:
                 print(f"Error: Mismatch in length for Initial ARI (seed {seed}).")

            if len(ground_truth_clusters) == len(refined_pred):
                try:
                    ari_refined = adjusted_rand_score(ground_truth_clusters, refined_pred)
                    print(f"ARI Refined (seed {seed}): {ari_refined:.4f}")
                    ari_scores_refined.append(ari_refined)
                except ValueError as e:
                    print(f"Error calculating Refined ARI for seed {seed}: {e}")
            else:
                 print(f"Error: Mismatch in length for Refined ARI (seed {seed}).")

            # Add results to the main DataFrame
            all_results_df[f'initial_pred_seed_{seed}'] = adata_run.obs['pred']
            all_results_df[f'refined_pred_seed_{seed}'] = adata_run.obs['refined_pred']

        except Exception as e:
            print(f"Error during SpaGCN run for seed {seed}: {e}. Skipping seed.")
            # Add NaN columns to indicate failure for this seed
            all_results_df[f'initial_pred_seed_{seed}'] = np.nan
            all_results_df[f'refined_pred_seed_{seed}'] = np.nan
            continue # Continue to the next seed

    # --- Post-Loop Analysis ---
    print(f"\n{'='*20} Multi-Seed Analysis Complete {'='*20}")

    # 5. Calculate and Print ARI Statistics (Initial and Refined)
    print("\nARI Statistics Across Seeds:")
    if ari_scores_initial:
        ari_initial_np = np.array(ari_scores_initial)
        print("  Initial Prediction ARI:")
        print(f"    Mean: {np.mean(ari_initial_np):.4f}, StdDev: {np.std(ari_initial_np):.4f}")
        print(f"    Min:  {np.min(ari_initial_np):.4f}, Max:  {np.max(ari_initial_np):.4f}")
        print(f"    Successful runs: {len(ari_initial_np)}")
    else:
        print("  Initial Prediction ARI: No successful runs.")

    if ari_scores_refined:
        ari_refined_np = np.array(ari_scores_refined)
        print("  Refined Prediction ARI:")
        print(f"    Mean: {np.mean(ari_refined_np):.4f}, StdDev: {np.std(ari_refined_np):.4f}")
        print(f"    Min:  {np.min(ari_refined_np):.4f}, Max:  {np.max(ari_refined_np):.4f}")
        print(f"    Successful runs: {len(ari_refined_np)}")
    else:
        print("  Refined Prediction ARI: No successful runs.")

    # 7. Generate and Save ARI Boxplot (Initial vs Refined)
    if ari_scores_initial or ari_scores_refined:
        plot_data = []
        plot_labels = []
        if ari_scores_initial:
            plot_data.append(ari_scores_initial)
            plot_labels.append('Initial')
        if ari_scores_refined:
            plot_data.append(ari_scores_refined)
            plot_labels.append('Refined')

        if plot_data:
            try:
                fig, ax = plt.subplots(figsize=(8, 6))
                bp = ax.boxplot(plot_data, patch_artist=True, labels=plot_labels)
                # Optional: customize box colors
                colors = ['lightblue', 'lightgreen']
                for patch, color in zip(bp['boxes'], colors):
                    patch.set_facecolor(color)
                ax.set_ylabel("Adjusted Rand Index (ARI)")
                ax.set_title(f"ARI Distribution Across {len(seeds_to_run)} Seeds")
                ax.grid(axis='y', linestyle='--', alpha=0.7)
                plt.tight_layout()
                boxplot_path = OUTPUT_DIR / "simulated_ari_boxplot_comparison.png"
                plt.savefig(boxplot_path, dpi=300)
                print(f"Saved ARI comparison boxplot to: {boxplot_path}")
                plt.close(fig)
            except Exception as e:
                print(f"Error generating or saving ARI comparison boxplot: {e}")

    # 8. Save combined results DataFrame to CSV
    results_csv_path = OUTPUT_DIR / "simulated_clustering_results_multi_seed.csv"
    print(f"\nSaving all clustering results to: {results_csv_path}")
    try:
        all_results_df.to_csv(results_csv_path)
    except Exception as e:
        print(f"Error saving CSV {results_csv_path}: {e}")

    print("\nMulti-seed simulated data analysis complete.")

    # Debug prints
    try:
        print("SpaGCN module file:", spg.__file__)
        print("search_res function file:", inspect.getsourcefile(spg.search_res))
    except Exception as e:
        print(f"Error getting module/function file info: {e}")


if __name__ == "__main__":
    main()
    # Close the log file
    if logfile and not logfile.closed: # Check if already closed by error handling
        print(f"--- Script End: Log file {SIM_LOG_FILENAME} closed ---")
        logfile.close() 