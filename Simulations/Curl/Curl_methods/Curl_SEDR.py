import os
os.environ['R_HOME'] = '/sw/eb/sw/R/4.4.2-gfbf-2024a/lib64/R'
os.environ['R_LIBS_USER'] = '/scratch/user/varogovchenko/Rlibs'

import scanpy as sc
import pandas as pd
from sklearn import metrics
import torch
import numpy as np # Import numpy

import sys
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import SEDR
import anndata # Import anndata explicitly

# --- Configuration ---
data_file = "/scratch/user/varogovchenko/BASTION_HPRC/Simulations/Curl/Curl_sim_data.csv"
output_dir_csv = "res_curl_3p_csv"
log_file_path = "SEDR_Curl_3p_MultiSeed_output_BestResult.log" # Updated log file name

# List of seeds to run
seeds_to_run = [
    141, 549,  75, 492, 676, 179, 587, 592, 601, 916, 518, 339, 921, 423, 330, 
    388, 273, 286,  61, 807, 283, 127, 165, 952, 311, 597, 473, 594, 605, 101
]

num_clusters = 10 # Based on 'cluster' column in data
graph_neighbors = 30
sedr_train_epochs = 1 # Keep low for quick testing, increase for better results

# --- Setup ---
logfile = open(log_file_path, "w")
sys.stdout = logfile
sys.stderr = logfile

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")
print(f"Using R_HOME: {os.environ.get('R_HOME')}")
print(f"Using R_LIBS_USER: {os.environ.get('R_LIBS_USER')}")

os.makedirs(output_dir_csv, exist_ok=True)

# --- Load Data (once) ---
print(f"Loading data from: {data_file}")
try:
    df = pd.read_csv(data_file)
except Exception as e:
    print(f"Error loading CSV: {e}")
    sys.exit(1)

print("Data loaded successfully:")
print(df.head())

# --- Create AnnData Object (once) ---
print("Creating AnnData object...")
features = df[['PC1', 'PC2', 'PC3']].values
spatial_coords = df[['X', 'Y']].values
observations = df[['super_cluster']].copy() # Keep original labels

adata = anndata.AnnData(X=features.copy())
adata.obsm['spatial'] = spatial_coords.copy()
adata.obs = observations.copy() # Use the copy directly
adata.obs.index = df.index.astype(str) # Use original index as obs_names
adata.var_names = ['Y1', 'Y2', 'Y3'] # Name the features
adata.var_names_make_unique()

print("AnnData object created:")
print(adata)

# --- Preprocessing (Commented out as per previous request) ---
print("Preprocessing data... SKIPPING") # Indicate skipping
# adata.layers['count'] = adata.X.copy() # Store raw features before normalization - No longer needed as X remains raw
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)
# print("Preprocessing complete:")
# print(adata)
# No HVG selection needed as we use Y1, Y2, Y3 directly

# --- Loop Over Seeds ---
all_aris = {} # Dictionary to store ARI for each seed

for seed in seeds_to_run:
    print(f"--- Processing with Random Seed: {seed} ---")
    SEDR.fix_seed(seed)

    # --- Construct Graph (inside loop for reproducibility) ---
    print(f"Constructing neighborhood graph with k={graph_neighbors}...")
    try:
        # Use spatial coordinates directly for graph construction if preferred, or adata if required
        # Assuming graph_construction only needs spatial coordinates and neighbors:
        # graph_dict = SEDR.graph_construction(adata.obsm['spatial'], graph_neighbors) # WRONG - Commented out
        # If it needs the full adata:
        graph_dict = SEDR.graph_construction(adata, graph_neighbors) # CORRECT - Pass full adata object
        print("Graph construction complete.")
    except Exception as e:
        print(f"Error during graph construction for seed {seed}: {e}")
        continue # Skip to next seed

    # --- Train SEDR (inside loop) ---
    # Use the raw features (Y1, Y2, Y3) directly
    sedr_input_features = adata.X 
    print(f"Training SEDR using {sedr_input_features.shape[1]} features for seed {seed}...")
    try:
        # Re-initialize network for each seed to ensure independence
        sedr_net = SEDR.Sedr(sedr_input_features, graph_dict, mode='clustering', device=device)
        using_dec = True # Use DEC by default
        if using_dec:
            print(f"Training with DEC for {sedr_train_epochs} N...")
            sedr_net.train_with_dec(N=sedr_train_epochs)
        else:
            print(f"Training without DEC for {sedr_train_epochs} N...")
            sedr_net.train_without_dec(N=sedr_train_epochs)
        
        sedr_feat, _, _, _ = sedr_net.process()
        # Don't overwrite adata.obsm['SEDR'] each time, maybe store if needed later?
        # For now, just use sedr_feat directly for clustering.
        print(f"✓ SEDR training and processing complete for seed {seed}.")
        # print(f"SEDR embedding shape: {sedr_feat.shape}") # Optional: Check shape
    except Exception as e:
        print(f"Error during SEDR training for seed {seed}: {e}")
        continue # Skip to next seed

    # --- Clustering (inside loop) ---
    print(f"Performing mclust clustering with k={num_clusters} using SEDR embedding for seed {seed}...")
    cluster_key_seed = f'mclust_seed_{seed}' 
    try:
        # We need to provide the SEDR features directly to mclust_R if we don't store them in obsm
        # Modifying mclust_R might be needed, OR we temporarily store it. Let's store temporarily.
        temp_sedr_key = 'SEDR_temp'
        adata.obsm[temp_sedr_key] = sedr_feat
        
        SEDR.mclust_R(adata, num_clusters, use_rep=temp_sedr_key, key_added=cluster_key_seed)
        
        # Clean up temporary embedding
        del adata.obsm[temp_sedr_key] 
        
        print(f"✓ Clustering complete. Results stored in adata.obs['{cluster_key_seed}']")
    except Exception as e:
        print(f"Error during mclust clustering for seed {seed}: {e}")
        # Clean up temporary embedding even if clustering failed
        if temp_sedr_key in adata.obsm:
            del adata.obsm[temp_sedr_key]
        print("Clustering failed. Check R setup and mclust installation.")
        continue # Skip to next seed

    # Calculate ARI if ground truth 'cluster' exists
    if 'super_cluster' in adata.obs.columns:
        try:
            true_labels = adata.obs['super_cluster']
            pred_labels = adata.obs[cluster_key_seed]
            ari = metrics.adjusted_rand_score(true_labels, pred_labels)
            all_aris[seed] = ari
            print(f"Adjusted Rand Index (ARI) vs 'cluster' for seed {seed}: {ari:.4f}")
        except Exception as e:
            print(f"Error calculating ARI for seed {seed}: {e}")
            all_aris[seed] = np.nan # Store NaN if ARI calculation fails
    else:
        if seed == seeds_to_run[0]: # Print only once
             print("Ground truth 'cluster' column not found in adata.obs. Skipping ARI calculation.")

# --- Combine and Save Results ---
print("Consolidating results...")

# Prepare final DataFrame
final_df = pd.DataFrame({
    'barcode': adata.obs.index,
    'x': adata.obsm['spatial'][:, 0],
    'y': adata.obsm['spatial'][:, 1],
    'original_cluster': adata.obs.get('super_cluster', pd.NA), # Use pd.NA for missing originals
})

# Add columns for each seed's clustering result
for seed in seeds_to_run:
    cluster_key_seed = f'mclust_seed_{seed}'
    if cluster_key_seed in adata.obs.columns:
        final_df[cluster_key_seed] = adata.obs[cluster_key_seed]
    else:
        # Add column with NAs if clustering failed for a seed
        final_df[cluster_key_seed] = pd.NA 
        print(f"Warning: Clustering results for seed {seed} not found. Filling with NA.")


csv_filename = os.path.join(output_dir_csv, f"sedr_curl_3p_multiseed_k{num_clusters}_results_BestResult.csv")
final_df.to_csv(csv_filename, index=False)
print(f"✓ Saved consolidated results to {csv_filename}")

# Print summary of ARIs
if all_aris:
    print("--- ARI Summary ---")
    for seed, ari in all_aris.items():
        print(f"Seed {seed}: {ari:.4f}")
    ari_series = pd.Series(all_aris)
    print(f"Mean ARI: {ari_series.mean():.4f}")
    print(f"Std Dev ARI: {ari_series.std():.4f}")
    print(f"Median ARI: {ari_series.median():.4f}")
    print(f"Min ARI: {ari_series.min():.4f} (Seed: {ari_series.idxmin()})")
    print(f"Max ARI: {ari_series.max():.4f} (Seed: {ari_series.idxmax()})")

# --- Finalize ---
print("=== Multi-seed Script finished ===")
# Rerun stdout/stderr to terminal BEFORE closing the logfile
sys.stdout = sys.__stdout__
sys.stderr = sys.__stderr__

logfile.close()

print(f"Output logged to {log_file_path}")

# ============================================================================
# 10 PC VERSION
# ============================================================================

# --- Configuration for 10 PCs ---
output_dir_csv_10p = "res_curl_10p_csv"
log_file_path_10p = "SEDR_Curl_10p_MultiSeed_output_BestResult.log"

# --- Setup for 10 PCs ---
logfile_10p = open(log_file_path_10p, "w")
sys.stdout = logfile_10p
sys.stderr = logfile_10p

print(f"Using device: {device}")
print(f"Using R_HOME: {os.environ.get('R_HOME')}")
print(f"Using R_LIBS_USER: {os.environ.get('R_LIBS_USER')}")

os.makedirs(output_dir_csv_10p, exist_ok=True)

# --- Create AnnData Object for 10 PCs ---
print("Creating AnnData object with 10 PCs...")
features_10p = df[['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']].values
spatial_coords_10p = df[['X', 'Y']].values
observations_10p = df[['super_cluster']].copy()

adata_10p = anndata.AnnData(X=features_10p.copy())
adata_10p.obsm['spatial'] = spatial_coords_10p.copy()
adata_10p.obs = observations_10p.copy()
adata_10p.obs.index = df.index.astype(str)
adata_10p.var_names = [f'Y{i}' for i in range(1, 11)]  # Y1, Y2, ..., Y10
adata_10p.var_names_make_unique()

print("AnnData object created:")
print(adata_10p)

# --- Preprocessing (Skipped) ---
print("Preprocessing data... SKIPPING")

# --- Loop Over Seeds for 10 PCs ---
all_aris_10p = {}  # Dictionary to store ARI for each seed

for seed in seeds_to_run:
    print(f"--- Processing with Random Seed: {seed} (10 PCs) ---")
    SEDR.fix_seed(seed)

    # --- Construct Graph ---
    print(f"Constructing neighborhood graph with k={graph_neighbors}...")
    try:
        graph_dict_10p = SEDR.graph_construction(adata_10p, graph_neighbors)
        print("Graph construction complete.")
    except Exception as e:
        print(f"Error during graph construction for seed {seed}: {e}")
        continue

    # --- Train SEDR ---
    sedr_input_features_10p = adata_10p.X
    print(f"Training SEDR using {sedr_input_features_10p.shape[1]} features for seed {seed}...")
    try:
        sedr_net_10p = SEDR.Sedr(sedr_input_features_10p, graph_dict_10p, mode='clustering', device=device)
        using_dec = True
        if using_dec:
            print(f"Training with DEC for {sedr_train_epochs} N...")
            sedr_net_10p.train_with_dec(N=sedr_train_epochs)
        else:
            print(f"Training without DEC for {sedr_train_epochs} N...")
            sedr_net_10p.train_without_dec(N=sedr_train_epochs)
        
        sedr_feat_10p, _, _, _ = sedr_net_10p.process()
        print(f"✓ SEDR training and processing complete for seed {seed}.")
    except Exception as e:
        print(f"Error during SEDR training for seed {seed}: {e}")
        continue

    # --- Clustering ---
    print(f"Performing mclust clustering with k={num_clusters} using SEDR embedding for seed {seed}...")
    cluster_key_seed_10p = f'mclust_seed_{seed}'
    try:
        temp_sedr_key_10p = 'SEDR_temp'
        adata_10p.obsm[temp_sedr_key_10p] = sedr_feat_10p
        
        SEDR.mclust_R(adata_10p, num_clusters, use_rep=temp_sedr_key_10p, key_added=cluster_key_seed_10p)
        
        del adata_10p.obsm[temp_sedr_key_10p]
        
        print(f"✓ Clustering complete. Results stored in adata_10p.obs['{cluster_key_seed_10p}']")
    except Exception as e:
        print(f"Error during mclust clustering for seed {seed}: {e}")
        if temp_sedr_key_10p in adata_10p.obsm:
            del adata_10p.obsm[temp_sedr_key_10p]
        print("Clustering failed. Check R setup and mclust installation.")
        continue

    # Calculate ARI
    if 'super_cluster' in adata_10p.obs.columns:
        try:
            true_labels_10p = adata_10p.obs['super_cluster']
            pred_labels_10p = adata_10p.obs[cluster_key_seed_10p]
            ari_10p = metrics.adjusted_rand_score(true_labels_10p, pred_labels_10p)
            all_aris_10p[seed] = ari_10p
            print(f"Adjusted Rand Index (ARI) vs 'super_cluster' for seed {seed}: {ari_10p:.4f}")
        except Exception as e:
            print(f"Error calculating ARI for seed {seed}: {e}")
            all_aris_10p[seed] = np.nan
    else:
        if seed == seeds_to_run[0]:
            print("Ground truth 'super_cluster' column not found in adata_10p.obs. Skipping ARI calculation.")

# --- Combine and Save Results for 10 PCs ---
print("Consolidating results for 10 PCs...")

final_df_10p = pd.DataFrame({
    'barcode': adata_10p.obs.index,
    'x': adata_10p.obsm['spatial'][:, 0],
    'y': adata_10p.obsm['spatial'][:, 1],
    'original_cluster': adata_10p.obs.get('super_cluster', pd.NA),
})

for seed in seeds_to_run:
    cluster_key_seed_10p = f'mclust_seed_{seed}'
    if cluster_key_seed_10p in adata_10p.obs.columns:
        final_df_10p[cluster_key_seed_10p] = adata_10p.obs[cluster_key_seed_10p]
    else:
        final_df_10p[cluster_key_seed_10p] = pd.NA
        print(f"Warning: Clustering results for seed {seed} not found. Filling with NA.")

csv_filename_10p = os.path.join(output_dir_csv_10p, f"sedr_curl_10p_multiseed_k{num_clusters}_results_BestResult.csv")
final_df_10p.to_csv(csv_filename_10p, index=False)
print(f"✓ Saved consolidated results to {csv_filename_10p}")

# Print summary of ARIs for 10 PCs
if all_aris_10p:
    print("--- ARI Summary (10 PCs) ---")
    for seed, ari in all_aris_10p.items():
        print(f"Seed {seed}: {ari:.4f}")
    ari_series_10p = pd.Series(all_aris_10p)
    print(f"Mean ARI: {ari_series_10p.mean():.4f}")
    print(f"Std Dev ARI: {ari_series_10p.std():.4f}")
    print(f"Median ARI: {ari_series_10p.median():.4f}")
    print(f"Min ARI: {ari_series_10p.min():.4f} (Seed: {ari_series_10p.idxmin()})")
    print(f"Max ARI: {ari_series_10p.max():.4f} (Seed: {ari_series_10p.idxmax()})")

# --- Finalize 10 PC version ---
print("=== Multi-seed Script finished (10 PCs) ===")
sys.stdout = sys.__stdout__
sys.stderr = sys.__stderr__

logfile_10p.close()

print(f"Output logged to {log_file_path_10p}") 