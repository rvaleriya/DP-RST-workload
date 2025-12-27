import sys
import os
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'  # Disable GPU
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'  # Reduce TensorFlow logging

# Add STAGATE_local_code to Python path
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(script_dir)
stagate_path = os.path.join(project_root, 'STAGATE_local_code')
sys.path.insert(0, stagate_path)

# Set CUDA paths
os.environ['CUDA_HOME'] = '/sw/eb/sw/CUDA/11.7.0'
os.environ['PATH'] = f"{os.environ['CUDA_HOME']}/bin:{os.environ['PATH']}"
os.environ['LD_LIBRARY_PATH'] = f"{os.environ['CUDA_HOME']}/lib64:{os.environ.get('LD_LIBRARY_PATH', '')}"

# Set R paths
os.environ['R_HOME'] = '/sw/eb/sw/R/4.4.2-gfbf-2024a'
os.environ['R_LIBS'] = f"/scratch/user/varogovchenko/Rlibs:/sw/eb/sw/R/4.4.2-gfbf-2024a/lib64/R/library"
os.environ['PATH'] = f"{os.environ['R_HOME']}/bin:{os.environ['PATH']}"

# Open log file in the same directory as the script
logfile = open("STAGATE_Curl_3p_10p_MultiSeed_output.log", "w")
sys.stdout = logfile
sys.stderr = logfile

import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import scanpy as sc
# import matplotlib.pyplot as plt # No plotting needed
import torch
import tensorflow as tf
import anndata
import scipy.sparse
from sklearn.metrics import adjusted_rand_score
import random # For setting random seeds
import gc  # For garbage collection

# Configure TensorFlow to use GPU and disable eager execution
tf.compat.v1.disable_eager_execution()
physical_devices = tf.config.list_physical_devices('GPU')
if physical_devices:
    try:
        for device in physical_devices:
            tf.config.experimental.set_memory_growth(device, True)
        print(f"Found {len(physical_devices)} GPU(s)")
    except RuntimeError as e:
        print(f"Error configuring GPU: {e}")
else:
    print("No GPU devices found")

from STAGATE.utils import Cal_Spatial_Net, Stats_Spatial_Net, mclust_R
from STAGATE.Train_STAGATE import train_STAGATE

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print(f"Device: {device}")

# --- Parameters ---
sim_data_path = "/scratch/user/varogovchenko/BASTION_HPRC/Simulations/Curl/Curl_sim_data.csv"
num_clusters = 10
input_features_3pc = ['PC1', 'PC2', 'PC3']
input_features_10pc = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']
spatial_coords = ['X', 'Y']
rad_cutoff_val = 3 # Set based on previous findings
ground_truth_col = 'super_cluster' # Column name for ground truth clusters in input CSV

seeds_to_run = [
    141, 549, 75, 492, 676, 179, 587, 592, 601, 916, 518, 339, 921, 423, 330,
    388, 273, 286, 61, 807, 283, 127, 165, 952, 311, 597, 473, 594, 605, 101
]

# Define output directories
base_dir = os.path.dirname(os.path.abspath(__file__))
csv_dir = os.path.join(base_dir, 'res_simdata_curl_3p_10p_csv')
# plot_dir = os.path.join(base_dir, 'res_simdata_plots') # No plotting needed

# Create directories if they don't exist
os.makedirs(csv_dir, exist_ok=True)
# os.makedirs(plot_dir, exist_ok=True) # No plotting needed

print(f"\n=== Processing Simulated Data: {os.path.basename(sim_data_path)} ===")
print(f"Running STAGATE for {len(seeds_to_run)} seeds with k={num_clusters}")
print(f"Will run for both 3 PCs and 10 PCs")

# --- Load Data Once ---
print(f"Loading data from: {sim_data_path}")
sim_df = pd.read_csv(sim_data_path)
print(f"Loaded data shape: {sim_df.shape}")
print(f"Available columns: {sim_df.columns.tolist()}")

# Check if ground truth column exists early
has_ground_truth = ground_truth_col in sim_df.columns
if not has_ground_truth:
    print(f"Warning: Ground truth column '{ground_truth_col}' not found in CSV. ARI will not be calculated.")

# --- Create AnnData objects for both 3PC and 10PC ---
# Check if required columns exist
required_cols_10pc = input_features_10pc + spatial_coords
missing_cols = [col for col in required_cols_10pc if col not in sim_df.columns]
if missing_cols:
    raise ValueError(f"Missing required columns in the CSV: {missing_cols}")

# Use a unique identifier for obs_names if available, otherwise use index
if 'barcode' in sim_df.columns:
    obs_names = sim_df['barcode'].astype(str)
elif 'cell_id' in sim_df.columns:
     obs_names = sim_df['cell_id'].astype(str)
else:
    print("Warning: No 'barcode' or 'cell_id' column found. Using DataFrame index as observation names.")
    obs_names = sim_df.index.astype(str)

if not obs_names.is_unique:
        print(f"Warning: Observation names are not unique. Making them unique.")
        obs_names = anndata.utils.make_index_unique(obs_names)

# Extract spatial coordinates (same for both)
spatial = sim_df[spatial_coords].values

# Create base AnnData objects for 3PC and 10PC
print("\n--- Creating AnnData objects ---")
features_3pc = sim_df[input_features_3pc].values
features_10pc = sim_df[input_features_10pc].values

adata_3pc_base = anndata.AnnData(X=features_3pc, obs=pd.DataFrame(index=obs_names))
adata_3pc_base.obsm['spatial'] = spatial
adata_10pc_base = anndata.AnnData(X=features_10pc, obs=pd.DataFrame(index=obs_names))
adata_10pc_base.obsm['spatial'] = spatial

# Add original columns to obs for reference (including potential ground truth)
for col in sim_df.columns:
    if col not in input_features_10pc and col not in spatial_coords:
        adata_3pc_base.obs[col] = sim_df[col].values
        adata_10pc_base.obs[col] = sim_df[col].values

print("Created AnnData object for 3 PCs:")
print(adata_3pc_base)
print("\nCreated AnnData object for 10 PCs:")
print(adata_10pc_base)

# --- Preprocessing (Minimal) ---
# Convert to sparse format
if not isinstance(adata_3pc_base.X, scipy.sparse.spmatrix):
    print("\nConverting 3PC adata.X to sparse CSR format...")
    adata_3pc_base.X = scipy.sparse.csr_matrix(adata_3pc_base.X)

if not isinstance(adata_10pc_base.X, scipy.sparse.spmatrix):
    print("Converting 10PC adata.X to sparse CSR format...")
    adata_10pc_base.X = scipy.sparse.csr_matrix(adata_10pc_base.X)

# --- Run STAGATE for multiple seeds ---
ari_scores_3pc = []
ari_scores_10pc = []

for i, seed in enumerate(seeds_to_run):
    print(f"\n{'='*60}")
    print(f"--- Running Seed {i+1}/{len(seeds_to_run)} ({seed}) ---")
    print(f"{'='*60}")

    # Set seeds for reproducibility
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    tf.compat.v1.set_random_seed(seed)
    if device.type == 'cuda':
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed) # if use multi-GPU
        # Ensure deterministic operations for CuDNN
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

    # === Run for 3 PCs ===
    print(f"\n--- Processing 3 PCs ---")
    # Create a copy for this run to avoid modifying the base
    adata_3pc = adata_3pc_base.copy()
    
    print("Calculating spatial network...")
    Cal_Spatial_Net(adata_3pc, rad_cutoff=rad_cutoff_val)
    
    print("Training STAGATE model...")
    adata_3pc_trained = train_STAGATE(adata_3pc, alpha=0)
    
    print(f"Running clustering with {num_clusters} clusters...")
    adata_3pc_clustered = mclust_R(adata_3pc_trained, used_obsm='STAGATE', num_cluster=num_clusters)
    
    # Store results
    cluster_col_name_3pc = f'stagate_3pc_seed_{seed}'
    sim_df[cluster_col_name_3pc] = adata_3pc_clustered.obs['mclust'].values
    print(f"Stored 3PC cluster labels in column: {cluster_col_name_3pc}")
    
    # Calculate ARI for 3PC
    if has_ground_truth:
        ground_truth_clusters = sim_df[ground_truth_col]
        predicted_clusters_3pc = sim_df[cluster_col_name_3pc]
        ari_3pc = adjusted_rand_score(ground_truth_clusters, predicted_clusters_3pc)
        ari_scores_3pc.append(ari_3pc)
        print(f"ARI (3PC) for seed {seed}: {ari_3pc:.4f}")
    
    # Clean up to free memory
    del adata_3pc, adata_3pc_trained, adata_3pc_clustered
    gc.collect()
    
    # === Run for 10 PCs ===
    print(f"\n--- Processing 10 PCs ---")
    # Create a copy for this run to avoid modifying the base
    adata_10pc = adata_10pc_base.copy()
    
    print("Calculating spatial network...")
    Cal_Spatial_Net(adata_10pc, rad_cutoff=rad_cutoff_val)
    
    print("Training STAGATE model...")
    adata_10pc_trained = train_STAGATE(adata_10pc, alpha=0)
    
    print(f"Running clustering with {num_clusters} clusters...")
    adata_10pc_clustered = mclust_R(adata_10pc_trained, used_obsm='STAGATE', num_cluster=num_clusters)
    
    # Store results
    cluster_col_name_10pc = f'stagate_10pc_seed_{seed}'
    sim_df[cluster_col_name_10pc] = adata_10pc_clustered.obs['mclust'].values
    print(f"Stored 10PC cluster labels in column: {cluster_col_name_10pc}")
    
    # Calculate ARI for 10PC
    if has_ground_truth:
        predicted_clusters_10pc = sim_df[cluster_col_name_10pc]
        ari_10pc = adjusted_rand_score(ground_truth_clusters, predicted_clusters_10pc)
        ari_scores_10pc.append(ari_10pc)
        print(f"ARI (10PC) for seed {seed}: {ari_10pc:.4f}")
    
    # Clean up to free memory
    del adata_10pc, adata_10pc_trained, adata_10pc_clustered
    gc.collect()

# --- Save combined results ---
print("\n--- Saving all results ---")
results_file = os.path.join(csv_dir, f"stagate_simdata_multiseed_{num_clusters}clusters_results.csv")
sim_df.to_csv(results_file, index=False)
print(f"Saved combined clustering results to: {results_file}")

# --- Report ARI Statistics ---
print(f"\n{'='*60}")
print("--- ARI Statistics ---")
print(f"{'='*60}")

if ari_scores_3pc:
    mean_ari_3pc = np.mean(ari_scores_3pc)
    std_ari_3pc = np.std(ari_scores_3pc)
    print(f"\n3 PCs ({len(ari_scores_3pc)} runs):")
    print(f"  Mean ARI: {mean_ari_3pc:.4f}")
    print(f"  Std Dev ARI: {std_ari_3pc:.4f}")
else:
    if has_ground_truth:
        print("\n3 PCs: ARI calculation was attempted but failed for all runs.")
    else:
        print("\n3 PCs: ARI calculation skipped as ground truth column was not found.")

if ari_scores_10pc:
    mean_ari_10pc = np.mean(ari_scores_10pc)
    std_ari_10pc = np.std(ari_scores_10pc)
    print(f"\n10 PCs ({len(ari_scores_10pc)} runs):")
    print(f"  Mean ARI: {mean_ari_10pc:.4f}")
    print(f"  Std Dev ARI: {std_ari_10pc:.4f}")
else:
    if has_ground_truth:
        print("\n10 PCs: ARI calculation was attempted but failed for all runs.")
    else:
        print("\n10 PCs: ARI calculation skipped as ground truth column was not found.")

# --- Plotting Section Removed --- 

print(f"\n✓ Finished multi-seed processing for simulated data (3 PCs and 10 PCs)")

# Close the log file at the end
logfile.close() 