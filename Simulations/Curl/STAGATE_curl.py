import sys
import os
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'  # Disable GPU
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'  # Reduce TensorFlow logging

# Set CUDA paths if necessary (adjust if your paths differ)
# os.environ['CUDA_HOME'] = '/sw/eb/sw/CUDA/11.7.0'
# os.environ['PATH'] = f"{os.environ['CUDA_HOME']}/bin:{os.environ['PATH']}"
# os.environ['LD_LIBRARY_PATH'] = f"{os.environ['CUDA_HOME']}/lib64:{os.environ.get('LD_LIBRARY_PATH', '')}"

# Set R paths if necessary (adjust if your paths differ)
# os.environ['R_HOME'] = '/sw/eb/sw/R/4.4.2-gfbf-2024a'
# os.environ['R_LIBS'] = f"/scratch/user/varogovchenko/Rlibs:/sw/eb/sw/R/4.4.2-gfbf-2024a/lib64/R/library"
# os.environ['PATH'] = f"{os.environ['R_HOME']}/bin:{os.environ['PATH']}"

# Open log file in the same directory as the script
logfile = open("STAGATE_Curl_10p_MultiSeed_output.log", "w")
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
input_features = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']
spatial_coords = ['X', 'Y']
rad_cutoff_val = 3 # Set based on previous findings
ground_truth_col = 'super_cluster' # Column name for ground truth clusters in input CSV

seeds_to_run = [
    141, 549, 75, 492, 676, 179, 587, 592, 601, 916, 518, 339, 921, 423, 330,
    388, 273, 286, 61, 807, 283, 127, 165, 952, 311, 597, 473, 594, 605, 101
]

# Define output directories
base_dir = os.path.dirname(os.path.abspath(__file__))
csv_dir = os.path.join(base_dir, 'res_simdata_curl_10p_csv')
# plot_dir = os.path.join(base_dir, 'res_simdata_plot
# s') # No plotting needed

# Create directories if they don't exist
os.makedirs(csv_dir, exist_ok=True)
# os.makedirs(plot_dir, exist_ok=True) # No plotting needed

print(f"\n=== Processing Simulated Data: {os.path.basename(sim_data_path)} ===")
print(f"Running STAGATE for {len(seeds_to_run)} seeds with k={num_clusters}")

# --- Load Data Once ---
print(f"Loading data from: {sim_data_path}")
sim_df = pd.read_csv(sim_data_path)
print(f"Loaded data shape: {sim_df.shape}")
print(f"Available columns: {sim_df.columns.tolist()}")

# Check if ground truth column exists early
has_ground_truth = ground_truth_col in sim_df.columns
if not has_ground_truth:
    print(f"Warning: Ground truth column '{ground_truth_col}' not found in CSV. ARI will not be calculated.")

# --- Create AnnData object --- (We'll reuse this structure)
# Check if required columns exist
required_cols = input_features + spatial_coords
missing_cols = [col for col in required_cols if col not in sim_df.columns]
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

# Extract features and spatial coordinates
features = sim_df[input_features].values
spatial = sim_df[spatial_coords].values

# Create base AnnData
adata = anndata.AnnData(X=features, obs=pd.DataFrame(index=obs_names))
adata.obsm['spatial'] = spatial
# Add original columns to obs for reference (including potential ground truth)
for col in sim_df.columns:
    if col not in input_features and col not in spatial_coords:
        adata.obs[col] = sim_df[col].values

print("Created base AnnData object:")
print(adata)

# --- Preprocessing (Minimal) ---
# Convert adata.X to sparse format once, as train_STAGATE expects it
if not isinstance(adata.X, scipy.sparse.spmatrix):
    print("Converting adata.X to sparse CSR format for train_STAGATE...")
    adata.X = scipy.sparse.csr_matrix(adata.X)

# --- Run STAGATE for multiple seeds ---
ari_scores = []

for i, seed in enumerate(seeds_to_run):
    print(f"\n--- Running Seed {i+1}/{len(seeds_to_run)} ({seed}) ---")

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

    # Create a temporary copy of adata for this run?
    # No, train_STAGATE modifies adata in place, including adding the representation.
    # mclust_R also adds to obs. Let's ensure these are handled.
    # It's safer to recalculate spatial net and train each time.

    print("Calculating spatial network...")
    # Note: Cal_Spatial_Net adds neighbors/distances to adata.uns['Spatial_Net']
    # It should be deterministic given the input coordinates and radius
    # Re-running it ensures the state is clean if it were modified somehow.
    Cal_Spatial_Net(adata, rad_cutoff=rad_cutoff_val)
    # Stats_Spatial_Net(adata) # Optional: prints stats

    print("Training STAGATE model...")
    # train_STAGATE adds representation to adata.obsm['STAGATE']
    # It will overwrite the previous run's representation
    adata_trained = train_STAGATE(adata, alpha=0)

    # Note: Check if train_STAGATE actually uses the seed parameter.
    # If not, the tf/torch/np seeds set above should control its randomness.

    print(f"Running clustering with {num_clusters} clusters...")
    # mclust_R adds 'mclust' to adata.obs
    adata_clustered = mclust_R(adata_trained, used_obsm='STAGATE', num_cluster=num_clusters)

    # Store results in the original DataFrame
    cluster_col_name = f'stagate_seed_{seed}'
    sim_df[cluster_col_name] = adata_clustered.obs['mclust'].values
    print(f"Stored cluster labels in column: {cluster_col_name}")

    # Calculate ARI if ground truth exists
    if has_ground_truth:
        ground_truth_clusters = sim_df[ground_truth_col] # Get from original df
        predicted_clusters = sim_df[cluster_col_name]
        ari = adjusted_rand_score(ground_truth_clusters, predicted_clusters)
        ari_scores.append(ari)
        print(f"Adjusted Rand Index (ARI) for seed {seed}: {ari:.4f}")

# --- Save combined results ---
print("\n--- Saving all results ---")
results_file = os.path.join(csv_dir, f"stagate_simdata_multiseed_{num_clusters}clusters_results.csv")
sim_df.to_csv(results_file, index=False)
print(f"Saved combined clustering results to: {results_file}")

# --- Report ARI Statistics ---
if ari_scores:
    mean_ari = np.mean(ari_scores)
    std_ari = np.std(ari_scores)
    print(f"\n--- ARI Statistics ({len(ari_scores)} runs) ---")
    print(f"Mean ARI: {mean_ari:.4f}")
    print(f"Std Dev ARI: {std_ari:.4f}")
elif has_ground_truth:
     print("\nARI calculation was attempted but failed for all runs.")
else:
     print("\nARI calculation skipped as ground truth column was not found.")

# --- Plotting Section Removed --- 

print(f"\n✓ Finished multi-seed processing for simulated data")

# Close the log file at the end
logfile.close() 