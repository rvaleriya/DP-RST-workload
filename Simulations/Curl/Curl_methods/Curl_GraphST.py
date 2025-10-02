print("STEP 0: Starting")
import scanpy as sc
import pandas as pd
import torch
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from pathlib import Path
import anndata
import random
from sklearn import metrics
import logging
print("STEP 1: All standard imports successful")

print("STEP 2: Importing GraphST")
from GraphST import GraphST
print("STEP 3: Imported GraphST")
sys.path.append('GraphST/GraphST')
print("STEP 4: Importing clustering")
from utils import clustering
print("STEP 5: All imports done")
print("STEP 1: Imports successful")

# --- Setup ---
print("Starting setup...")
log_file_path = "GraphST_Curl_10p_MultiSeed_output.log"
output_dir_csv = "res_curl_10p_csv"
seeds_to_run = [
    141, 549, 75, 492, 676, 179, 587, 592, 601, 916, 518, 339, 921, 423, 330,
    388, 273, 286, 61, 807, 283, 127, 165, 952, 311, 597, 473, 594, 605, 101
]
num_clusters = 10
graph_radius = 50

print("Creating output directory...")
os.makedirs(output_dir_csv, exist_ok=True)

print("Setting up device...")
device = torch.device('cpu')
print("Device set to CPU")

# Set up logging
logging.basicConfig(
    filename=log_file_path,
    level=logging.INFO,
    format='%(asctime)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
console = logging.StreamHandler()
console.setLevel(logging.INFO)
logging.getLogger('').addHandler(console)
print(f"Logging configured to: {log_file_path}")

# --- Load Data ---
data_file = "/scratch/user/varogovchenko/BASTION_HPRC/Simulations/Curl/Curl_sim_data.csv"
logging.info(f"Loading data from: {data_file}")
df = pd.read_csv(data_file)

print(df.head())

logging.info("Data loaded")

# --- Prepare AnnData ---
logging.info("Creating AnnData object...")
features = df[['PC1', 'PC2', 'PC3','PC4', 'PC5', 'PC6','PC7', 'PC8', 'PC9','PC10']].values
spatial_coords = df[['X', 'Y']].values
observations = df[['super_cluster']].copy()

# Create AnnData with preprocessed data
adata_base = anndata.AnnData(X=features.copy())
adata_base.obsm['spatial'] = spatial_coords.copy()
adata_base.obs = observations.copy()
adata_base.obs.index = df.index.astype(str)
adata_base.var_names = ['Y1', 'Y2', 'Y3', 'Y4', 'Y5', 'Y6', 'Y7', 'Y8', 'Y9', 'Y10']
adata_base.var_names_make_unique()

logging.info("AnnData object created:")
logging.info(str(adata_base))

# --- Loop Over Seeds ---
all_aris = {}
for seed in seeds_to_run:
    logging.info(f"\n--- Running seed {seed} ---")
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)

    adata = adata_base.copy()
    
    try:
        logging.info("Initializing GraphST model...")
        model = GraphST.GraphST(adata, device=device)
        
        logging.info("Starting model training...")
        try:
            adata = model.train()
            logging.info("Model training completed")
        except Exception as train_e:
            logging.error(f"Training failed: {str(train_e)}")
            raise
            
        logging.info(f"Starting clustering with radius {graph_radius}...")
        clustering(adata, n_clusters=num_clusters, radius=graph_radius, method='mclust', refinement=True)
        
        # Store only refined results
        cluster_key = f"refined_mclust_seed_{seed}"
        adata.obs[cluster_key] = adata.obs["domain"]

        if 'super_cluster' in adata.obs.columns:
            ari = metrics.adjusted_rand_score(adata.obs['super_cluster'], adata.obs[cluster_key])
            all_aris[seed] = ari
            logging.info(f"✓ ARI vs true (radius={graph_radius}): {ari:.4f}")
        else:
            all_aris[seed] = np.nan
            logging.info("No ground truth. Skipping ARI.")

        adata_base.obs[cluster_key] = adata.obs[cluster_key]

    except Exception as e:
        logging.error(f"Failed at seed {seed}: {e}")
        logging.error(f"Error type: {type(e)}")
        logging.error(f"Error details: {str(e)}")
        adata_base.obs[f"refined_mclust_seed_{seed}"] = pd.NA
        all_aris[seed] = np.nan
        continue  # Skip to next seed on error

# --- Save Final Results ---
logging.info("Consolidating and saving results...")
final_df = pd.DataFrame({
    'barcode': adata_base.obs.index,
    'x': adata_base.obsm['spatial'][:, 0],
    'y': adata_base.obsm['spatial'][:, 1],
    'original_cluster': adata_base.obs.get('super_cluster', pd.NA),
})

for seed in seeds_to_run:
    col = f'refined_mclust_seed_{seed}'
    final_df[col] = adata_base.obs.get(col, pd.NA)

output_path = os.path.join(output_dir_csv, f"graphst_curl_10p_multiseed_k{num_clusters}_radius{graph_radius}_results.csv")
final_df.to_csv(output_path, index=False)
logging.info(f"✓ Results saved to {output_path}")

# --- Print ARI Summary ---
if all_aris:
    logging.info("\n--- ARI Summary (radius={}) ---".format(graph_radius))
    ari_series = pd.Series(all_aris)
    for seed, ari in all_aris.items():
        logging.info(f"Seed {seed}: {ari:.4f}")
    logging.info(f"Mean ARI: {ari_series.mean():.4f}")
    logging.info(f"Std Dev ARI: {ari_series.std():.4f}")
    logging.info(f"Median ARI: {ari_series.median():.4f}")
    logging.info(f"Min ARI: {ari_series.min():.4f} (Seed: {ari_series.idxmin()})")
    logging.info(f"Max ARI: {ari_series.max():.4f} (Seed: {ari_series.idxmax()})")

print(f"✓ Finished. Output logged to {log_file_path}")
