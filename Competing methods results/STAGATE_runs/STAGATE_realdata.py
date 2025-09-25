import sys
import os
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'  # Disable GPU
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'  # Reduce TensorFlow logging

# Set CUDA paths
os.environ['CUDA_HOME'] = '/sw/eb/sw/CUDA/11.7.0'
os.environ['PATH'] = f"{os.environ['CUDA_HOME']}/bin:{os.environ['PATH']}"
os.environ['LD_LIBRARY_PATH'] = f"{os.environ['CUDA_HOME']}/lib64:{os.environ.get('LD_LIBRARY_PATH', '')}"

# Set R paths
os.environ['R_HOME'] = '/sw/eb/sw/R/4.4.2-gfbf-2024a'
os.environ['R_LIBS'] = f""
os.environ['PATH'] = f"{os.environ['R_HOME']}/bin:{os.environ['PATH']}"

# Open log file in the same directory as the script
logfile = open("STAGATE_RealData_output.log", "w")
sys.stdout = logfile
sys.stderr = logfile

import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import torch
import os
import tensorflow as tf

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

from STAGATE.STAGATE.utils import Cal_Spatial_Net, Stats_Spatial_Net, mclust_R
from STAGATE.STAGATE.Train_STAGATE import train_STAGATE

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print(f"Device: {device}")

# Define datasets and cluster numbers to try per dataset
dataset_config = {
    "Brain": [5, 7, 9],
    "Breast": [3, 5, 7],
    "Gut": [3, 5, 7, 9, 12],
    "Gut_reduced": [3, 5, 7, 9, 12],
    "Prostate": [3, 5, 7]  
}

# Define output directories
base_dir = os.path.dirname(os.path.abspath(__file__))
csv_dir = os.path.join(base_dir, 'res_realdata_csv')
plot_dir = os.path.join(base_dir, 'res_realdata_plots')

# Create directories if they don't exist
os.makedirs(csv_dir, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)

base_path = "../ST_Datasets"

for name, cluster_list in dataset_config.items():
    print(f"\n=== Processing {name} ===")
    dataset_path = os.path.join(base_path, name)
    print(dataset_path)

    # Load data
    adata = sc.read_visium(dataset_path, count_file="filtered_feature_bc_matrix.h5", load_images=True)
    adata.var_names_make_unique()
    print(adata)

    # If Gut_reduced: subset spots based on CSV
    if name == "Gut_reduced":
        print("Subsetting Gut_reduced to selected spots...")
        subset_csv = os.path.join(dataset_path, "gut_df_wt_muscle_rownames.csv")
        spot_ids = pd.read_csv(subset_csv, header=None).squeeze().astype(str)
        adata = adata[adata.obs_names.isin(spot_ids)].copy()
        print(f"✓ Subsetted to {adata.n_obs} spots")
        print(adata)

    # Preprocessing
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable].copy()
    print(adata)

    # Run STAGATE
    Cal_Spatial_Net(adata, rad_cutoff=150)
    Stats_Spatial_Net(adata)

    adata = train_STAGATE(adata, alpha=0)

    sc.pp.neighbors(adata, use_rep='STAGATE')
    sc.tl.umap(adata)

    # Clustering at 3, 5, 7 clusters
    for k in cluster_list:
        print(f"Running clustering with {k} clusters...")
        adata = mclust_R(adata, used_obsm='STAGATE', num_cluster=k)
        # Rename the mclust column to include the number of clusters
        adata.obs[f'mclust_{k}'] = adata.obs['mclust'].copy()
        del adata.obs['mclust']  # Remove the original column to avoid confusion

    # Save clustering results
    obs_cols = [f"mclust_{k}" for k in cluster_list]
    df = adata.obs[obs_cols].copy()
    df['barcode'] = adata.obs_names
    df[['x', 'y']] = adata.obsm['spatial']
    
    # Save results with full path
    results_file = os.path.join(csv_dir, f"stagate_{name.lower()}_results.csv")
    df.to_csv(results_file, index=False)
    print(f"Saved clustering results to: {results_file}")

    # Plot clustering results side-by-side
    fig, axs = plt.subplots(1, len(cluster_list), figsize=(5 * len(cluster_list), 5))
    for i, k in enumerate(cluster_list):
        sc.pl.spatial(
            adata,
            color=f"mclust_{k}",
            ax=axs[i],
            show=False,
            title=f"{name} - {k} clusters",
            size=1.5
        )
    plt.tight_layout()
    
    # Save plot with full path
    plot_file = os.path.join(plot_dir, f"stagate_{name.lower()}_multicluster.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved plot to: {plot_file}")

    print(f"✓ Finished {name}")

# Close the log file at the end
logfile.close()
