import os
os.environ['R_HOME'] = ''
os.environ['R_LIBS_USER'] = ''

import scanpy as sc
import pandas as pd
from sklearn import metrics
import torch
from sklearn.decomposition import PCA # Import PCA here

import matplotlib.pyplot as plt
import seaborn as sns

import os, sys
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import SEDR

logfile = open("SEDR_RealData_output.log", "w")
sys.stdout = logfile
sys.stderr = logfile

random_seed = 2023
SEDR.fix_seed(random_seed)

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

# Define PC numbers to test
pc_list = [3, 10] 

base_path = "../ST_Datasets"

# Ensure result directories exist
os.makedirs("res_realdata_csv", exist_ok=True)
os.makedirs("res_realdata_plots", exist_ok=True)

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

    # Preprocessing (do only once per dataset)
    adata.layers['count'] = adata.X.toarray() # Keep raw counts for PCA
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable].copy()
    print("Preprocessing completed:")
    print(adata)

    # Constructing neighborhood graph (assuming it doesn't depend on PCA)
    graph_dict = SEDR.graph_construction(adata, 12)
    print("Graph construction completed.")

    # Loop over different numbers of PCs
    for n_pcs in pc_list:
        print(f"--- Processing {name} with n_pcs = {n_pcs} ---")
        
        # PCA calculation
        pca_key = f'X_pca_{n_pcs}'
        print(f"Calculating PCA with {n_pcs} components...")
        # Use raw counts from the layer for PCA
        adata_X = PCA(n_components=n_pcs, random_state=random_seed).fit_transform(adata.layers['count']) 
        adata.obsm[pca_key] = adata_X
        print(f"✓ PCA stored in adata.obsm['{pca_key}']")

        # Training SEDR using the specific PCA result
        sedr_key = f'SEDR_{n_pcs}'
        print(f"Training SEDR using {pca_key}...")
        sedr_net = SEDR.Sedr(adata.obsm[pca_key], graph_dict, mode='clustering', device=device)
        using_dec = True
        if using_dec:
            sedr_net.train_with_dec(N=1)
        else:
            sedr_net.train_without_dec(N=1)
        sedr_feat, _, _, _ = sedr_net.process()
        adata.obsm[sedr_key] = sedr_feat
        print(f"✓ SEDR embedding stored in adata.obsm['{sedr_key}']")

        # Clustering using the specific SEDR embedding
        print(f"Clustering using {sedr_key} for k={cluster_list}...")
        # os.environ['R_USER'] = '' # This should be fine if set once at the top or inherited
        for k in cluster_list:
            cluster_key = f'mclust_pca{n_pcs}_{k}'
            SEDR.mclust_R(adata, k, use_rep=sedr_key, key_added=cluster_key)
        print(f"✓ Clustering completed.")

        # Save clustering results for this PC run
        obs_cols = [f"mclust_pca{n_pcs}_{k}" for k in cluster_list]
        df = adata.obs[obs_cols].copy()
        df['barcode'] = adata.obs_names
        df[['x', 'y']] = adata.obsm['spatial']
        csv_filename = f"res_realdata_csv/sedr_{name.lower()}_pca{n_pcs}_results.csv"
        df.to_csv(csv_filename, index=False)
        print(f"✓ Saved results for n_pcs={n_pcs} to {csv_filename}")

    # Plotting (after loop for both PC runs for the current dataset)
    print(f"Generating comparison plot for {name}...")
    fig, axs = plt.subplots(len(pc_list), len(cluster_list), figsize=(5 * len(cluster_list), 5 * len(pc_list))) 
    
    for pc_idx, n_pcs in enumerate(pc_list):
        for k_idx, k in enumerate(cluster_list):
            cluster_key = f"mclust_pca{n_pcs}_{k}"
            ax = axs[pc_idx, k_idx] if len(pc_list) > 1 else axs[k_idx] # Handle single PC case if needed
            sc.pl.spatial(
                adata,
                color=cluster_key,
                ax=ax,
                show=False,
                title=f"{name} - {k} clus (PCA={n_pcs})",
                size=1.5
            )
            
    plt.tight_layout()
    plot_filename = f"res_realdata_plots/sedr_{name.lower()}_multicluster_pca_comparison.png"
    plt.savefig(plot_filename)
    plt.close(fig) # Close the figure to free memory
    print(f"✓ Saved comparison plot to {plot_filename}")

    print(f"✓ Finished {name}")

logfile.close()
print("\n=== Script finished ===")
