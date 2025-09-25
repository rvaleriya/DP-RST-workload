import scanpy as sc
import pandas as pd
import torch
import sys
import os
import matplotlib.pyplot as plt
from GraphST import GraphST
sys.path.append('GraphST/GraphST')
from utils import clustering

logfile = open("GraphST_RealData_output.log", "w")
sys.stdout = logfile
sys.stderr = logfile

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

    # Run GraphST
    model = GraphST.GraphST(adata, device=device)
    adata = model.train()

    # Clustering at 3, 5, 7 clusters
    for k in cluster_list:
        clustering(adata, n_clusters=k, radius=50, method='mclust', refinement=True)
        adata.obs[f"mclust_{k}"] = adata.obs["mclust"]

    # Save clustering results
    obs_cols = [f"mclust_{k}" for k in cluster_list]
    df = adata.obs[obs_cols].copy()
    df['barcode'] = adata.obs_names
    df[['x', 'y']] = adata.obsm['spatial']
    df.to_csv(f"res_realdata_csv/graphst_{name.lower()}_results.csv", index=False)

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
    plt.savefig(f"res_realdata_plots/graphst_{name.lower()}_multicluster.png")
    plt.close()

    print(f"✓ Finished {name}")