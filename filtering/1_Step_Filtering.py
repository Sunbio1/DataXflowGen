# -*- coding: utf-8 -*-
"""
Created on Mon May 26 17:38:04 2025

@author: Samantha Crouch
"""

import pandas as pd
import scanpy as sc
import numpy as np

# 1. Loading h5ad file
adata = sc.read_h5ad('GSM7068366_ETO_day_10.h5ad') # change h5ad file name 
#print (adata.X)
#print(adata.var)

##### 2. Filtering step

# Number of cells and genes before filtering
n_cells_before = adata.n_obs
n_genes_before = adata.n_vars

print(f"Cells before filtering: {n_cells_before}")
print(f"Genes before filtering: {n_genes_before}")

# Perform filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Number of cells and genes after filtering
n_cells_after = adata.n_obs
n_genes_after = adata.n_vars

print(f"Cells after filtering: {n_cells_after} (Difference: {n_cells_before - n_cells_after})")
print(f"Genes after filtering: {n_genes_after} (Difference: {n_genes_before - n_genes_after})")

# === 1. Mark mitochondrial genes (starting with "MT-") ===
adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')

# === 2. Mark ribosomal genes (starting with "RPS" or "RPL") ===
adata.var['rb'] = adata.var_names.str.upper().str.startswith(('RPS', 'RPL'))

# === 3. Calculate QC metrics ===
# (specifying the columns for the QC flags)
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'rb'], inplace=True)

# Violin plots of QC metrics
sc.pl.violin(
    adata,
    ['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'pct_counts_rb'],
    jitter=0.4,
    multi_panel=True
)

# Number of genes starting with "MT-"
mt_genes = adata.var_names.str.upper().str.startswith("MT-")
print(f"Found mitochondrial genes: {mt_genes.sum()} out of {adata.shape[1]}")

# Calculate QC metrics if not already done
sc.pp.calculate_qc_metrics(adata, inplace=True)

# Calculate the mean of total counts per cell
mean_total_counts = adata.obs['total_counts'].mean()

print(f"🔢 Average library size (total_counts per cell): {mean_total_counts:.2f}")

# Normalize each cell to the same library size
sc.pp.normalize_total(adata, target_sum=3e4)

#print(adata.var)

# Extract gene names from the 'gene_name' column
gene_names = adata.var['gene_name'].values  # or .tolist()

# Extract the matrix (Note: Transposed! Rows = genes, columns = cells)
X_matrix = adata.X.T.toarray() if hasattr(adata.X, "toarray") else adata.X.T

# Create a DataFrame with the correct gene names as index
df_counts = pd.DataFrame(
    X_matrix,
    index=gene_names,
    columns=adata.obs_names
)

# Optional: If you want the index (gene_name) as a separate column 'gene_name'
df_counts.reset_index(inplace=True)
df_counts.rename(columns={'index': 'gene_name'}, inplace=True)

# Preview
#print(df_counts.head())

# Save as CSV
df_counts.to_csv("GSM7068366_ETO_day_10_count.csv", index=False) # Change name based on input
