# -*- coding: utf-8 -*-
"""
Created on Sat May 31 11:01:02 2025

@author: Sam
"""

import pandas as pd
import os

os.chdir("YourFolder")

# Read the files
adj_matrix_file = "path/loop_modified_adjacency_matrix.csv" # When first loop adj_matrix_file = "path/modified_adjacency_matrix.csv"
new_interactions_file = "cleaned_gene_interactions.csv"

# Read the existing adjacency matrix
adj_matrix = pd.read_csv(adj_matrix_file, index_col=0)

# Read the new interactions
new_interactions = pd.read_csv(new_interactions_file)

# Check the column names
print("Column names in the new interactions file:", new_interactions.columns)

# Assume the third column is the interaction value
value_column = new_interactions.columns[2]

# Extract genes from the new interactions
new_genes = set(new_interactions['Gene1']).union(set(new_interactions['Gene2']))

# Combine all genes from both files
all_genes = sorted(set(adj_matrix.index).union(new_genes))

# Initialize a new matrix with all genes, filling with 0
new_matrix = pd.DataFrame(0, index=all_genes, columns=all_genes)

# Copy existing values into the new matrix
for i in adj_matrix.index:
    for j in adj_matrix.columns:
        new_matrix.loc[i, j] = adj_matrix.loc[i, j]

# Add the new interactions to the matrix
for _, row in new_interactions.iterrows():
    gene1, gene2, value = row['Gene1'], row['Gene2'], row[value_column]
    new_matrix.loc[gene1, gene2] = value

# Save the extended adjacency matrix
new_matrix.to_csv("extended_adjacency_matrix.csv")
print("Extended adjacency matrix has been saved as 'extended_adjacency_matrix.csv'.")
