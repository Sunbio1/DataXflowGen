# -*- coding: utf-8 -*-
"""
Created on Mon May 12 15:54:21 2025

@author: Samantha Crouch
"""

import numpy as np
import pandas as pd
import os

# Set the working directory
os.chdir("YourFolder")

def create_weighted_adjacency_matrix(data, nodes):
    """
    Create a weighted adjacency matrix for a graph.

    Parameters:
        data (pandas.DataFrame): DataFrame with columns ['Gene1', 'Gene2', 'Interaction'].
        nodes (list): List of unique nodes in the graph.

    Returns:
        numpy.ndarray: Weighted adjacency matrix of the graph.
    """
    # Create a mapping from node to index
    node_to_index = {node: i for i, node in enumerate(nodes)}

    # Initialize a matrix with zeros
    adjacency_matrix = np.zeros((len(nodes), len(nodes)), dtype=int)

    # Populate the adjacency matrix based on the DataFrame
    for _, row in data.iterrows():
        u = node_to_index[row['Gene1']]
        v = node_to_index[row['Gene2']]
        weight = row['Interaction']
        adjacency_matrix[u][v] = weight

    return adjacency_matrix

def save_adjacency_matrix_to_csv(matrix, nodes, output_file):
    """Save the adjacency matrix to a CSV file."""
    df = pd.DataFrame(matrix, index=nodes, columns=nodes)
    df.to_csv(output_file)
    print(f"Adjacency matrix saved to {output_file}")

def display_adjacency_matrix(matrix, nodes):
    """Print the adjacency matrix in a readable format."""
    print("Adjacency Matrix:")
    print(" ", " ".join(nodes))
    for node, row in zip(nodes, matrix):
        print(node, " ".join(map(str, row)))

if __name__ == "__main__":
    csv_file = "cleaned_gene_interactions.csv"
    data = pd.read_csv(csv_file)

    # Check column names and rename if necessary
    if len(data.columns) == 3:
        data.columns = ['Gene1', 'Gene2', 'Interaction']
    else:
        print("Warning: Unexpected column structure:", data.columns)

    nodes = sorted(set(data['Gene1']).union(data['Gene2']))

    adjacency_matrix = create_weighted_adjacency_matrix(data, nodes)
    output_file = "adjacency_matrix.csv"
    save_adjacency_matrix_to_csv(adjacency_matrix, nodes, output_file)
    display_adjacency_matrix(adjacency_matrix, nodes)
