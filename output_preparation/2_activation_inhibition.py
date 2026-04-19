# -*- coding: utf-8 -*-

"""
Created on Sat May 31 09:45:19 2025

@author: Samantha Crouch
"""
import pandas as pd
import os

# Change working directory (if necessary)

os.chdir("YourFolder")


# File names
main_file = 'f_filtered_data.csv'  # e.g., 'data.csv'
gene_file = 'full_gene_with_highest_chi2.csv'  # e.g., 'gene.csv'
gene_file_with_regulation = 'gene_with_regulation.csv'  # Output file

# Read the main data file
main_df = pd.read_csv(main_file)

# Read the gene file containing the gene with the highest Chi² value
gene_df = pd.read_csv(gene_file)
gene_name = gene_df.iloc[0]['Variable']  # Get the gene name from the first entry

# Create a new 'Gene' column in the main data file (strip '_obs' from Variable column)
main_df['Gene'] = main_df['Variable'].str.replace('_obs', '')

# Filter the main data for the specified gene
gene_data = main_df[main_df['Gene'] == gene_name]

# Calculate mean values
mean_simulated = gene_data['SimulatedData'].mean()
mean_exp = gene_data['ExpData'].mean()

# Calculate the ratio
ratio = mean_simulated / mean_exp

# Determine regulation type
regulation = 'Inhibition' if ratio > 1 else 'Activation'

# Add the regulation type as a new column to the gene file
gene_df['Regulation'] = regulation

# Save the updated gene file with the regulation info
gene_df.to_csv(gene_file_with_regulation, index=False)

# Print results for verification
print(f"Gene: {gene_name}")
print(f"Mean SimulatedData: {mean_simulated}")
print(f"Mean ExpData: {mean_exp}")
print(f"Ratio (Simulated/Exp): {ratio}")
print(f"Regulation: {regulation}")
print(f"Gene file with regulation saved as: {gene_file_with_regulation}")

#### END