# -*- coding: utf-8 -*-
"""
Created on Sat May 31 08:41:14 2025

@author: Samantha Crouch
"""

######## Calculate Chi² sums
import pandas as pd
import os
os.chdir("YourFolder")

# Path to the input file
input_file = "path/chi2_reduced.csv" 
output_file = "chi2_sums_per_variable_reduced.csv"

# Read CSV
df = pd.read_csv(input_file)

# Calculate Chi² sums per variable
chi2_sums = df.groupby("Variable")["Chi2"].sum().reset_index()
# Remove '_obs' from variable names
chi2_sums['Variable'] = chi2_sums['Variable'].str.replace('_obs', '', regex=False)

# Save the result to a new CSV
chi2_sums.to_csv(output_file, index=False)

print(f"Done! New file saved as '{output_file}'")

######## Find the gene with the highest Chi² value
import pandas as pd

# Path to the input file
input_file = "path/chi2_full.csv" 
output_file = "full_gene_with_highest_chi2.csv"

# Read CSV
df = pd.read_csv(input_file)

# Calculate Chi² sums per variable
chi2_sums = df.groupby("Variable")["Chi2"].sum().reset_index()

# Remove '_obs' from variable names
chi2_sums['Variable'] = chi2_sums['Variable'].str.replace('_obs', '', regex=False)

# Find the gene with the highest Chi² value
top_gene = chi2_sums.sort_values("Chi2", ascending=False).head(1)

# Save to a new CSV
top_gene.to_csv(output_file, index=False)

print(f"Done! The gene with the highest Chi² value was saved as '{output_file}'")

###### Extract filtered data (Simulated vs Exp Data)

# Path to the input file
input1_file = "path/chi2_full.csv" 
output1_file = "data_per_variable_full.csv"

# Read CSV
df = pd.read_csv(input1_file)

# Select only desired columns
filtered_df = df[['Variable', 'ExpData', 'SimulatedData', 'DataIndex']]

# Remove '_obs' from variable names
filtered_df = filtered_df.copy()
filtered_df['Variable'] = filtered_df['Variable'].str.replace('_obs', '', regex=False)

# Save the result as a new CSV
filtered_df.to_csv("f_filtered_data.csv", index=False)

print(f"Done! New file saved as '{output1_file}'")

##### END
