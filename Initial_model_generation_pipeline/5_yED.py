# -*- coding: utf-8 -*-
"""
Created on Mon May 12 15:56:29 2025

@author: Samantha Crouch
"""

import pandas as pd
import os

# Set working directory
os.chdir('YourFolder') 

# ---- 1. Load the adjacency matrix from CSV file ----
input_file = 'modified_adjacency_matrix.csv'  # Adjust the file name if needed
print(input_file)

df = pd.read_csv(input_file, index_col=0)  # Set first column as index

# ---- 2. Transform the matrix into (Gene1, Value, Gene2) format ----
interactions = []

# Iterate through the matrix
for gene1 in df.index:
    for gene2 in df.columns:
        value = df.at[gene1, gene2]
        if value != 0:  # Only save relevant interactions
            interactions.append([gene1, value, gene2])

# ---- 3. Convert to DataFrame ----
interaction_df = pd.DataFrame(interactions, columns=["Gene1", "Value", "Gene2"])

# ---- 4. Replace values in the "Value" column with + / - ----
interaction_df["Value"] = interaction_df["Value"].apply(lambda x: "+" if x > 0 else "-")

# ---- 5. Create column 5 (unique genes from Gene1 & Gene2) ----
column5_list = sorted(set(interaction_df["Gene1"]).union(set(interaction_df["Gene2"])))  # Remove duplicates & sort
column5_df = pd.DataFrame({"Column5": column5_list})  # Convert to DataFrame

# ---- 6. Add an empty column ----
interaction_df["EmptyColumn"] = None  # Empty column (will be filled in Excel later)

# ---- 7. Insert Column5 into the table ----
max_rows = max(len(interaction_df), len(column5_df))  # Determine maximum number of rows
interaction_df = interaction_df.reindex(range(max_rows))  # Fill the DataFrame
interaction_df["Column5"] = None  # Placeholder for Column5
interaction_df.loc[: len(column5_df) - 1, "Column5"] = column5_df["Column5"].values  # Insert values

# ---- 8. Reorder columns (move empty column to 4th position) ----
interaction_df = interaction_df[["Gene1", "Value", "Gene2", "EmptyColumn", "Column5"]]

print(interaction_df)
# ---- 10. Save as Excel file ----
output_file = "topo.xlsx"
interaction_df.to_excel(output_file, index=False, header=False, engine="openpyxl")

print(f"✅ File successfully saved as {output_file}!")

######### END
