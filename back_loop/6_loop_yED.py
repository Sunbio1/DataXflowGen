# -*- coding: utf-8 -*-
"""
Created on Mon May 12 15:56:29 2025

@author: Samantha Crouch
"""

import pandas as pd
import os

# Set the directory
os.chdir("YourFolder")


# 1. Load the adjacency matrix from the CSV file
input_file = "loop_modified_adjacency_matrix.csv"  # Adjust the filename if needed
print(input_file)

# df = pd.read_csv(input_file, index_col=0, sep=';')  # If your CSV uses semicolons

df = pd.read_csv(input_file, index_col=0)  # Use the first column as index

# 2. Transform the matrix into (Gene1, Value, Gene2)
interactions = []

# Iterate through the matrix
for gene1 in df.index:
    for gene2 in df.columns:
        value = df.at[gene1, gene2]
        if value != 0:  # Store only relevant interactions
            interactions.append([gene1, value, gene2])

# 3. Convert to DataFrame
interaction_df = pd.DataFrame(interactions, columns=["Gene1", "Value", "Gene2"])

# 4. Replace numeric values in "Value" with + / - 
interaction_df["Value"] = interaction_df["Value"].apply(lambda x: "+" if x > 0 else "-")

# 5. Create column 5 (unique genes from Gene1 & Gene2)
col5_list = sorted(set(interaction_df["Gene1"]).union(set(interaction_df["Gene2"])))  # remove duplicates & sort
col5_df = pd.DataFrame({"Column5": col5_list})

# 6. Add an empty column
interaction_df["EmptyColumn"] = None  # empty column (for later use in Excel)

# 7. Insert Column5 into the table
max_rows = max(len(interaction_df), len(col5_df))
interaction_df = interaction_df.reindex(range(max_rows))  # pad the DataFrame
interaction_df["Column5"] = None  # placeholder
interaction_df.loc[: len(col5_df) - 1, "Column5"] = col5_df["Column5"].values

# 8. Reorder columns (put EmptyColumn in 4th position)
interaction_df = interaction_df[["Gene1", "Value", "Gene2", "EmptyColumn", "Column5"]]

print(interaction_df)

# 9. Save as an Excel file
output_file = "loop_topo.xlsx"
interaction_df.to_excel(output_file, index=False, header=False, engine="openpyxl")

print(f"✅ File successfully saved as {output_file}!")

###### END
