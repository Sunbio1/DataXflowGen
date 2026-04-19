# -*- coding: utf-8 -*-
"""
Created on Sat May 31 11:14:14 2025

@author: Samantha Crouch
"""

import pandas as pd
import re
import csv
import os

# Set the working directory
os.chdir(r"YourFolder")

# ------------------------------------
# 1. Load local gene list
# ------------------------------------
gene_list_file = (
    "path/loop_gene_list.csv")

gene_list_df = pd.read_csv(gene_list_file)
# assumes genes are in the first column
local_genes = gene_list_df.iloc[:, 0].tolist()

# ------------------------------------
# 2. Load API response as raw text
# ------------------------------------
api_filename = (
    "path/back_loop_api_response.csv" # change path!!
)

with open(api_filename, mode="r", newline="") as file:
    reader = csv.reader(file)
    header = next(reader)  # skip header
    row = next(reader)     # first data row
    api_response = row[0]

# ------------------------------------
# 3. Extract first gene (gene1) from API block
# ------------------------------------
# Remove markdown backticks and split into lines
cleaned = api_response.replace("```", "").strip()
lines = cleaned.splitlines()

# Take the first non-empty line that is not literally "csv"
gene1 = next(
    (
        line.strip()
        for line in lines
        if line.strip() and line.strip().lower() != "csv"
    ),
    None,
)

print(f"API gene: {gene1}")

# ------------------------------------
# 4. Extract full gene list from API block
# ------------------------------------
# Find text between triple backticks (DOTALL so newlines are included)
matches = re.findall(r"```(.*?)```", api_response, re.DOTALL)

if matches:
    extracted_text = matches[0]

    # Build list: drop empty lines and the marker line "csv"
    api_genes = [
        line.strip()
        for line in extracted_text.splitlines()
        if line.strip() and line.strip().lower() != "csv"
    ]

    # Merge with local gene list, remove duplicates, keep alphabetical order
    combined_genes = sorted(set(local_genes).union(api_genes))

    # ------------------------------------
    # 5. Save combined list
    # ------------------------------------
    output_file = (
        "path/loop_gene_list.csv" # change path!!
    )

    pd.DataFrame({"Gene": combined_genes}).to_csv(output_file, index=False)

    print(f"Combined gene list saved to '{output_file}'.")
else:
    print("No text found between triple backticks.")
