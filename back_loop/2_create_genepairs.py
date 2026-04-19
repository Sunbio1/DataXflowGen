# -*- coding: utf-8 -*-
"""
Created on Tue Dec 23 17:49:59 2025
@author: Sam
"""
import os
import csv
import re

# Set the working directory
os.chdir(r"YourFolder")

# Step 1: Read the API gene from the CSV file
api_filename = "back_loop_api_response.csv"
with open(api_filename, mode="r", newline="") as file:
    reader = csv.reader(file)
    header = next(reader)  # Skip header
    row = next(reader)     # Read the first row
    api_response = row[0]

    # Prepare the string
    cleaned = api_response.replace("```", "").strip()  # Removes all ```
    lines = cleaned.splitlines()  # Split into lines

# Debug:
print(f"Lines: {lines}")

# Filter out anything that is "csv", ignore empty lines
gene1 = next((line.strip() for line in lines if line.strip().lower() != "csv" and line.strip()), None)

print(f"API gene: {gene1}")

# Step 2: Read the list of genes (from loop_gene_list.csv)
gene_list = []
gene_list_filename = r'path/loop_gene_list.csv'  # change here!
with open(gene_list_filename, mode='r', newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        gene_list.append(row['Gene'])

print("Genes from list:", gene_list)

# Step 3: Create gene pairs (API gene with each gene in the list, plus API gene with itself)
gene_pairs = []

for g in gene_list:
    if g != gene1:
        gene_pairs.append((gene1, g))
        gene_pairs.append((g, gene1))

# Add API gene with itself as a pair
gene_pairs.append((gene1, gene1))

# Step 4: Write the gene pairs to a CSV file
output_filename = "gene_pairs.csv"
with open(output_filename, mode="w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["GeneA", "GeneB"])  # Header
    writer.writerows(gene_pairs)

print(f"CSV file '{output_filename}' was successfully created.")


#### END