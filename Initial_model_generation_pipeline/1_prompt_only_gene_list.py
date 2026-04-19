# -*- coding: utf-8 -*-
"""
Created on Sun Dec 14 12:23:21 2025

@author: Samantha Crouch
"""

from openai import OpenAI
import re
import csv
import itertools
import os
import pandas as pd

os.chdir("YourFolder")
YOUR_API_KEY = "yourKey"
client = OpenAI(api_key=YOUR_API_KEY, base_url="https://api.perplexity.ai")

# Read the CSV file
gene_list = []
with open("path/gene_list_old.csv", "r", newline="") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        gene_list.append(row["Gene"])

# Create a comma-separated gene list as a string
genes_str = ",".join(gene_list)

messages = [
    {
        "role": "system",
        "content": (
            "You are a knowledgeable biologist and medical expert. "
            "Only provide information that is well-established in the scientific literature. "
            "Avoid speculation or assumptions."
        ),
    },
    {
        "role": "user",
        "content": (
            f"I have the following differentially expressed genes (this list is authoritative): {genes_str}\n"
            "Task: Build a senescence-relevant gene interaction network for lung cancer using ONLY genes from my list.\n"
            "CRITICAL CONSTRAINTS:\n"
            "1) Do NOT add any genes that are not in my list.\n"
            "2) Do NOT omit any gene from my list. Every input gene must appear in your output exactly once.\n"
            "3) Interactions should be based only on well-established evidence reported in multiple peer-reviewed publications.\n"
            "   If some genes have no such well-established interactions, they must still be included as isolated nodes.\n"
            "\n"
            "Output format: EXACTLY one line of comma-separated gene symbols containing all input genes, no duplicates, no extra text.\n"
            "Enclose the result in triple backticks like this: ```...```"
        ),
    },
]

# Send API request and extract response
response = client.chat.completions.create(
    model="sonar",
    messages=messages,
    temperature=0,  # Deterministic output
    top_p=0,        # Restrict to high probability results
    max_tokens=2000
)

# Extract the content from the response
api_response = response.choices[0].message.content
print("API Response:", api_response)

# Save response to CSV
df = pd.DataFrame({"response": [api_response]})
df.to_csv("api_response.csv", index=False)
print("API Response saved to 'api_response.csv'")

# Extract genes from triple backticks
gene_pattern = r"```([\s\S]+?)```"
match = re.search(gene_pattern, api_response, re.DOTALL)

if match:
    csv_genes = match.group(1)
    gene_list = [gene.strip() for gene in csv_genes.split(",")]
    print("Relevant genes:", gene_list)
else:
    print("No relevant genes found.")

result_string = re.findall(gene_pattern, api_response)

# If genes were extracted, create pairwise combinations and save
if result_string:
    genes = [gene.strip() for gene in result_string[0].split(",")]
    print("Extracted genes:", genes)

    gene_pairs = list(itertools.combinations(genes, 2))
    print("All possible pairs:", gene_pairs)

    # Save genes to CSV
    with open("Gene_list.csv", mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Gene"])
        for gene in genes:
            writer.writerow([gene])
    print("Genes saved to 'Gene_list.csv'.")
else:
    print("No genes found in the extracted section.")

# Function to standardize gene names
def standardize_gene_name(gene: str) -> str:
    gene = gene.replace("NFκB1", "NFKB1")
    gene = gene.replace("NF-κB1", "NFKB1")
    gene = gene.replace("NF-κB (RELA)", "NFKB1")
    gene = gene.replace("κ", "K")
    gene = gene.replace("β", "b")
    gene = gene.replace("TGF-β", "TGFb")
    gene = gene.replace("Cystatin-C", "CST3")
    gene = gene.replace("C/EBP", "CEBPB")
    gene = gene.replace("Aurora kinase B", "AURKB")
    return gene

if result_string:
    genes = [standardize_gene_name(gene.strip()) for gene in result_string[0].split(",")]
    print("Standardized genes:", genes)

    gene_pairs = list(itertools.combinations(genes, 2))
    print("All possible pairs:", gene_pairs)

    # Save standardized genes
    with open("relevant_genes.csv", mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Gene"])
        for gene in genes:
            writer.writerow([gene])

    # Save pairs
    with open("gene_pairs.csv", mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Gene1", "Gene2"])
        for pair in gene_pairs:
            writer.writerow(pair)

    print("Genes and pairs saved.")
else:
    print("No genes found in the extracted section.")

###### Extend gene pairs
def extend_gene_pairs(input_file: str, output_file: str) -> None:
    with open(input_file, mode="r", newline="", encoding="utf-8") as infile:
        reader = csv.reader(infile)
        next(reader)  # Skip header
        gene_pairs = [row for row in reader]

    unique_gene_pairs = set()
    for row in gene_pairs:
        if len(row) >= 2:
            gene1, gene2 = row[:2]
            unique_gene_pairs.add((gene1, gene2))
            unique_gene_pairs.add((gene2, gene1))
            unique_gene_pairs.add((gene1, gene1))
            unique_gene_pairs.add((gene2, gene2))

    sorted_gene_pairs = sorted(unique_gene_pairs)

    with open(output_file, mode="w", newline="", encoding="utf-8") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["Gene1", "Gene2"])
        writer.writerows(sorted_gene_pairs)

    print(f"Extended and cleaned gene pairs saved to {output_file}.")

input_file = "gene_pairs.csv"
output_file = "extended_gene_pairs.csv"
extend_gene_pairs(input_file, output_file)

############################### END
