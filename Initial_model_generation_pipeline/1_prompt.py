# -*- coding: utf-8 -*-
"""
Created on Mon May 12 15:39:23 2025

@author: Samantha Crouch
"""

from openai import OpenAI
import re
import csv
import itertools
import os
import pandas as pd

##### API access!
os.chdir('yourFolder')
YOUR_API_KEY = 'yourKey'
client = OpenAI(api_key=YOUR_API_KEY, base_url='https://api.perplexity.ai')

# Read CSV file
gene_list = []
with open('/Gene_list_old.csv', 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        gene_list.append(row['Gene'])

# Create a CSV list in string format
genes_str = ','.join(gene_list)

# First prompt - Find possible genes (nodes)
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
            f"I have the following differentially expressed genes: {genes_str}. "
            "I would like to set a gene interaction network describing senescence in fibroblast cells. "
            "Only include my genes with 3 additional genes that are strongly connected to several of them, "
            "based on multiple peer-reviewed publications. "
            "Only use interactions that are reported in multiple peer-reviewed scientific publications. "
            "Output format: a single line of comma-separated gene names, without duplicates, without any explanation or extra text. "
            "Enclose the result in triple backticks like this: ```"
        ),
    },
]

# Send API request and extract response
response = client.chat.completions.create(
    model="sonar",
    messages=messages,
    temperature=0,  # Deterministic output
    top_p=0         # Restrict to high probability results
)

# Extract the content from the 'choices' attribute
api_response = response.choices[0].message.content
print("API Response:", api_response)

# Put into a list (useful if there are multiple lines)
data = {'response': [api_response]}

# Convert into DataFrame
df = pd.DataFrame(data)

# Save CSV
df.to_csv("api_response.csv", index=False)

print("API Response saved to 'api_response.csv'")

# Flexible regex to extract the genes from the response
gene_pattern = r'```([\s\S]+?)```'
match = re.search(gene_pattern, api_response, re.DOTALL)

if match:
    # Extract CSV content
    csv_genes = match.group(1)
    gene_list = [gene.strip() for gene in csv_genes.split(',')]
    print("Relevant genes:", gene_list)
else:
    print("No relevant genes found.")

# Extract the gene string inside triple backticks
result_string = re.findall(gene_pattern, api_response)

###### Create Gene List
# If genes were extracted, remove spaces and save them as a list
if result_string:
    genes = [gene.strip() for gene in result_string[0].split(',')]
    print("Extracted genes:", genes)
    gene_pairs = list(itertools.combinations(genes, 2))
    print("All possible pairs:", gene_pairs)

    # Save genes to CSV
    with open('Gene_list.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Gene'])
        for gene in genes:
            writer.writerow([gene])

    print("Genes saved to 'Gene_list.csv'.")
else:
    print("No genes found in the 'Results' section.")

import csv
import itertools

# Function to standardize gene names
def standardize_gene_name(gene):
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
    genes = [standardize_gene_name(gene.strip()) for gene in result_string[0].split(',')]
    print("Standardized genes:", genes)
    gene_pairs = list(itertools.combinations(genes, 2))
    print("All possible pairs:", gene_pairs)

    # Save single genes to CSV
    with open('relevant_genes.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Gene'])
        for gene in genes:
            writer.writerow([gene])

    # Save pairs to CSV
    with open('gene_pairs.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Gene1', 'Gene2'])
        for pair in gene_pairs:
            writer.writerow(pair)

    print("Genes and pairs saved.")
else:
    print("No genes found in the 'Results' section.")

###### Extend gene pairs
import csv

def extend_gene_pairs(input_file, output_file):
    with open(input_file, mode='r', newline='', encoding='utf-8') as infile:
        reader = csv.reader(infile)
        next(reader)
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
    with open(output_file, mode='w', newline='', encoding='utf-8') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['Gene1', 'Gene2'])
        writer.writerows(sorted_gene_pairs)

    print(f"Extended and cleaned gene pairs saved to {output_file}.")

input_file = 'gene_pairs.csv'
output_file = 'extended_gene_pairs.csv'
extend_gene_pairs(input_file, output_file)

############################### END
