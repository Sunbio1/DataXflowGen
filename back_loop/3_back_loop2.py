# -*- coding: utf-8 -*-
"""
Created on Sat May 31 10:44:13 2025

@author: Samantha Crouch
"""

from openai import OpenAI
import re
import csv
import itertools
import os
import pandas as pd

os.chdir("YourFolder")

YOUR_API_KEY = 'YourKey'
client = OpenAI(api_key=YOUR_API_KEY, base_url="https://api.perplexity.ai")

######### Query AI for interaction per gene pair
input_file = 'gene_pairs.csv'
output_file = 'gene_interactions.csv'

with open(input_file, mode='r') as infile:
    reader = csv.reader(infile)
    next(reader)

    with open(output_file, mode='w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['Gene1', 'Gene2', 'Interaction'])

        for row in reader:
            gene1, gene2 = row[0], row[1]
            messages = [
                {
                    "role": "system",
                    "content": (
                        "You are a knowledgeable biologist and expert in senescence. "
                        "Only provide information that is directly supported by peer-reviewed scientific publications. "
                        "Do not speculate or infer missing information."
                    ),
                },
                {
                    "role": "user",
                    "content": (
                        f"For the gene pair ({gene1}, {gene2}), provide only a single number indicating the interaction of {gene1} on {gene2} in senescence:\n"
                        "-1 for direct inhibition,\n"
                        "-2 for indirect inhibition,\n"
                        "1 for direct activation,\n"
                        "2 for indirect activation,\n"
                        "0 for no interaction or no known interaction.\n\n"
                        "Only use interactions that are reported in multiple peer-reviewed scientific publications. "
                        "Do not guess or infer. If the interaction is not clearly supported by repeated literature evidence, return 0.\n\n"
                        "Do not include explanations, references, or any additional text. Output only the number."
                    ),
                },
            ]

            response1 = client.chat.completions.create(
                model="sonar",
                messages=messages,
                temperature=0,
                top_p=0
            )

            interaction_matrix = response1.choices[0].message.content
            writer.writerow([gene1, gene2, interaction_matrix])
            print(f"Processed: {gene1}, {gene2}")

import csv
import re

def extract_number(cell):
    """
    Extracts the first (optionally negative) number from the cell,
    regardless of whether it contains quotes or other text.
    """
    match = re.search(r'-?\d+', cell)
    return match.group(0) if match else ''

def clean_csv(input_file, output_file):
    with open(input_file, 'r', encoding='utf-8', newline='') as infile:
        reader = csv.reader(infile)
        cleaned_rows = []

        for row in reader:
            # Merge all columns for processing (only for removing brackets and **)
            full_row = ','.join(row)
            # Remove content inside [] and ()
            full_row = re.sub(r'\[.*?\]', '', full_row)
            full_row = re.sub(r'\(.*?\)', '', full_row)
            # Remove line breaks and trim whitespace
            full_row = full_row.replace('\n', ' ').replace('\r', ' ').strip()
            # Remove **
            full_row = full_row.replace('**', '')
            # Split back into columns
            split_row = full_row.split(',', 2)
            if len(split_row) < 3:
                split_row += [''] * (3 - len(split_row))
            # Clean the third column: keep only the first number (including minus sign if present)
            split_row[2] = extract_number(split_row[2])
            cleaned_rows.append(split_row[:3])

    with open(output_file, 'w', encoding='utf-8', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerows(cleaned_rows)

    print(f"Cleaned file saved as {output_file}")

# Adjust file paths
input_file = 'gene_interactions.csv'
output_file = 'cleaned_gene_interactions.csv'
clean_csv(input_file, output_file)

### END
