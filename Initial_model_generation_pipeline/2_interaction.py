# -*- coding: utf-8 -*-
"""
Created on Wed Sep 24 11:29:35 2025

@author: Samantha Crouch
"""

from openai import OpenAI
import re
import csv
import itertools
import os
import pandas as pd

os.chdir("YourFolder")
YOUR_API_KEY = "YourKey"

client = OpenAI(api_key=YOUR_API_KEY, base_url="https://api.perplexity.ai")

######### Query AI for interaction per gene pair
input_file = "extended_gene_pairs.csv"
output_file = "gene_interactions.csv"

with open(input_file, mode="r", newline="") as infile:
    reader = csv.reader(infile)
    next(reader)  

    with open(output_file, mode="w", newline="") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["Gene1", "Gene2", "Interaction"])

        for row in reader:
            gene1, gene2 = row[0], row[1]

            messages = [
                {
                    "role": "system",
                    "content": (
                        "You are a knowledgeable biologist and an expert in senescence. "
                        "Only provide information that is directly supported by multiple peer-reviewed scientific publications. "
                        "Do not speculate or infer missing information."
                    ),
                },
                {
                    "role": "user",
                    "content": (
                        f"For the gene pair ({gene1}, {gene2}), provide only a single number indicating the interaction of {gene1} on {gene2} "
                        "in senescence in fibroblast cells:\n"
                        "-1 for direct inhibition,\n"
                        "-2 for indirect inhibition,\n"
                        "1 for direct activation,\n"
                        "2 for indirect activation,\n"
                        "0 for no known or no reproducible interaction.\n\n"
                        "Only use interactions that are confirmed by at least two independent peer-reviewed publications "
                        "and that apply specifically to the described context. If no clear evidence exists, return 0.\n\n"
                        "Respond with a single number only. No explanations or references."
                    ),
                },
            ]

            response1 = client.chat.completions.create(
                model="sonar",
                messages=messages,
                temperature=0,
                top_p=0.0,
            )

            interaction_matrix = response1.choices[0].message.content
            writer.writerow([gene1, gene2, interaction_matrix])
            print(f"Processed: {gene1}, {gene2}")

import csv
import re

def extract_number(cell: str) -> str:
    """
    Extracts the first (optionally negative) integer from a cell,
    regardless of whether it contains quotes or additional text.
    """
    match = re.search(r"-?\d+", cell)
    return match.group(0) if match else ""

def clean_csv(input_file: str, output_file: str) -> None:
    with open(input_file, "r", encoding="utf-8", newline="") as infile:
        reader = csv.reader(infile)
        cleaned_rows = []

        for row in reader:
            # Join all columns for cleaning steps (bracket/parenthesis removal, etc.)
            full_row = ",".join(row)

            # Remove content inside [] and ()
            full_row = re.sub(r"\[.*?\]", "", full_row)
            full_row = re.sub(r"\(.*?\)", "", full_row)

            # Remove newlines and trim whitespace
            full_row = full_row.replace("\n", " ").replace("\r", " ").strip()

            # Remove markdown bold markers
            full_row = full_row.replace("**", "")

            # Split back into up to 3 columns
            split_row = full_row.split(",", 2)
            if len(split_row) < 3:
                split_row += [""] * (3 - len(split_row))

            # Clean the third column: keep only the first integer (including minus sign)
            split_row[2] = extract_number(split_row[2])

            cleaned_rows.append(split_row[:3])

    with open(output_file, "w", encoding="utf-8", newline="") as outfile:
        writer = csv.writer(outfile)
        writer.writerows(cleaned_rows)

    print(f"Cleaned file saved as {output_file}")

# file paths
input_file = "gene_interactions.csv"
output_file = "cleaned_gene_interactions.csv"
clean_csv(input_file, output_file)
