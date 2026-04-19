# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 12:23:21 2025

@author: Samantha Crouch
"""

import os
import pandas as pd
import re
from openai import OpenAI

# Set working directory
os.chdir(r"YourFolder")

YOUR_API_KEY = "YourKey"
client = OpenAI(api_key=YOUR_API_KEY, base_url="https://api.perplexity.ai")

# -----------------------------
# Load input files
# -----------------------------
gene_input_file = r"path/gene_with_regulation.csv"
excluded_gene_file = r"path/loop_gene_list.csv"
# excluded_gene_file = r"path/Gene_list.csv"  # if starting a new loop

# Find the target gene with the highest Chi²
df = pd.read_csv(gene_input_file)

# Robustly normalize Regulation (Activation/ACtivation/activation/whitespace all treated the same)
df["Regulation_norm"] = df["Regulation"].astype(str).str.strip().str.casefold()

# Keep only valid regulation values
df = df[df["Regulation_norm"].isin(["activation", "inhibition"])].copy()
if df.empty:
    raise ValueError("No valid rows found: Regulation must be Activation/Inhibition (after normalization).")

top_gene_row = df.sort_values(by="Chi2", ascending=False).iloc[0]

gene_name = str(top_gene_row["Variable"]).strip()
chi2_val = float(top_gene_row["Chi2"])
reg_norm = str(top_gene_row["Regulation_norm"]).strip()          # "activation" / "inhibition"
reg_type = reg_norm.capitalize()                                 # "Activation" / "Inhibition" (clean for logging)

# Load excluded genes (robust: drop empty lines/NaN)
excluded_df = pd.read_csv(excluded_gene_file, header=None)
excluded_genes = (
    excluded_df.iloc[:, 0]
    .astype(str)
    .str.strip()
)
excluded_genes = excluded_genes[(excluded_genes != "") & (excluded_genes.str.lower() != "nan")]
excluded_genes_str = ", ".join(excluded_genes.tolist())

# -----------------------------
#  unrestricted prompt
# -----------------------------
system_msg = {
    "role": "system",
    "content": (
        "You are a molecular biology expert specialized in gene regulation and senescence. "
        "Use only reliable sources such as STRING, BioGRID, or recent peer-reviewed scientific literature."
        "Be strict about directionality. Do not guess."
    )
}

user_msg = {
    "role": "user",
    "content": (
        f"I want to identify an **upstream gene** that causes this {reg_type.lower()} effect on the target gene {gene_name}.\n\n"
        f"The following genes are already part of the network and must not be suggested again:\n"
        f"{excluded_genes_str}\n\n"
        f"Find **exactly 1 human gene (GeneCards Symbol)** that is a **documented upstream {reg_role}** of {gene_name}. "
        f"An upstream {reg_role} is defined as a gene or gene product that directly or indirectly 
		causes {reg_direction} regulation "
        f"of the target gene via transcriptional control, signaling, or post-translational modification.\n\n"
        f"Do not return genes that are downstream targets, co-activators, or non-causal interactors of {gene_name}.\n"
        f"Only consider causal, directionally correct, upstream {reg_role}s based on reliable evidence in the
		context of senescence in fibroblast cells.\n\n"
        "If no such gene is clearly documented, return an empty CSV output.\n\n"
        "Return only the gene symbol (HGNC symbol), enclosed in triple backticks like this:\n```csv\nGENE_SYMBOL\n```.\n"
        "Do not include any extra text, comments, headers, formatting, or backslashes. Only the gene symbol, 
		exactly as shown.\n\n"
        "After the triple backticks, add a short explanation (2–4 sentences) explaining why this gene was chosen."
    )
}

# Save prompt locally
with open(f"prompt_for_{gene_name}.txt", "w", encoding="utf-8") as f:
    f.write(user_msg["content"])
print(f"Prompt saved to 'prompt_for_{gene_name}.txt'")

messages = [system_msg, user_msg]

# -----------------------------
# API call
# -----------------------------
response = client.chat.completions.create(
    model="sonar",
    messages=messages,
    temperature=0,
    top_p=0
)

raw_response = (response.choices[0].message.content or "").strip()

# -----------------------------
# Extraction: code block ``` ... ``` (with or without language tag)
# -----------------------------
# Take the first triple-backtick block (either ``` or ```csv etc.)
m = re.search(r"```(?:\w+)?\s*([\s\S]*?)\s*```", raw_response)
gene_symbol = m.group(1).strip() if m else ""

# Optional: take only the first line (in case the model outputs multiple lines)
if gene_symbol:
    gene_symbol = gene_symbol.splitlines()[0].strip()

# Format exactly as desired
formatted_response = f"```\n{gene_symbol}\n```"

print(f"API Response for {gene_name}:\n{formatted_response}")

# -----------------------------
# Save result
# -----------------------------
df_out = pd.DataFrame([{
    "TargetGene": gene_name,
    "Regulation": reg_type,   # clean "Activation"/"Inhibition"
    "Chi2": chi2_val,
    "Response": formatted_response
}])

df_out.to_csv("back_loop_api_response.csv", index=False)
print("API Response saved to 'back_loop_api_response.csv'")

#### END