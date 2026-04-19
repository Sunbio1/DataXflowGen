# -*- coding: utf-8 -*-
"""
Created on Sun Dec 14 12:23:21 2025

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
# Build the prompt
# -----------------------------
system_msg = {
    "role": "system",
    "content": (
        "You are a molecular biology expert specialized in human gene regulation and senescence.\n"
        "You MUST be strict about directionality and causality.\n"
        "Prefer regulators that are broadly expressed and commonly measured across bulk/scRNA-seq datasets, including fibroblasts.\n"
        "Avoid lineage-restricted immune-cell transcription factors (e.g., IRF8, SPI1, IKZF1, PAX5, GATA3) unless there is direct fibroblast evidence.\n"
        "If multiple valid candidates exist, do NOT always choose the same one; select a different valid candidate than the most common default.\n"
        "Prioritize evidence from peer-reviewed literature or curated pathway resources; do not guess."
    )
}

sign = "POSITIVE" if reg_norm == "activation" else "NEGATIVE"
verb = "INCREASE" if reg_norm == "activation" else "DECREASE"
forbid = "inhibit" if reg_norm == "activation" else "activate"

direction_block = (
    f"Internally identify one human upstream {sign} regulator of {gene_name}.\n"
    f"It must TEND to {verb} {gene_name} (mRNA, protein, or activity) in human or mammalian cells.\n"
    "Treat as valid evidence any peer-reviewed study where perturbing the regulator "
    "changes the target in the required direction under stress, senescence, proliferation, "
    "differentiation, or oncogenic signalling.\n"
    f"Do NOT return genes that primarily {forbid} {gene_name}.\n"
    "Only output EMPTY if you cannot identify any such regulator.\n"
    "If at least one reasonably well supported regulator exists, you MUST return one of them instead of EMPTY.\n"
)

user_msg = {
    "role": "user",
    "content": (
        f"Target gene: {gene_name} (human HGNC symbol)\n"
        f"Required effect on the target: {reg_norm.upper()}\n\n"
        f"The following genes are already part of the network and must NOT be suggested:\n"
        f"{excluded_genes_str}\n\n"
        f"{direction_block}\n"
        "Hard constraints:\n"
        f"- Must be upstream (not primarily a downstream target of {gene_name}).\n"
        "- Prefer, but do not require, direct causal evidence for the required effect.\n"
        "- Do NOT return pure co-factors without regulatory impact on the target.\n"
        "- Context: human gene regulation with preference for stress/senescence; fibroblast data are helpful but not mandatory.\n"
        "\n"
        "Output format (STRICT):\n"
        "Return ONLY the gene symbol inside triple backticks. No other text.\n"
        "Valid example:\n"
        "```\nHSF1\n```\n"
        "Empty example:\n"
        "```\nEMPTY\n```\n"
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