# -*- coding: utf-8 -*-
"""
Created on Mon Dec  1 13:48:27 2025

@author: Samantha Crouch
"""


import pandas as pd
import numpy as np

# -------------------------------
# Paths
# -------------------------------
data_path = "path_to_data"
gene_list_path = "path/loop_gene_list.csv"

# -------------------------------
# Files
# -------------------------------
files = [
    ('0', 'GSM7068354_CTRL_2_count_cleaned.csv'),
    ('0.5', 'GSM7068361_ETO_day_0_count_cleaned.csv'),
    ('1', 'GSM7068362_ETO_day_1_count_cleaned.csv'),
    ('2', 'GSM7068363_ETO_day_2_count_cleaned.csv'),
    ('4', 'GSM7068364_ETO_day_4_count_cleaned.csv'),
    ('7', 'GSM7068365_ETO_day_7_count_cleaned.csv'),
    ('10', 'GSM7068366_ETO_day_10_count_cleaned.csv')
]

# -------------------------------
# Read gene list (expects a column named 'Gene')
# -------------------------------
gene_list = pd.read_csv(gene_list_path)
genes = gene_list['Gene'].tolist()

# -------------------------------
# Load CSVs and transpose
# -------------------------------
all_data = {}
for label, file in files:
    filepath = data_path + file
    df = pd.read_csv(filepath, index_col=0)
    df = df.T  # we want one column per gene
    all_data[label] = df

###################################
# SANITY CHECKS / PRESENCE CHECKS #
###################################

####### Collect raw values per gene across all timepoints
gene_values = {gene: [] for gene in genes}
missing_genes_total = set(genes)

for label, df in all_data.items():
    missing_genes_in_file = []
    for gene in genes:
        if gene in df.columns:
            # to numeric, drop NaN/inf
            row = pd.to_numeric(df[gene], errors='coerce').replace([np.inf, -np.inf], np.nan).dropna()
            # allow zeros here because we will use log1p later
            row_nonneg = row[row >= 0]
            if not row_nonneg.empty:
                gene_values[gene].extend(row_nonneg.values)  # keep linear values here
            else:
                print(f"All values for gene {gene} at time {label} are < 0 (or empty) and were removed.")
            missing_genes_total.discard(gene)
        else:
            print(f"Gene {gene} is missing in file with time {label}.")
            missing_genes_in_file.append(gene)

    if missing_genes_in_file:
        print(f"Missing genes in file {label}: {missing_genes_in_file}")

if missing_genes_total:
    print(f"Genes missing in ALL files: {list(missing_genes_total)}")
else:
    print("All genes appear in at least one file.")

####### Per-time presence matrices
exists_rows = []
haspos_rows = []
summary_rows = []

for label, df in all_data.items():
    row_exists = {'t': label}
    row_haspos = {'t': label}
    missing = []

    for gene in genes:
        exists = (gene in df.columns)
        row_exists[gene] = int(exists)

        if exists:
            s = pd.to_numeric(df[gene], errors='coerce').replace([np.inf, -np.inf], np.nan).dropna()
            has_pos = (s > 0).any()
            row_haspos[gene] = int(has_pos)
        else:
            row_haspos[gene] = 0
            missing.append(gene)

    exists_rows.append(row_exists)
    haspos_rows.append(row_haspos)

    present_cnt = sum(v for g, v in row_exists.items() if g != 't')
    pos_cnt     = sum(v for g, v in row_haspos.items() if g != 't')

    summary_rows.append({
        't': label,
        'all_genes_present': (present_cnt == len(genes)),
        'present_count': present_cnt,
        'missing_count': len(genes) - present_cnt,
        'with_pos_values_count': pos_cnt,
        'missing_genes': ','.join(missing) if missing else ''
    })

    if missing:
        print(f"{label}: Missing genes: {missing} | present: {present_cnt}/{len(genes)}, >0: {pos_cnt}")
    else:
        print(f"{label}: All genes present | >0 after cleaning: {pos_cnt}/{len(genes)}")

# Build presence DataFrames and sort columns (t first)
exists_df  = pd.DataFrame(exists_rows).sort_values('t').reset_index(drop=True)
haspos_df  = pd.DataFrame(haspos_rows).sort_values('t').reset_index(drop=True)
summary_df = pd.DataFrame(summary_rows).sort_values('t').reset_index(drop=True)

def sort_cols(df):
    return df[['t'] + [c for c in df.columns if c != 't']]

exists_df  = sort_cols(exists_df)
haspos_df  = sort_cols(haspos_df)
summary_df = sort_cols(summary_df)

# ======================================================================
# MAIN ANALYSIS
# 1) Apply a log transform (log1p by default).
# 2) For each time point and each gene, find the maximum value in log space.
# 3) Divide each log-transformed value by that max -> results are between 0 and 1.
# 4) For each time point and gene, compute the mean of these normalized values.
# 5) For each time point and gene, compute the (inflated) standard error using n_eff.
# ======================================================================

####### Choose log mode
LOG1P = True

####### Effective sample size for the standard error (SE)
N_EFF = 10   # you can adjust this (e.g., 10, 20, 30)


per_time_max_rows = []   # stores per-time per-gene max in LOG space for reference
output_rows = []         # stores mean and SE_eff per-time per-gene

for label, _ in files:  # iterate in defined order
    df = all_data[label]
    max_row = {'t': label}
    out_row = {'t': label}

    for gene in genes:
        if gene not in df.columns:
            # No column -> everything is NaN for this gene at this time
            max_row[f'{gene}_max_log'] = np.nan
            out_row[f'{gene}_obs'] = np.nan
            out_row[f"{gene}_obs_std"] = np.nan
            continue

        # 1---- Clean values
        s = pd.to_numeric(df[gene], errors='coerce').replace([np.inf, -np.inf], np.nan).dropna()

        # 2---- Apply log. With LOG1P we can keep zeros; with plain log we must drop zeros.
        if LOG1P:
            s = s[s > 0]        # allow zeros for log1p
            x = np.log1p(s)      # log(1 + s)
        else:
            s = s[s > 0]         # strictly positive for plain log
            x = np.log(s)        # log(s)

        n = int(x.size)
        if n == 0:
            # nothing left after cleaning
            max_row[f'{gene}_max_log'] = np.nan
            out_row[f'{gene}_obs'] = np.nan
            out_row[f"{gene}_obs_std"] = np.nan
            continue

        # 3---- Per-time, per-gene MAX in log space; then normalize to [0, 1]
        tmax = float(x.max())
        max_row[f'{gene}_max_log'] = tmax

        if np.isfinite(tmax) and tmax > 0:
            normalized = x / tmax      # values now in [0, 1]
            # 4---- Mean (per time & gene)
            mu = normalized.mean()
            # 5---- SD und SE mit effektiver Stichprobengröße N_EFF
            sd = normalized.std(ddof=1) if n > 1 else np.nan
            if np.isfinite(sd):
                se_eff = sd / np.sqrt(N_EFF)
            else:
                se_eff = np.nan
        else:
            mu, se_eff = np.nan, np.nan

        out_row[f'{gene}_obs'] = mu
        out_row[f'{gene}_obs_std'] = se_eff   # Name bleibt *_obs_std für die Tabelle

    per_time_max_rows.append(max_row)
    output_rows.append(out_row)

# -------------------------------
# Build result DataFrames
# -------------------------------
per_time_max_df = pd.DataFrame(per_time_max_rows)
per_time_max_df = per_time_max_df[['t'] + sorted([c for c in per_time_max_df.columns if c != 't'])]

output_df = pd.DataFrame(output_rows)
output_df = output_df[['t'] + sorted([c for c in output_df.columns if c != 't'])]

# -------------------------------
# Final tidy ordering for result_df:
# put 't' first, then per gene: *_obs then *_obs_std
# -------------------------------
result_df = pd.DataFrame(output_rows)

def sort_key(col: str):
    if col == 't':
        return ('', -1)
    if col.endswith('_obs'):
        gene = col[:-4]
        order = 0
    elif col.endswith('_obs_std'):
        gene = col[:-7]
        order = 1
    else:
        gene = col
        order = 99
    return (gene, order)

fixed_cols = ['t']
other_cols = [c for c in result_df.columns if c != 't']
sorted_cols = sorted(other_cols, key=sort_key)
result_df = result_df[fixed_cols + sorted_cols]

# -------------------------------
# Save result
# -------------------------------
# %%
output_file = r"output_path"

result_df.to_csv(output_file, index=False)
print("✅ CSV saved to:", output_file)
