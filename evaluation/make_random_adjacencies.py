# -*- coding: utf-8 -*-

"""
Created on Sun Dec 14 12:23:21 2025

@author: Samantha Crouch
"""
import os
import random
from collections import Counter

import numpy as np
import pandas as pd


def read_matrix_csv(path: str, sep: str = ",") -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0, sep=sep)
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()
    df = df.applymap(lambda x: int(str(x).strip()))
    if df.shape[0] != df.shape[1]:
        raise ValueError("Matrix is not square.")
    if list(df.index) != list(df.columns):
        raise ValueError("Index/columns are not identical (same order required).")
    return df


def count_zeros(df: pd.DataFrame) -> int:
    return int((df.to_numpy() == 0).sum())


def is_isolated(arr: np.ndarray, k: int) -> bool:
    """Isolated = no edge to other nodes (incoming/outgoing excluding diagonal). Self-loop does NOT count."""
    row_wo = np.delete(arr[k, :], k)
    col_wo = np.delete(arr[:, k], k)
    indeg = int((col_wo != 0).sum())
    outdeg = int((row_wo != 0).sum())
    return indeg == 0 and outdeg == 0


def shuffle_matrix_with_constraints(
    df: pd.DataFrame,
    n_no_input: int,
    n_self_reg: int,
    n_overlap: int,
    seed: int | None = None,
    max_tries: int = 200000,
) -> tuple[pd.DataFrame, dict]:
    rng = random.Random(seed)

    genes = list(df.index)
    n = len(genes)

    if not (0 <= n_no_input <= n and 0 <= n_self_reg <= n and 0 <= n_overlap <= n):
        raise ValueError("Constraints must be within [0..n].")
    if n_overlap > min(n_no_input, n_self_reg):
        raise ValueError("n_overlap must not be larger than min(n_no_input, n_self_reg).")
    if (n_no_input + n_self_reg - n_overlap) > n:
        raise ValueError("Impossible: too many nodes occupied by constraints.")

    arr = df.to_numpy()
    counts = Counter(arr.flatten().tolist())
    zeros_before = count_zeros(df)

    def forced_zero_positions(no_input_nodes: set[int], self_nodes: set[int]) -> set[tuple[int, int]]:
        forced = set()
        # Diagonal: only self_nodes may be non-zero
        for k in range(n):
            if k not in self_nodes:
                forced.add((k, k))
        # No-input: entire column excluding diagonal is 0
        for j in no_input_nodes:
            for i in range(n):
                if i != j:
                    forced.add((i, j))
        return forced

    def build_once():
        all_idx = list(range(n))

        overlap_idx = set(rng.sample(all_idx, n_overlap))
        remaining = [k for k in all_idx if k not in overlap_idx]

        no_only = set(rng.sample(remaining, n_no_input - n_overlap))
        remaining2 = [k for k in remaining if k not in no_only]

        self_only = set(rng.sample(remaining2, n_self_reg - n_overlap))

        no_input_idx = set(overlap_idx) | set(no_only)
        self_idx = set(overlap_idx) | set(self_only)

        fz = forced_zero_positions(no_input_idx, self_idx)

        # enough zeros available?
        if counts.get(0, 0) < len(fz):
            return None

        c = counts.copy()
        out = np.empty((n, n), dtype=int)
        out[:] = 999999

        # forced zeros
        for (i, j) in fz:
            out[i, j] = 0
            c[0] -= 1
            if c[0] < 0:
                return None

        # self-reg diagonals: non-zero values
        nonzero_vals = []
        for v, k in c.items():
            if v != 0 and k > 0:
                nonzero_vals.extend([v] * k)
        if len(nonzero_vals) < n_self_reg:
            return None
        rng.shuffle(nonzero_vals)

        for k in self_idx:
            v = nonzero_vals.pop()
            out[k, k] = v
            c[v] -= 1
            if c[v] < 0:
                return None

        # fill rest
        remaining_vals = []
        for v, k in c.items():
            if k > 0:
                remaining_vals.extend([v] * k)
        rng.shuffle(remaining_vals)

        empties = [(i, j) for i in range(n) for j in range(n) if out[i, j] == 999999]
        if len(remaining_vals) != len(empties):
            return None

        for (i, j), v in zip(empties, remaining_vals):
            out[i, j] = v

        out_df = pd.DataFrame(out, index=genes, columns=genes)

        # validate zeros preserved
        if count_zeros(out_df) != zeros_before:
            return None

        # validate self-reg count
        diag = np.diag(out_df.to_numpy())
        if int((diag != 0).sum()) != n_self_reg:
            return None

        # validate no-input count (excluding diagonal)
        a2 = out_df.to_numpy()
        no_input_found = []
        for j in range(n):
            col_wo = np.delete(a2[:, j], j)
            if int((col_wo != 0).sum()) == 0:
                no_input_found.append(genes[j])
        if len(no_input_found) != n_no_input:
            return None

        # validate overlap count
        self_found = [genes[i] for i in range(n) if a2[i, i] != 0]
        overlap_found = set(no_input_found).intersection(self_found)
        if len(overlap_found) != n_overlap:
            return None

        # validate: no isolated nodes (self-loop does NOT count as a connection)
        for k in range(n):
            if is_isolated(a2, k):
                return None

        meta = {
            "no_input_nodes": sorted(no_input_found),
            "self_reg_nodes": sorted(self_found),
            "overlap_nodes": sorted(overlap_found),
        }
        return out_df, meta

    for _ in range(max_tries):
        res = build_once()
        if res is not None:
            return res

    raise RuntimeError("No valid shuffle solution found (constraints too strict / too many forced zeros).")


def main():
    # --- Paths ---
    workdir = r"YourFolder"
    src_csv = r"path/loop_modified_adjacency_matrix.csv"
    out_dir = os.path.join(workdir, "random_adj")

    os.makedirs(out_dir, exist_ok=True)

    # --- CSV Separator ---
    sep = ","  # if semicolon: ";"

    # --- Settings (fair like original) ---
    # If your original has e.g.: n_no_input=2, n_self_reg=2, overlap=1 (SPI1)
    n_no_input = 2
    n_self_reg = 2
    n_overlap = 1

    # how many randoms
    N = 10  # e.g. 200 / 500 / 1000

    # --- Load Original ---
    df = read_matrix_csv(src_csv, sep=sep)
    zeros_before = count_zeros(df)
    print("Zeros in original:", zeros_before)

    # --- Generate ---
    seen = set()
    base_seed = 1000

    for i in range(1, N + 1):
        seed = base_seed + i
        shuffled, info = shuffle_matrix_with_constraints(
            df,
            n_no_input=n_no_input,
            n_self_reg=n_self_reg,
            n_overlap=n_overlap,
            seed=seed,
            max_tries=200000,
        )

        key = shuffled.to_numpy().tobytes()
        if key in seen:
            # very rare, but if it happens: try a new seed
            seed = base_seed + 100000 + i
            shuffled, info = shuffle_matrix_with_constraints(
                df,
                n_no_input=n_no_input,
                n_self_reg=n_self_reg,
                n_overlap=n_overlap,
                seed=seed,
                max_tries=200000,
            )
            key = shuffled.to_numpy().tobytes()

        seen.add(key)

        out_path = os.path.join(out_dir, f"random_{i:03d}.csv")
        shuffled.to_csv(out_path)

        print(
            f"[{i:03d}] saved: {out_path} | no-input={info['no_input_nodes']} | "
            f"self={info['self_reg_nodes']} | overlap={info['overlap_nodes']}"
        )

    print("\nDONE.")
    print("Output dir:", out_dir)


if __name__ == "__main__":
    main()
