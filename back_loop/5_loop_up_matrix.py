# -*- coding: utf-8 -*-

"""
Created on Tue Dec 23 17:49:59 2025

@author: Samantha Crouch
"""

"""
Adjacency matrix cleanup:
- Replace  2 with 0 if there is a path i->j using only edges ==  1
- Replace -2 with 0 if there is a path i->j using only edges == -1
"""

import csv
from pathlib import Path
from collections import deque
import numpy as np


BASE_DIR = Path(r"YourFolder")
INPUT_FILE = "extended_adjacency_matrix.csv"

OUT_WITH_HEADERS = "loop_modified_adjacency_matrix.csv"
OUT_MATRIX_ONLY = "loop_updated_adjacency_matrix.csv"


############### CSV helpers #################

def detect_delimiter(path: Path) -> str:
    """Detect delimiter ';' or ',' with a safe fallback."""
    sample = path.read_text(encoding="utf-8", errors="replace")[:4096]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=";,")
        return dialect.delimiter
    except csv.Error:
        return ";"


def read_matrix_with_headers(path: Path, delimiter: str):
    """
    Reads CSV layout:
      - first row: column headers (first cell may be empty)
      - first col: row headers
      - rest: integers
    Returns: headers(list[str]), row_headers(list[str]), A(np.ndarray int64)
    """
    with path.open("r", newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f, delimiter=delimiter))

    if len(rows) < 2:
        raise ValueError("CSV seems empty or has no data rows.")

    headers = rows[0]
    data_rows = rows[1:]

    row_headers = []
    numeric = []

    for r_idx, row in enumerate(data_rows, start=2):
        if len(row) < 2:
            raise ValueError(
                f"Row {r_idx} has <2 columns. Delimiter wrong? Row: {row}"
            )
        row_headers.append(row[0])

        try:
            numeric.append([int(x) for x in row[1:]])
        except ValueError as e:
            raise ValueError(
                f"Non-integer value in row {r_idx}. Row: {row}"
            ) from e

    A = np.array(numeric, dtype=np.int64)

    # Basic sanity checks
    if A.shape[0] != A.shape[1]:
        raise ValueError(f"Matrix is not square: {A.shape}. Check input file structure.")

    # Header length should match: 1 index header + N column headers
    if len(headers) != A.shape[1] + 1:
        # Not fatal, but warn by raising: most downstream name-mapping expects it.
        raise ValueError(
            f"Header length mismatch: len(headers)={len(headers)} but expected {A.shape[1] + 1} "
            f"(1 index column + {A.shape[1]} data columns)."
        )

    if len(row_headers) != A.shape[0]:
        raise ValueError("Row header count does not match matrix size.")

    return headers, row_headers, A


def write_matrix_with_headers(path: Path, delimiter: str, headers, row_headers, A: np.ndarray):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=delimiter)
        w.writerow(headers)
        for i, rh in enumerate(row_headers):
            w.writerow([rh] + A[i].tolist())


def write_matrix_only(path: Path, delimiter: str, A: np.ndarray):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=delimiter)
        w.writerows(A.tolist())


############### Graph logic #################
# A is an adjacency matrix. A[i, j] is the directed effect of gene i on gene j.
#
# Codes:
#   1  = direct activation
#   2  = indirect activation
#  -1  = direct inhibition
#  -2  = indirect inhibition
#   0  = no / unknown interaction
#
# Goal:
# Remove redundant indirect edges (2 or -2) if the same connection can already be reached
# by chaining ONLY direct edges of the same sign.

def build_adjlist(A: np.ndarray, weight: int):
    """
    Build adjacency lists for edges where A[i,j] == weight.
    Returns: list[list[int]]
    """
    n = A.shape[0]
    adj = [[] for _ in range(n)]
    # np.where gives all edges for given weight
    src, dst = np.where(A == weight)
    for i, j in zip(src.tolist(), dst.tolist()):
        adj[i].append(j)
    return adj


def reachability_from_adjlist(adj):
    """
    Compute reachability matrix R where R[i,j]=True if j reachable from i.
    Uses BFS from each node: O(n*(n+m)) which is usually way faster than BFS per edge.
    """
    n = len(adj)
    R = np.zeros((n, n), dtype=bool)

    for s in range(n):
        visited = R[s]  # view into row
        q = deque([s])
        visited[s] = True

        while q:
            u = q.popleft()
            for v in adj[u]:
                if not visited[v]:
                    visited[v] = True
                    q.append(v)

    return R


def apply_rule(A: np.ndarray, edge_weight_to_check: int, path_weight: int) -> np.ndarray:
    """
    Replace edge_weight_to_check with 0 if there is a path i->j using only edges == path_weight.
    """
    U = A.copy()

    adj = build_adjlist(A, path_weight)
    R = reachability_from_adjlist(adj)

    mask = (U == edge_weight_to_check)
    # For all positions where mask True, check reachability
    idx = np.argwhere(mask)
    for i, j in idx:
        if R[i, j]:
            U[i, j] = 0

    return U


def process_adjacency_matrix(A: np.ndarray) -> np.ndarray:
    A1 = apply_rule(A,  2,  1)
    A2 = apply_rule(A1, -2, -1)
    return A2


##################  Main ############

def main():
    in_path = BASE_DIR / INPUT_FILE
    if not in_path.exists():
        raise FileNotFoundError(f"Input file not found: {in_path}")

    delim = detect_delimiter(in_path)
    print("Detected delimiter:", repr(delim))

    headers, row_headers, A = read_matrix_with_headers(in_path, delim)
    print("Original matrix shape:", A.shape)

    U = process_adjacency_matrix(A)
    print("Updated matrix shape:", U.shape)

    #  Change detection 
    changes = np.argwhere(A != U)
    if changes.size:
        print("Matrix was updated. Changes detected:")
        for i, j in changes:
            col_name = headers[j + 1]
            print(f" - {row_headers[i]} -> {col_name}: {A[i, j]} -> {U[i, j]}")
    else:
        print("Matrix unchanged: No updates were made.")

    # Write outputs
    write_matrix_with_headers(BASE_DIR / OUT_WITH_HEADERS, delim, headers, row_headers, U)
    write_matrix_only(BASE_DIR / OUT_MATRIX_ONLY, delim, U)

    # Print final file content
    print("\nFinal matrix (with headers):")
    with (BASE_DIR / OUT_WITH_HEADERS).open("r", encoding="utf-8") as f:
        for line in f:
            print(line.rstrip("\n"))


if __name__ == "__main__":
    main()

############ END