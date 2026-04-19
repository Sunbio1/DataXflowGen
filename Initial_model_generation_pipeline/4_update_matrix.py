# -*- coding: utf-8 -*-
"""
Created on Sun Jan  4 15:52:37 2026

@author: Samantha Crouch
"""

import csv
import os
import numpy as np
from collections import deque

# Set the working directory
os.chdir("YourFolder")  


################## CSV helpers ###################

def detect_delimiter(path):
    """Detect whether the file is comma- or semicolon-separated."""
    with open(path, "r", encoding="utf-8", newline="") as f:
        sample = f.read(4096)
    return csv.Sniffer().sniff(sample, delimiters=";,").delimiter


def read_matrix_with_headers(file_path, delimiter):
    """
    Reads a square adjacency matrix CSV with:
      - first row = column headers (first cell may be empty)
      - first column = row headers
      - remaining cells = integers
    Returns: headers(list), row_headers(list), matrix(np.ndarray)
    """
    with open(file_path, newline="", encoding="utf-8") as csvfile:
        reader = csv.reader(csvfile, delimiter=delimiter)
        data = list(reader)

    if len(data) < 2:
        raise ValueError("CSV seems empty or has no data rows.")

    headers = data[0]
    row_headers = [row[0] for row in data[1:]]

    # Numeric part
    numeric = []
    for r, row in enumerate(data[1:], start=2):
        if len(row) < 2:
            raise ValueError(
                f"Row {r} has <2 columns. Delimiter likely wrong. Row content: {row}"
            )
        try:
            numeric.append([int(cell) for cell in row[1:]])
        except ValueError as e:
            raise ValueError(
                f"Non-integer value in row {r}. If you have floats, adjust int() to float(). "
                f"Row content: {row}"
            ) from e

    matrix = np.array(numeric, dtype=np.int64)

    # Sanity check: should be square
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"Matrix is not square: {matrix.shape}. Check the input CSV structure.")

    return headers, row_headers, matrix


def write_matrix_with_headers(file_path, delimiter, headers, row_headers, matrix):
    """Writes adjacency matrix with same header/index layout as input."""
    rows = []
    rows.append(headers)  # header row

    for i, gene in enumerate(row_headers):
        rows.append([gene] + [str(x) for x in matrix[i].tolist()])

    with open(file_path, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile, delimiter=delimiter)
        writer.writerows(rows)


############# Graph logic ##############
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

def bfs_check_path(A, start, end, weight):
    """BFS to check if there is a path from start to end using only edges == weight."""
    n = A.shape[0]
    visited = np.zeros(n, dtype=bool)
    q = deque([start])
    visited[start] = True

    while q:
        node = q.popleft()
        if node == end:
            return True

        # neighbors with the required edge weight
        neighs = np.where(A[node] == weight)[0]
        for neigh in neighs:
            if not visited[neigh]:
                visited[neigh] = True
                q.append(neigh)

    return False


def update_matrix_for_weight(A, edge_weight_to_check, path_weight):
    """
    Replace edge_weight_to_check with 0 if there exists an alternative path from i to j
    using only edges with path_weight.
    """
    U = A.copy()
    rows, cols = U.shape

    for i in range(rows):
        for j in range(cols):
            if U[i, j] == edge_weight_to_check:
                if bfs_check_path(A, i, j, path_weight):
                    U[i, j] = 0
    return U


def process_adjacency_matrix(A):
    """Apply your two update rules."""
    A1 = update_matrix_for_weight(A,  2,  1)
    A2 = update_matrix_for_weight(A1, -2, -1)
    return A2


###################### Main ####################

if __name__ == "__main__":
    file_path = 'adjacency_matrix.csv'

    delim = detect_delimiter(file_path)
    print("Detected delimiter:", delim)

    headers, row_headers, A = read_matrix_with_headers(file_path, delim)
    print("Original matrix shape:", A.shape)

    U = process_adjacency_matrix(A)
    print("Updated matrix shape:", U.shape)

    # Change detection
    changes_idx = np.argwhere(A != U)
    if changes_idx.size:
        print("Matrix was updated. Changes detected:")
        for i, j in changes_idx:
            # headers[0] is typically empty / index label, so columns start at headers[1]
            col_name = headers[j + 1] if (j + 1) < len(headers) else f"col{j}"
            print(f" - {row_headers[i]} -> {col_name}: {A[i, j]} -> {U[i, j]}")
    else:
        print("Matrix unchanged: No updates were made!")

    # Save outputs (same delimiter as input)
    write_matrix_with_headers('modified_adjacency_matrix.csv', delim, headers, row_headers, U)

    # Optional: save matrix-only (no headers)
    with open(#updated_adjacency_matrix.csv', "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=delim)
        w.writerows(U.tolist())

    # Print final matrix (with delimiter)
    print("\nFinal matrix (with headers):")
    with open('modified_adjacency_matrix.csv', "r", encoding="utf-8") as f:
        for line in f:
            print(line.rstrip("\n"))

######## END