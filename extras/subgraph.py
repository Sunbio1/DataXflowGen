# -*- coding: utf-8 -*-

"""
Created on Sun Dec 14 12:23:21 2025

@author: Samantha Crouch
"""

"""
Iterate through all gene subsets and classify each induced subgraph component as:
  1) isolated node
  2) directed simple path
  3) single directed cycle
  4) other (connected but not path/cycle)

Input:
  - modified_adjacency_matrix.csv (header row + header column with gene names)
Output:
  - subgraph.csv (semicolon-separated)
"""

from __future__ import annotations

import csv
import itertools
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Sequence, Set, Tuple

import numpy as np

input_path = Path("path/modified_adjacency_matrix.csv")
output_path = Path("path/subgraph.csv")
# -----------------------------
# I/O
# -----------------------------
def read_csv_matrix(file_path: Path) -> Tuple[List[str], np.ndarray]:
    with file_path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        rows = list(reader)

    if not rows or len(rows) < 2:
        raise ValueError("CSV looks empty or malformed.")

    genes = rows[0][1:]
    matrix = [[int(v) for v in row[1:]] for row in rows[1:]]

    G = np.asarray(matrix, dtype=int)
    if G.shape[0] != G.shape[1]:
        raise ValueError(f"Adjacency matrix must be square, got {G.shape}.")
    if len(genes) != G.shape[0]:
        raise ValueError("Header gene count does not match matrix size.")

    return genes, G



############ Graph utilities ##########################

def build_induced_subgraph(
    G: np.ndarray, nodes: Sequence[int]
) -> Tuple[Dict[int, List[int]], Dict[int, Set[int]]]:
    """
    Build:
      - subG: directed adjacency list (u -> [v,...]) for edges with G[u,v] != 0
      - undirected_subG: undirected adjacency (u - v) if G[u,v] != 0 OR G[v,u] != 0
    """
    subG: DefaultDict[int, List[int]] = defaultdict(list)
    undirected_subG: DefaultDict[int, Set[int]] = defaultdict(set)

    node_set = set(nodes)
    for u in nodes:
        # directed edges
        for v in nodes:
            if u == v:
                continue
            if G[u, v] != 0:
                subG[u].append(v)
            if G[u, v] != 0 or G[v, u] != 0:
                undirected_subG[u].add(v)
                undirected_subG[v].add(u)

        # ensure keys exist even if isolated in one structure
        subG[u] = subG[u]
        undirected_subG[u] = undirected_subG[u]

    # Remove any accidental neighbors not in subset (defensive)
    for u in list(subG.keys()):
        subG[u] = [v for v in subG[u] if v in node_set and v != u]
    for u in list(undirected_subG.keys()):
        undirected_subG[u] = {v for v in undirected_subG[u] if v in node_set and v != u}

    return dict(subG), dict(undirected_subG)


def find_undirected_components(undirected_subG: Dict[int, Set[int]], nodes: Sequence[int]) -> List[Set[int]]:
    visited: Set[int] = set()
    components: List[Set[int]] = []

    for start in nodes:
        if start in visited:
            continue
        stack = [start]
        visited.add(start)
        comp = {start}

        while stack:
            cur = stack.pop()
            for nbr in undirected_subG.get(cur, set()):
                if nbr not in visited:
                    visited.add(nbr)
                    comp.add(nbr)
                    stack.append(nbr)

        components.append(comp)

    return components


def contains_multi_node_cycle(subG: Dict[int, List[int]], nodes: Set[int]) -> bool:
    # 0=unvisited, 1=visiting, 2=done
    color = {n: 0 for n in nodes}

    def dfs(u: int) -> bool:
        color[u] = 1
        for v in subG.get(u, []):
            if v not in nodes or v == u:
                continue
            if color[v] == 0:
                if dfs(v):
                    return True
            elif color[v] == 1:
                return True
        color[u] = 2
        return False

    for n in nodes:
        if color[n] == 0 and dfs(n):
            return True
    return False


def is_directed_simple_path(subG: Dict[int, List[int]], comp: Set[int]) -> bool:
    if len(comp) < 2:
        return False

    in_deg = {n: 0 for n in comp}
    out_deg = {n: 0 for n in comp}

    for u in comp:
        for v in subG.get(u, []):
            if v in comp and v != u:
                out_deg[u] += 1
                in_deg[v] += 1

    starts = [n for n in comp if in_deg[n] == 0]
    ends = [n for n in comp if out_deg[n] == 0]
    if len(starts) != 1 or len(ends) != 1:
        return False

    for n in comp:
        if n not in starts and n not in ends:
            if in_deg[n] != 1 or out_deg[n] != 1:
                return False

    if contains_multi_node_cycle(subG, comp):
        return False

    return True


def order_directed_path(subG: Dict[int, List[int]], comp: Set[int]) -> List[int]:
    in_deg = {n: 0 for n in comp}
    out_deg = {n: 0 for n in comp}

    for u in comp:
        for v in subG.get(u, []):
            if v in comp and v != u:
                out_deg[u] += 1
                in_deg[v] += 1

    starts = [n for n in comp if in_deg[n] == 0]
    if len(starts) != 1:
        return sorted(comp)

    start = starts[0]
    order: List[int] = []
    visited: Set[int] = set()
    cur = start

    while True:
        order.append(cur)
        visited.add(cur)
        nxt = None
        for v in subG.get(cur, []):
            if v in comp and v not in visited and v != cur:
                nxt = v
                break
        if nxt is None:
            break
        cur = nxt

    return order


def is_single_directed_cycle(subG: Dict[int, List[int]], comp: Set[int]) -> bool:
    if len(comp) < 2:
        return False

    in_deg = {n: 0 for n in comp}
    out_deg = {n: 0 for n in comp}

    for u in comp:
        for v in subG.get(u, []):
            if v in comp and v != u:
                out_deg[u] += 1
                in_deg[v] += 1

    # every node must have exactly 1 in and 1 out
    for n in comp:
        if in_deg[n] != 1 or out_deg[n] != 1:
            return False

    # verify it's one single cycle (not multiple disjoint cycles inside comp)
    start = min(comp)
    cur = start
    visited: Set[int] = set()

    for _ in range(len(comp)):
        if cur in visited:
            return False
        visited.add(cur)

        succ = [v for v in subG.get(cur, []) if v in comp and v != cur]
        if len(succ) != 1:
            return False
        cur = succ[0]

    return cur == start and len(visited) == len(comp)


def order_single_cycle(subG: Dict[int, List[int]], comp: Set[int]) -> List[int]:
    start = min(comp)
    order = [start]
    visited = {start}
    cur = start

    while True:
        succ = [v for v in subG.get(cur, []) if v in comp and v != cur]
        if not succ:
            break
        nxt = succ[0]
        if nxt == start:
            break
        if nxt in visited:
            break
        order.append(nxt)
        visited.add(nxt)
        cur = nxt

    return order


def classify_component(subG: Dict[int, List[int]], comp: Set[int]) -> Tuple[int, List[int]]:
    comp_list = sorted(comp)

    if len(comp) == 1:
        return 1, comp_list

    if is_directed_simple_path(subG, comp):
        return 2, order_directed_path(subG, comp)

    if is_single_directed_cycle(subG, comp):
        return 3, order_single_cycle(subG, comp)

    return 4, comp_list


def process_graph(G: np.ndarray, nodes: Sequence[int]) -> Dict[int, List[List[int]]]:
    subG, undirected_subG = build_induced_subgraph(G, nodes)
    components = find_undirected_components(undirected_subG, nodes)

    results: Dict[int, List[List[int]]] = {1: [], 2: [], 3: [], 4: []}
    for comp in components:
        cls_id, ordering = classify_component(subG, comp)
        results[cls_id].append(ordering)

    return results


def translate_results(classification: Dict[int, List[List[int]]], index_to_gene: Sequence[str]) -> Dict[int, List[List[str]]]:
    out: Dict[int, List[List[str]]] = {1: [], 2: [], 3: [], 4: []}
    for cls_id, comps in classification.items():
        out[cls_id] = [[index_to_gene[i] for i in comp] for comp in comps]
    return out


def format_components(components: List[List[str]]) -> str:
    """
    CSV-friendly formatting:
      - each component as "A->B->C"
      - components separated by " | "
    """
    return " | ".join("->".join(comp) for comp in components)


################## Main ##########################

def main() -> None:
    # NOTE: iterating *all* subsets is exponential (2^N - 1). This will explode fast.


    genes, G = read_csv_matrix(input_path)

    N = len(genes)
    all_nodes = list(range(N))
    rows_out: List[List[str]] = []

    for size in range(1, N + 1):
        for subset in itertools.combinations(all_nodes, size):
            nodes = list(subset)
            classification = process_graph(G, nodes)
            named = translate_results(classification, genes)

            subset_names = [genes[i] for i in nodes]
            subset_str = ", ".join(subset_names)

            for cls_id in (1, 2, 3, 4):
                comps = named[cls_id]
                if comps:
                    rows_out.append([subset_str, str(cls_id), format_components(comps)])

    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter=";")
        writer.writerow(["Subset", "Class", "Ordering"])
        writer.writerows(rows_out)


if __name__ == "__main__":
    main()

############## END