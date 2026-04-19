# -*- coding: utf-8 -*-
"""
Created on Tue Dec 23 17:49:59 2025

@author: Samantha Crouch

Upstream root-cause analysis for a single problem gene (from gene_with_regulation.csv),
including a robust 2-hop analysis even if the required side (Activation/Inhibition)
is missing in the network.

Core choices in this version:
- Adjacency weights are treated as SIGN only (presence + sign). No strength assumption.
- Optional state-gating from chi2_full long-format file (0..1):
  Time,Variable,ExpData,SimulatedData,Std,Chi2,ModelIndex,DataIndex
  Snapshot time is chosen as the time where the TARGET has maximal Chi2.
- Balance is shown as counts (#activators/#inhibitors), not sum of weights.
- In prints, show both: w_raw (from file) and w_eff (binary used for scoring).

Score:
  Score(R -> T) = I(w_raw != 0) * chi2_sum(R) * state_snapshot(R)
  (if state_snapshot missing -> defaults to 1.0)
"""

import numpy as np
import csv


ADJ_FILE = r"path/extended_adjacency_matrix.csv"
CHI2_SUM_FILE = r"path/chi2_sums_per_variable_full.csv"
GENE_WITH_REG_FILE = r"path/gene_with_regulation.csv"

# Long-format file containing ExpData/SimulatedData per timepoint (your example)
CHI2_FULL_FILE = r"path/chi2_full.csv"
STATE_SOURCE = "ExpData"   # "ExpData" or "SimulatedData" (must match column names in CHI2_FULL_FILE)

TOP_K = 15
DO_TWO_HOP = True
TWO_HOP_TOP_UPSTREAM = 5      # how many top upstream candidates to expand in 2-hop
TWO_HOP_TOP_PER_NODE = 15     # how many upstream edges per 2-hop node to print

SHOW_SELF_LOOPS = True             # print self-loop lines (if present)
EXCLUDE_SELF_FROM_2HOP_RANK = True # keep self-loop, but exclude it from ranked upstream in 2-hop


################## HELPERS ##################
def sniff_dialect(f, sample_bytes=4096):
    pos = f.tell()
    sample = f.read(sample_bytes)
    f.seek(pos)
    try:
        return csv.Sniffer().sniff(sample, delimiters=',;')
    except csv.Error:
        return csv.excel


def is_number(s: str) -> bool:
    try:
        float(s.strip().replace(',', '.'))
        return True
    except:
        return False


def norm_gene(name: str) -> str:
    return (name or "").replace("_obs", "").strip()


def is_need_activation(need: str) -> bool:
    s = (need or "").lower()
    return ("activ" in s) or s.startswith("act") or (s == "a")


def is_need_inhibition(need: str) -> bool:
    s = (need or "").lower()
    return ("inhib" in s) or s.startswith("inh") or (s == "i")


################# Input/Output ###############
def read_adjacency(path):
    """
    Expected format:
      Header: [empty], G1, G2, ...
      Rows:   Gi, 0, 1, -1, ...
    Assumes row order matches header order (common for adjacency exports).
    """
    with open(path, 'r', newline='', encoding='utf-8-sig') as f:
        dialect = sniff_dialect(f)
        reader = csv.reader(f, dialect)
        rows = [[c.strip() for c in row] for row in reader if row and any(c.strip() for c in row)]

    if not rows or len(rows[0]) < 2:
        raise RuntimeError("Adjacency file looks empty or has no valid header format.")

    genes = [norm_gene(g) for g in rows[0][1:]]  # normalize header genes too
    data = rows[1:]

    M = np.zeros((len(data), len(genes)), dtype=int)
    for i, r in enumerate(data):
        if len(r) < 2:
            continue
        vals = r[1:]
        try:
            M[i, :] = np.array([int(x) for x in vals], dtype=int)
        except ValueError:
            M[i, :] = np.array([int(float(x.replace(',', '.'))) for x in vals], dtype=int)

    return genes, M


def read_chi2_sums(path, valid_genes):
    """
    Reads chi² per gene from chi2_sums_per_variable_full.csv.
    Robust:
    - detects header vs no header
    - tries common column names
    - else: name column = 0, chi² column = first numeric column != name
    """
    chi2 = {}
    with open(path, 'r', newline='', encoding='utf-8-sig') as f:
        dialect = sniff_dialect(f)
        reader = csv.reader(f, dialect)
        rows = [row for row in reader if row and any(c.strip() for c in row)]
        if not rows:
            return chi2

        header = [c.strip() for c in rows[0]]
        data_start = 1
        if any(is_number(c) for c in header):
            header = None
            data_start = 0

        name_idx = 0
        chi_idx = None

        if header:
            low = [h.lower() for h in header]
            for cand in ["gene", "name", "variable", "var", "node", "id"]:
                if cand in low:
                    name_idx = low.index(cand)
                    break
            for cand in ["chi2", "chi2_sum", "chisq", "chisquare", "sum_chi2", "chi2s", "chi2_sums"]:
                if cand in low:
                    chi_idx = low.index(cand)
                    break

        if chi_idx is None:
            probe = rows[data_start]
            for j, c in enumerate(probe):
                if j == name_idx:
                    continue
                if is_number(c):
                    chi_idx = j
                    break

        if chi_idx is None:
            return chi2

        for row in rows[data_start:]:
            if len(row) <= max(name_idx, chi_idx):
                continue
            name = norm_gene(row[name_idx])
            val_str = row[chi_idx].strip().replace(',', '.')
            try:
                val = float(val_str)
            except ValueError:
                continue
            if name in valid_genes:
                chi2[name] = val

    return chi2


def read_single_target_need(path):
    """
    Reads exactly ONE target from gene_with_regulation.csv:
    - takes the first data row
    - need is taken from the last column (typically Activation/Inhibition)
    """
    with open(path, 'r', newline='', encoding='utf-8-sig') as f:
        dialect = sniff_dialect(f)
        reader = csv.reader(f, dialect)
        for row in reader:
            if not row or len(row) < 2:
                continue
            name = row[0].strip()
            if name.lower() in ("name", "gene", "variable", "var"):
                continue
            need = row[-1].strip()
            gene = norm_gene(name)
            return gene, need
    raise RuntimeError("Could not read a target from gene_with_regulation.csv (no data row found).")


def read_state_snapshot_from_chi2_full(path, valid_genes, target_gene, source="ExpData"):
    """
    Reads long-format chi2_full with columns:
    Time,Variable,ExpData,SimulatedData,Std,Chi2,ModelIndex,DataIndex

    Strategy:
    1. Find time where target_gene has maximal Chi2 (largest mismatch point).
    2. Return dict gene -> value (ExpData or SimulatedData) at that time.
    """
    with open(path, 'r', newline='', encoding='utf-8-sig') as f:
        dialect = sniff_dialect(f)
        reader = csv.reader(f, dialect)
        rows = [row for row in reader if row and any(c.strip() for c in row)]
        if not rows:
            raise RuntimeError("chi2_full file is empty.")

    header = [c.strip() for c in rows[0]]
    low = [h.lower() for h in header]
    data = rows[1:]

    def idx(name):
        key = name.lower()
        if key not in low:
            raise RuntimeError(f"chi2_full missing column: {name}")
        return low.index(key)

    i_time = idx("Time")
    i_var  = idx("Variable")
    i_src  = idx(source)
    i_chi  = idx("Chi2")

    tgt = norm_gene(target_gene)

    # 1. find time of max chi2 for target
    best_time = None
    best_chi = -1.0
    for r in data:
        if len(r) <= max(i_time, i_var, i_chi):
            continue
        var = norm_gene(r[i_var])
        if var != tgt:
            continue
        t_str = r[i_time].strip().replace(',', '.')
        chi_str = r[i_chi].strip().replace(',', '.')
        if not (is_number(t_str) and is_number(chi_str)):
            continue
        t = float(t_str)
        chi = float(chi_str)
        if chi > best_chi:
            best_chi = chi
            best_time = t

    if best_time is None:
        raise RuntimeError(f"Target '{target_gene}' not found in chi2_full (Variable column).")

    # 2. build snapshot at that time
    snapshot = {}
    for r in data:
        if len(r) <= max(i_time, i_var, i_src):
            continue
        t_str = r[i_time].strip().replace(',', '.')
        if not is_number(t_str):
            continue
        t = float(t_str)
        if t != best_time:
            continue

        var = norm_gene(r[i_var])
        if var not in valid_genes:
            continue

        v_str = r[i_src].strip().replace(',', '.')
        if not is_number(v_str):
            continue
        snapshot[var] = float(v_str)

    return snapshot, best_time, best_chi


##################### ANALYSIS #############################
def w_eff_from_raw(w_raw: int) -> int:
    return 1 if w_raw != 0 else 0


def score(reg, w_raw, chi2, state=None):
    """
    Edge strength unknown -> do NOT use abs(w_raw) as strength.
    Use binary existence (w_raw != 0) and optional state gating (0..1).

    Score = I(w_raw!=0) * chi2(reg) * state(reg)
    """
    v = chi2.get(reg, None)
    if v is None:
        return None

    w_eff = w_eff_from_raw(w_raw)

    if state is None:
        s = 1.0
    else:
        # If missing state entry, default to 1.0 (do not penalize missing values)
        s = state.get(reg, 1.0)

    return w_eff * v * s


def sort_by_score_desc(items):
    # items: (reg, w_raw, reg_chi, score)
    def key(x):
        s = x[3]
        return -1.0 if s is None else s
    return sorted(items, key=key, reverse=True)


def print_candidates(title, cands, top_k, state=None):
    print(title)
    if not cands:
        print("  -")
        return

    for reg, w_raw, reg_chi, s in cands[:top_k]:
        reg_chi_str = "-" if reg_chi is None else f"{reg_chi:.3f}"
        s_str = "-" if s is None else f"{s:.3f}"

        st_str = "-"
        if state is not None:
            st = state.get(reg, None)
            st_str = "-" if st is None else f"{st:.3f}"

        w_eff = w_eff_from_raw(w_raw)
        sign = "ACT" if w_raw > 0 else "INH"

        print(
            f"  {sign}  {reg:12s}  "
            f"w_raw={w_raw:+d} w_eff={w_eff:d}  "
            f"state={st_str:>7s}  ownχ²={reg_chi_str:>10s}  score={s_str:>10s}"
        )


def rank_upstream_candidates(target_gene, need, genes, M, g2i, chi2, state=None, top_k=15):
    if target_gene not in g2i:
        raise RuntimeError(f"Target '{target_gene}' not found in adjacency gene list.")

    t = g2i[target_gene]
    col = M[:, t]

    # Direct 1-hop incoming edges into target
    ups_act = [(genes[r], int(w)) for r, w in enumerate(col) if w > 0]
    ups_inh = [(genes[r], int(w)) for r, w in enumerate(col) if w < 0]

    cand_act = [(reg, w_raw, chi2.get(reg, None), score(reg, w_raw, chi2, state)) for reg, w_raw in ups_act]
    cand_inh = [(reg, w_raw, chi2.get(reg, None), score(reg, w_raw, chi2, state)) for reg, w_raw in ups_inh]

    cand_act = sort_by_score_desc(cand_act)
    cand_inh = sort_by_score_desc(cand_inh)

    tgt_chi = chi2.get(target_gene, None)
    tgt_chi_str = "-" if tgt_chi is None else f"{tgt_chi:.3f}"

    # Regulation balance: counts (no strength assumption)
    pos_sum = int(np.sum(col > 0))
    neg_sum = int(np.sum(col < 0))
    net = pos_sum - neg_sum

    print("\n==============================")
    print(f"Target: {target_gene}")
    print(f"Need : {need}")
    print(f"Target χ² (from chi2_sums): {tgt_chi_str}")
    print("==============================")
    print(f"Regulation balance (counts): +{pos_sum} / -{neg_sum}  => net={net}")
    if neg_sum == 0 and pos_sum > 0:
        print("Note: Only 'gas' (activators) and no 'brake' (inhibitors) -> runaway/overdrive possible (structure).")
    if pos_sum == 0 and neg_sum > 0:
        print("Note: Only inhibitors and no activators -> may choke the target (structure).")

    want_act = is_need_activation(need)
    want_inh = is_need_inhibition(need)

    if want_inh:
        if not cand_inh:
            print("\nNeed=Inhibition, but there are 0 inhibitors in the adjacency -> missing/wrong edge (structure problem).")
            print("Additionally possible: existing activations are driving the target (overdrive).")
        print_candidates("\nTop upstream INHIBITORS (I(w!=0)*ownχ²*state):", cand_inh, top_k, state=state)
        print_candidates("\nOpposite side: upstream ACTIVATORS (can drive target / overdrive):", cand_act, top_k, state=state)

        cands_by_absw = sorted(cand_act, key=lambda x: abs(x[1]), reverse=True)
        print_candidates("\nActivators by abs(|w_raw|) (informational only):", cands_by_absw, top_k, state=state)

    elif want_act:
        if not cand_act:
            print("\nNeed=Activation, but there are 0 activators in the adjacency -> missing/wrong edge (structure problem).")
            print("Additionally possible: existing inhibition is suppressing the target.")
        print_candidates("\nTop upstream ACTIVATORS (I(w!=0)*ownχ²*state):", cand_act, top_k, state=state)
        print_candidates("\nOpposite side: upstream INHIBITORS (can suppress target):", cand_inh, top_k, state=state)

        cands_by_absw = sorted(cand_inh, key=lambda x: abs(x[1]), reverse=True)
        print_candidates("\nInhibitors by abs(|w_raw|) (informational only):", cands_by_absw, top_k, state=state)

    else:
        print("\nNeed is not cleanly recognizable as Activation/Inhibition. Listing both.")
        print_candidates("\nActivators (I(w!=0)*ownχ²*state):", cand_act, top_k, state=state)
        print_candidates("\nInhibitors (I(w!=0)*ownχ²*state):", cand_inh, top_k, state=state)

    return cand_act, cand_inh


def two_hop_analysis(top_upstream_regs, genes, M, g2i, chi2, state=None, top_k=15):
    """
    For each seed regulator R: list its upstream sources S (S -> R) and rank by I(w!=0)*own chi²*state.
    Keeps self-loops (optionally prints them), but can exclude them from the ranked upstream list.
    """
    print("\n==============================")
    print("2-hop analysis: Who regulates the top upstream candidates?")
    print("==============================")

    for reg_gene in top_upstream_regs:
        if reg_gene not in g2i:
            continue

        t = g2i[reg_gene]
        col = M[:, t]

        ups_all = [(genes[r], int(w)) for r, w in enumerate(col) if w != 0]

        # Self-loop separat halten
        self_loop = [(src, w_raw) for (src, w_raw) in ups_all if src == reg_gene]

        # Optional: Self-loop aus Ranking entfernen, aber nicht aus Anzeige
        if EXCLUDE_SELF_FROM_2HOP_RANK:
            ups = [(src, w_raw) for (src, w_raw) in ups_all if src != reg_gene]
        else:
            ups = ups_all

        cand = [(src, w_raw, chi2.get(src, None), score(src, w_raw, chi2, state)) for (src, w_raw) in ups]
        cand = sort_by_score_desc(cand)

        reg_chi = chi2.get(reg_gene, None)
        reg_chi_str = "-" if reg_chi is None else f"{reg_chi:.3f}"

        pos_sum = int(np.sum(col > 0))
        neg_sum = int(np.sum(col < 0))
        net = pos_sum - neg_sum

        print(f"\nRegulator: {reg_gene}   ownχ²={reg_chi_str}   Balance(counts) +{pos_sum}/-{neg_sum} net={net}")

        # Self-loop optional separat ausgeben
        if SHOW_SELF_LOOPS and self_loop:
            for src, w_raw in self_loop:
                src_chi = chi2.get(src, None)
                src_chi_str = "-" if src_chi is None else f"{src_chi:.3f}"

                st_str = "-"
                if state is not None:
                    st = state.get(src, None)
                    st_str = "-" if st is None else f"{st:.3f}"

                s_val = score(src, w_raw, chi2, state)
                s_str = "-" if s_val is None else f"{s_val:.3f}"

                w_eff = w_eff_from_raw(w_raw)
                sign = "ACT" if w_raw > 0 else "INH"

                print(
                    f"  (self) {sign}  {src:12s}  "
                    f"w_raw={w_raw:+d} w_eff={w_eff:d}  "
                    f"state={st_str:>7s}  ownχ²={src_chi_str:>10s}  score={s_str:>10s}"
                )

        if not cand:
            print("  No upstream edges (excluding self-loop)." if EXCLUDE_SELF_FROM_2HOP_RANK else "  No upstream edges.")
            continue

        for src, w_raw, src_chi, s in cand[:top_k]:
            src_chi_str = "-" if src_chi is None else f"{src_chi:.3f}"
            s_str = "-" if s is None else f"{s:.3f}"

            st_str = "-"
            if state is not None:
                st = state.get(src, None)
                st_str = "-" if st is None else f"{st:.3f}"

            w_eff = w_eff_from_raw(w_raw)
            sign = "ACT" if w_raw > 0 else "INH"

            print(
                f"  {sign}  {src:12s}  "
                f"w_raw={w_raw:+d} w_eff={w_eff:d}  "
                f"state={st_str:>7s}  ownχ²={src_chi_str:>10s}  score={s_str:>10s}"
            )


def choose_twohop_seeds(need, cand_act, cand_inh, n, state=None, min_act_state=0.05):
    """
    Smarter seeds for 2-hop:
    - If Need=Activation but activators are essentially OFF (max state < min_act_state),
      then seed inhibitors instead ("who drives the brake?").
    - If Need=Inhibition but inhibitors missing -> seed activators (drivers/overdrive).
    - Otherwise seed the needed side.
    """
    want_act = is_need_activation(need)
    want_inh = is_need_inhibition(need)

    def max_state(cands):
        if state is None or not cands:
            return None
        ms = None
        for reg, w_raw, reg_chi, s in cands:
            v = state.get(reg, None)
            if v is None:
                continue
            ms = v if ms is None else max(ms, v)
        return ms

    if want_inh and len(cand_inh) == 0:
        primary = cand_act

    elif want_act:
        ms = max_state(cand_act)
        if ms is not None and ms < min_act_state and len(cand_inh) > 0:
            print(f"\n2-hop seed switch: Need=Activation but activators are OFF (max state={ms:.3f} < {min_act_state}).")
            print("Seeding INHIBITORS instead (who drives the brake?).")
            primary = cand_inh
        else:
            primary = cand_act

    elif want_inh:
        primary = cand_inh

    else:
        primary = cand_act + cand_inh

    primary_sorted = sort_by_score_desc(primary)
    seeds = [x[0] for x in primary_sorted if x[3] is not None][:n]
    return seeds

def print_summary_chain(target_gene, need, state_time, state_chi, cand_act, cand_inh, seeds,
                        genes, M, g2i, chi2, state=None, top_n=3):
    """
    Prints a simple, human-readable chain summary:
    - Target + snapshot info
    - Direct (1-hop) top activators/inhibitors
    - If seeded inhibitors (e.g. MTOR): show who regulates the blocker (2-hop)
    """
    print("\n==============================")
    print("SUMMARY (simple chain)")
    print("==============================")

    tgt_chi = chi2.get(target_gene, None)
    tgt_chi_str = "-" if tgt_chi is None else f"{tgt_chi:.3f}"
    print(f"Target: {target_gene} | Need: {need}")
    print(f"Snapshot: t={state_time} (target max pointwise Chi2={state_chi:.3f}) | Target Chi2(sum)={tgt_chi_str}")

    def fmt_item(x):
        reg, w_raw, reg_chi, sc = x
        st = None if state is None else state.get(reg, None)
        st_str = "-" if st is None else f"{st:.3f}"
        reg_chi_str = "-" if reg_chi is None else f"{reg_chi:.3f}"
        sc_str = "-" if sc is None else f"{sc:.3f}"
        sign = "ACT" if w_raw > 0 else "INH"
        return f"{reg} [{sign}, state={st_str}, chi2={reg_chi_str}, score={sc_str}]"

    # 1-hop direct influences on target
    top_act = [x for x in cand_act if x[3] is not None][:top_n]
    top_inh = [x for x in cand_inh if x[3] is not None][:top_n]

    if top_act:
        print("\nDirect drivers (ACT -> target):")
        for x in top_act:
            print("  - " + fmt_item(x))
    else:
        print("\nDirect drivers (ACT -> target): none")

    if top_inh:
        print("\nDirect blockers (INH -> target):")
        for x in top_inh:
            print("  - " + fmt_item(x))
    else:
        print("\nDirect blockers (INH -> target): none")

    # If we seeded something, summarize 2-hop for the first seed only (usually the main blocker)
    if not seeds:
        print("\nNo 2-hop seeds -> no upstream chain.")
        return

    main_seed = seeds[0]
    if main_seed not in g2i:
        print(f"\nSeed '{main_seed}' not in network -> cannot expand chain.")
        return

    # Determine relation of seed to target (was it in inhibitors list or activators list?)
    seed_role = None
    for reg, w_raw, *_ in cand_inh:
        if reg == main_seed:
            seed_role = "BLOCKER (INH -> target)"
            break
    if seed_role is None:
        for reg, w_raw, *_ in cand_act:
            if reg == main_seed:
                seed_role = "DRIVER (ACT -> target)"
                break
    if seed_role is None:
        seed_role = "SEED"

    print(f"\nExpanded 2-hop for seed: {main_seed} ({seed_role})")

    # 2-hop upstream of seed (sources -> seed)
    t = g2i[main_seed]
    col = M[:, t]
    ups_all = [(genes[r], int(w)) for r, w in enumerate(col) if w != 0]

    # split into activators/inhibitors of the seed
    ups_act = []
    ups_inh = []
    for src, w_raw in ups_all:
        if EXCLUDE_SELF_FROM_2HOP_RANK and src == main_seed:
            continue
        item = (src, w_raw, chi2.get(src, None), score(src, w_raw, chi2, state))
        if w_raw > 0:
            ups_act.append(item)
        else:
            ups_inh.append(item)

    ups_act = sort_by_score_desc([x for x in ups_act if x[3] is not None])
    ups_inh = sort_by_score_desc([x for x in ups_inh if x[3] is not None])

    if ups_act:
        print(f"\nWho drives {main_seed} (ACT -> {main_seed}):")
        for x in ups_act[:top_n]:
            print("  - " + fmt_item(x))
    else:
        print(f"\nWho drives {main_seed} (ACT -> {main_seed}): none")

    if ups_inh:
        print(f"\nWho brakes {main_seed} (INH -> {main_seed}):")
        for x in ups_inh[:top_n]:
            print("  - " + fmt_item(x))
    else:
        print(f"\nWho brakes {main_seed} (INH -> {main_seed}): none")

    # Final one-liner
    print("\nOne-line interpretation:")
    if seed_role.startswith("BLOCKER"):
        print(f"  {target_gene} needs activation likely because {main_seed} is blocking it; "
              f"{main_seed} is driven by {', '.join([x[0] for x in ups_act[:top_n]]) or 'n/a'} "
              f"and braked by {', '.join([x[0] for x in ups_inh[:top_n]]) or 'n/a'}.")
    else:
        print(f"  {target_gene} is influenced by seed {main_seed}; check upstream drivers/brakes above.")


############### MAIN ##############
if __name__ == "__main__":
    genes, M = read_adjacency(ADJ_FILE)
    g2i = {g: i for i, g in enumerate(genes)}

    chi2 = read_chi2_sums(CHI2_SUM_FILE, set(genes))
    if not chi2:
        raise RuntimeError(
            "No chi² values loaded. Check chi2_sums_per_variable_full.csv (columns/delimiter)."
        )

    target_gene, need = read_single_target_need(GENE_WITH_REG_FILE)

    # load state snapshot at time where TARGET has max mismatch (max Chi2) in chi2_full
    state = None
    state_time = None
    state_chi = None
    try:
        state, state_time, state_chi = read_state_snapshot_from_chi2_full(
            CHI2_FULL_FILE, set(genes), target_gene, source=STATE_SOURCE
        )
        print(f"\nState snapshot loaded: {STATE_SOURCE} at time={state_time} (target max Chi2={state_chi:.3f})")
    except Exception as e:
        print(f"\nState snapshot not loaded (continuing without state gating): {e}")

    cand_act, cand_inh = rank_upstream_candidates(
        target_gene=target_gene,
        need=need,
        genes=genes,
        M=M,
        g2i=g2i,
        chi2=chi2,
        state=state,
        top_k=TOP_K
    )

    seeds = []
    if DO_TWO_HOP:
        seeds = choose_twohop_seeds(
            need, cand_act, cand_inh, TWO_HOP_TOP_UPSTREAM, state=state
        )
        print("2-hop seeds:", seeds)

        if seeds:
            two_hop_analysis(
                seeds, genes, M, g2i, chi2, state=state, top_k=TWO_HOP_TOP_PER_NODE
            )
        else:
            print("\n2-hop skipped: No suitable upstream candidates with their own chi²/score found.")

    # Summary chain at the end (works even if 2-hop is disabled or seeds are empty)
    # If state snapshot was not loaded, state_time/state_chi may be None -> print_summary_chain should handle it.
    if state_time is None:
        # fallback values to avoid formatting errors if your summary assumes floats
        state_time = float("nan")
    if state_chi is None:
        state_chi = float("nan")

    print_summary_chain(
        target_gene=target_gene,
        need=need,
        state_time=state_time,
        state_chi=state_chi,
        cand_act=cand_act,
        cand_inh=cand_inh,
        seeds=seeds,
        genes=genes,
        M=M,
        g2i=g2i,
        chi2=chi2,
        state=state,
        top_n=3
    )
