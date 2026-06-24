#!/usr/bin/env python3
"""Per-epitope apex feature table for the structural-basis supplement (Appendix A5, Deliverable 1).

For each benchmark epitope (NLV, LLW, LLL, YLQ, GLC; all HLA-A*02:01) this derives, from the
epitope's VDJdb2026 TRB reference CDR3s (``_feat_probe.ref_table``), a small set of *apex* descriptors
that explain why some epitopes are detectable by identity scoring and others need a structural
(Miyazawa-Jernigan) potential:

  apex_hydropathy : mean Kyte-Doolittle hydropathy over the CDR3 apex (central 5 residues),
                    averaged across reference CDR3s. High + featureless -> "hydrophobic ridge".
  pep_hydropathy  : Kyte-Doolittle hydropathy of the peptide's TCR-facing core (positions 4-6),
                    characterising the epitope surface the CDR3 apex contacts (public peptide aa).
  mj_strength     : mean seqtree ``structural`` (MJ sidechain-volume/hydropathy) similarity of each
                    apex residue to the canonical hydrophobic core {A,I,L,M,F,V,W}, averaged over
                    reference CDR3s. High = apex sits in the hydrophobic cluster the MJ matrix groups.
  apex_entropy    : mean per-position Shannon entropy (bits) over apex columns (length-pooled).
                    Low entropy = "featured" (conserved motif); high entropy = "featureless".
  dom_v_frac      : fraction of reference CDR3s using the single most common TRBV gene.
                    High = a dominant public V (featured, identity-recognisable).
  best_roc        : vdjmatch best-chain detection ROC-AUC for the task (from the manuscript
                    detection_bychain.tsv / detection_TRB.tsv), to tie features to detectability.

NUMBERS / AGGREGATES ONLY -- no raw CDR3 sequences are written. VDJdb reference sequences are public.

Run from repo root with the project venv:
    ./.venv/bin/python bench/epitope_features.py
"""
from __future__ import annotations

import math
import sys
from collections import Counter
from pathlib import Path

BENCH = Path(__file__).resolve().parent
sys.path.insert(0, str(BENCH))
import _feat_probe as fp  # noqa: E402

# Manuscript repo (write targets + detection ROC source).
MS = Path("/Users/mikesh/vcs/manuscripts/2026-vdjmatch")
RESULTS = MS / "benchmarks" / "results"
OUT_TSV = RESULTS / "epitope_features.tsv"
OUT_TEX = MS / "appendix" / "_tab_epitope_features.tex"
STRUCTURAL_INC = "/Users/mikesh/vcs/code/seqtree/src/structural.inc"

TASKS = ["NLV", "LLW", "LLL", "YLQ", "GLC"]

# Kyte-Doolittle hydropathy scale (J. Mol. Biol. 1982; 157:105-132).
KD = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8,
    "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2,
    "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}
# Hydrophobic core used to summarise MJ interaction strength of an apex residue.
HYDRO_CORE = ("A", "I", "L", "M", "F", "V", "W")
APEX = 5  # central window width


def structural_sim(inc: str) -> dict[tuple[str, str], int]:
    """Parse seqtree's ``kStructural[24*24]`` similarity matrix (0-10) -> {(a,b): sim}."""
    import re
    order = "ARNDCQEGHILKMFPSTWYVBZX*"
    txt = Path(inc).read_text().split("kStructural[24 * 24] = {", 1)[-1]
    nums = [int(x) for x in re.findall(r"-?\d+", txt.split("};")[0])][: 24 * 24]
    return {(order[i], order[j]): nums[i * 24 + j] for i in range(24) for j in range(24)}


def apex_residues(cdr3: str, width: int = APEX) -> str:
    """Central ``width`` residues of a CDR3 (the apex loop); shorter CDR3s returned whole."""
    n = len(cdr3)
    if n <= width:
        return cdr3
    start = (n - width) // 2
    return cdr3[start: start + width]


def shannon_bits(counts) -> float:
    """Shannon entropy (bits) of a residue-count mapping."""
    tot = sum(counts.values())
    if tot <= 0:
        return 0.0
    h = 0.0
    for c in counts.values():
        if c:
            p = c / tot
            h -= p * math.log2(p)
    return h


def epitope_row(task: str, sim: dict, best_roc: float) -> dict:
    """Compute the apex descriptor row for one task from its VDJdb2026 TRB reference."""
    ref = fp.ref_table(task, "TRB")
    cdr3s = ref["cdr3"].to_list()
    vs = ref["v"].to_list()
    n = len(cdr3s)

    hydro_vals, mj_vals = [], []
    # Length-pooled per-position apex columns for entropy: align apex windows left-to-right.
    columns: list[Counter] = [Counter() for _ in range(APEX)]
    for s in cdr3s:
        ap = apex_residues(s)
        hyd = [KD[a] for a in ap if a in KD]
        if hyd:
            hydro_vals.append(sum(hyd) / len(hyd))
        # mean structural similarity of each apex residue to the hydrophobic core
        per = []
        for a in ap:
            sims = [sim.get((a, h), 0) for h in HYDRO_CORE]
            per.append(sum(sims) / len(sims))
        if per:
            mj_vals.append(sum(per) / len(per))
        # entropy columns (only count windows of full APEX width to keep columns comparable)
        if len(ap) == APEX:
            for i, a in enumerate(ap):
                columns[i][a] += 1

    apex_hydropathy = sum(hydro_vals) / len(hydro_vals)
    mj_strength = sum(mj_vals) / len(mj_vals)
    # peptide TCR-facing core (positions 4-6, 0-indexed 3:6) Kyte-Doolittle mean
    pep = fp.E[task]
    pcore = pep[3:6]
    pep_hydropathy = sum(KD[a] for a in pcore if a in KD) / len(pcore)
    col_h = [shannon_bits(c) for c in columns if sum(c.values()) > 0]
    apex_entropy = sum(col_h) / len(col_h) if col_h else 0.0
    vc = Counter(vs)
    dom_v_frac = vc.most_common(1)[0][1] / n if n else 0.0

    return {
        "epitope": task,
        "epitope_aa": fp.E[task],
        "n_ref": n,
        "apex_hydropathy": round(apex_hydropathy, 3),
        "pep_hydropathy": round(pep_hydropathy, 3),
        "mj_strength": round(mj_strength, 3),
        "apex_entropy": round(apex_entropy, 3),
        "dom_v_frac": round(dom_v_frac, 3),
        "best_roc": best_roc,
    }


def best_roc_table() -> dict[str, float]:
    """vdjmatch best-chain detection ROC-AUC per task, from the manuscript detection TSVs.
    Uses detection_bychain.tsv (paired/TRA/TRB) where available, else detection_TRB.tsv."""
    best = {t: 0.0 for t in TASKS}

    def scan(path: Path, method_col: str):
        if not path.exists():
            return
        lines = path.read_text().strip().splitlines()
        header = lines[0].split("\t")
        ti = header.index("task")
        mi = header.index(method_col)
        ri = header.index("roc_auc")
        for ln in lines[1:]:
            f = ln.split("\t")
            t = f[ti]
            if t in best and f[mi] == "vdjmatch":
                best[t] = max(best[t], float(f[ri]))

    scan(RESULTS / "detection_TRB.tsv", "method")
    scan(RESULTS / "detection_bychain.tsv", "method")  # paired chains for YLQ/GLC
    return best


def main():
    sim = structural_sim(STRUCTURAL_INC)
    rocs = best_roc_table()
    rows = [epitope_row(t, sim, rocs[t]) for t in TASKS]

    cols = ["epitope", "epitope_aa", "n_ref", "apex_hydropathy", "pep_hydropathy",
            "mj_strength", "apex_entropy", "dom_v_frac", "best_roc"]
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TSV.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"wrote {OUT_TSV}")
    for r in rows:
        print("  " + "  ".join(f"{c}={r[c]}" for c in cols))

    write_tex(rows)
    print(f"wrote {OUT_TEX}")


def write_tex(rows):
    """booktabs LaTeX table; plain \\caption{} (the supplement auto-numbers)."""
    lines = [
        r"\begin{table}[h]\centering\small",
        (r"\caption{Apex descriptors of the five HLA-A*02:01 benchmark epitopes, derived from their "
         r"VDJdb2026 TRB reference repertoires, alongside the best-chain \textsc{vdjmatch} detection "
         r"ROC-AUC. \emph{Apex hydrop.} is the mean Kyte--Doolittle value over the central five "
         r"CDR3$\beta$ residues; \emph{Pep.\ hydrop.} is the same scale over the peptide's TCR-facing "
         r"core (positions 4--6); \emph{MJ str.} is the mean Miyazawa--Jernigan (\texttt{structural}) "
         r"similarity of apex residues to the hydrophobic core; \emph{Apex $H$} is the mean "
         r"per-position Shannon entropy (bits) over the apex window; \emph{Dom.\ V} is the fraction of "
         r"references using the most common TRBV gene. A hydrophobic peptide core with high apex "
         r"entropy and a diffuse V usage marks a featureless hydrophobic ridge (e.g.\ NLV), where the "
         r"structural potential helps; a polar core with lower apex entropy marks a featured, "
         r"identity-recognisable motif (e.g.\ YLQ).}"),
        r"\begin{tabular}{@{}llccccccc@{}}\toprule",
        (r"Epitope & Sequence & $n_{\mathrm{ref}}$ & Apex hydrop. & Pep.\ hydrop. & MJ str. "
         r"& Apex $H$ (bits) & Dom.\ V & Best ROC \\\midrule"),
    ]
    for r in rows:
        lines.append(
            f"{r['epitope']} & \\texttt{{{r['epitope_aa']}}} & {r['n_ref']} & "
            f"{r['apex_hydropathy']:.2f} & {r['pep_hydropathy']:.2f} & {r['mj_strength']:.2f} & "
            f"{r['apex_entropy']:.2f} & {r['dom_v_frac']:.2f} & {r['best_roc']:.3f} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}\end{table}", ""]
    OUT_TEX.write_text("\n".join(lines))


if __name__ == "__main__":
    main()
