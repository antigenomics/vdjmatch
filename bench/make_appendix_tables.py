#!/usr/bin/env python3
"""Generate Appendix B: cross-method benchmark tables (LaTeX) from the benchmark outputs.

Tables: (1) datasets/conditions + purpose + n=, (2) methods: algorithm + classifier/clustering,
(3) results: mean +- sd across epitopes, plus by-epitope tables for C1, C4, C5 (with the confusion
matrix at each method's significance call). Reads only the numeric benchmark summaries — never the
manuscript test data. Writes appendix/appendix_b_benchmark.tex.

    python bench/make_appendix_tables.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))

OUT = Path("bench/out/benchmark")
TEX = Path("appendix/appendix_b_benchmark.tex")

# (id, purpose, reference, mode, loci)  — n= filled from the run's dataset_n.tsv
CONDITIONS = [
    ("C1", "Discriminate two co-presented HLA-A*02 epitopes (YFV LLW vs BST-2 LLL) in a sorted sample",
     "VDJdb2026 (exact-filtered)", "classify", "TRB"),
    ("C2", "Detect CMV NLV-specific TCRs vs tetramer-negative controls in a sorted sample",
     "VDJdb2026 (exact-filtered)", "classify", "TRA, TRB"),
    ("C3", "Rank functionally-validated (padj<1e-5) vs non-validated TCRvdb clonotypes, GLC + YLQ",
     "VDJdb2026 (exact-filtered)", "classify", "TRA, TRB, paired"),
    ("C4", "Per-epitope one-vs-rest within the VDJdb2025 reference-replicated shortlist",
     "VDJdb2025 shortlist (>=2 refs), self", "classify", "TRA, TRB, paired"),
    ("C5", "Per-epitope one-vs-rest on clonotypes NEW in VDJdb2026 (temporal hold-out)",
     "VDJdb2025 shortlist (temporal)", "classify", "TRA, TRB, paired"),
    ("C6", "False-positive rate: OLGA-generated TCRs (no true specificity), any hit is a FP",
     "VDJdb2025", "FP-only", "TRA (s5), TRB (s4)"),
    ("Cl", "Clustering quality on the VDJdb2026 shortlist (retention / purity)",
     "VDJdb2026 shortlist, self", "cluster", "TRA, TRB, paired"),
]

# (method, type, chains, algorithm)
METHODS = [
    ("vdjmatch", "classify + cluster", "A / B / paired",
     "First-hit control-calibrated E-value: nearest reference within an adaptive edit ball, epitope "
     "enrichment as a Poisson tail against an OLGA/repertoire background; neighbour-vote clustering."),
    ("tcrdist3", "classify + cluster", "A / B / paired",
     "V-gene-aware weighted CDR1/2/2.5/3 Hamming distance (TCRdist); 1-NN / k-NN label transfer for "
     "classification, distance-threshold single-linkage for clustering. Run in parallel (parasail)."),
    ("tcrmatch", "classify", "B only",
     "BLOSUM62 k-mer similarity kernel between trimmed CDR3$\\beta$ and a labelled reference set; "
     "score 0--1, nearest-reference epitope transfer."),
    ("MixTCRpred", "classify", "paired (A+B)",
     "Per-pMHC transformer over CDR1/2/3 of both chains; one pretrained model per epitope (146 pMHCs; "
     "covers NLV/LLW/GLC/YLQ/GIL, not LLL). Paired input required."),
    ("NetTCR-2.2", "classify", "B (+paired)",
     "Convolutional network on CDR3 (and optionally CDR1/2) + peptide; pan-specific binding "
     "probability for arbitrary (TCR, peptide) pairs."),
    ("GLIPH2", "cluster", "B",
     "Groups CDR3$\\beta$ by shared global and local (motif) similarity with V-gene/length/HLA "
     "enrichment; outputs specificity clusters (no per-pair score)."),
    ("imw-DETECT", "classify", "as provided",
     "ImmuneWatch DETECT (commercial); pre-computed per-sample binding scores supplied as .ods."),
    ("ERGO-II", "classify", "B (+paired)",
     "Autoencoder/LSTM TCR embedding + peptide; supervised TCR--peptide binding probability."),
]


def esc(s):
    return s.replace("&", "\\&").replace("_", "\\_").replace("%", "\\%")


def tbl(header, rows, colspec, caption, label, small=True):
    out = ["\\begin{table}[htbp]\\centering", f"\\caption{{{caption}}}\\label{{{label}}}"]
    if small:
        out.append("\\footnotesize")
    out.append(f"\\begin{{tabular}}{{{colspec}}}\\toprule")
    out.append(" & ".join(header) + " \\\\\\midrule")
    for r in rows:
        out.append(" & ".join(str(x) for x in r) + " \\\\")
    out += ["\\bottomrule\\end{tabular}", "\\end{table}", ""]
    return "\n".join(out)


def fmt(df, cond, loci=None):
    """summary.tsv rows for a condition -> method x metric mean(sd) cells."""
    d = df.filter(pl.col("condition") == cond)
    if loci:
        d = d.filter(pl.col("locus").is_in(loci))
    rows = []
    for (loc, meth), g in d.group_by(["locus", "method"], maintain_order=True):
        cells = {m: f"{v:.2f}" for m, v in zip(g["metric"], g["mean"])}
        rows.append((esc(meth), loc, cells.get("roc_auc", "--"), cells.get("pr_auc", "--"),
                     cells.get("f1", "--"), cells.get("fp", "--")))
    return rows


def main():
    n_df = pl.read_csv(OUT / "dataset_n.tsv", separator="\t") if (OUT / "dataset_n.tsv").exists() else None
    summ = pl.read_csv(OUT / "summary.tsv", separator="\t") if (OUT / "summary.tsv").exists() else None
    by = pl.read_csv(OUT / "by_epitope.tsv", separator="\t") if (OUT / "by_epitope.tsv").exists() else None

    parts = [r"\section{Datasets and conditions}", ""]
    # Table 1: conditions + n=
    nrows = []
    for cid, purpose, ref, mode, loci in CONDITIONS:
        n = ""
        if n_df is not None:
            sel = n_df.filter(pl.col("condition") == cid)
            if sel.height:
                n = "; ".join(f"{r['locus']}: {r['n_queries']}q/{r['n_epitopes']}ep"
                              for r in sel.iter_rows(named=True))
        nrows.append((cid, esc(loci), esc(ref), mode, esc(n) or "--", esc(purpose)))
    parts.append(tbl(["ID", "Loci", "Reference", "Mode", "$n$ (pos/neg)", "Purpose"], nrows,
                     "llp{2.6cm}lp{2.4cm}p{4.2cm}", "Benchmark conditions, references and sample sizes.",
                     "tab:conditions"))

    # Table 2: methods
    mrows = [(esc(m), t, esc(c), esc(a)) for m, t, c, a in METHODS]
    parts.append(r"\section{Methods}")
    parts.append(tbl(["Method", "Mode", "Chains", "Algorithm"], mrows,
                     "llp{2.1cm}p{8.2cm}", "Compared methods: algorithm and classifier/clustering mode.",
                     "tab:methods"))

    # Table 3: results mean+-sd (classification conditions)
    parts.append(r"\section{Results}")
    if summ is not None:
        for cid in ("C1", "C2", "C3", "C4", "C5"):
            rows = fmt(summ, cid)
            if rows:
                parts.append(tbl(["Method", "Locus", "ROC", "PR", "F1", "FP"], rows,
                                 "llllll", f"Condition {cid}: mean across epitopes (ROC/PR/F1/FP).",
                                 f"tab:res{cid}"))
        # by-epitope for C1, C4, C5
        if by is not None:
            for cid in ("C1", "C4", "C5"):
                d = by.filter((pl.col("condition") == cid) &
                              pl.col("metric").is_in(["roc_auc", "pr_auc", "f1", "fp"]))
                if not d.height:
                    continue
                piv = (d.pivot(values="value", index=["method", "locus", "epitope"], on="metric",
                               aggregate_function="first").sort(["method", "locus", "epitope"]))
                cols = [c for c in ("roc_auc", "pr_auc", "f1", "fp") if c in piv.columns]
                rows = [(esc(r["method"]), esc(r["epitope"]), r["locus"],
                         *[f"{r[c]:.2f}" if r[c] is not None else "--" for c in cols])
                        for r in piv.iter_rows(named=True)]
                parts.append(tbl(["Method", "Epitope", "Locus"] + [c.split("_")[0].upper() for c in cols],
                                 rows, "lll" + "l" * len(cols),
                                 f"Condition {cid}: by-epitope per method (FP = false-positive rate at "
                                 f"the method's significance call).", f"tab:by{cid}"))
    if summ is not None and "C6" in summ["condition"].to_list() or (OUT / "fp_rates.tsv").exists():
        fp = pl.read_csv(OUT / "fp_rates.tsv", separator="\t")
        rows = [(esc(r["method"]), r["locus"], f"{r['fp_rate']*100:.2f}\\%", r["n"]) for r in fp.iter_rows(named=True)]
        parts.append(tbl(["Method", "Locus", "FP rate", "$n$"], rows, "llll",
                         "Condition C6: OLGA false-positive rate (any hit is a FP; lower is better).",
                         "tab:resC6"))

    body = "\n".join(parts)
    doc = (r"\documentclass[11pt]{article}" "\n"
           r"\usepackage[margin=2.2cm]{geometry}\usepackage{booktabs}\usepackage{amsmath}" "\n"
           r"\renewcommand{\thesection}{B.\arabic{section}}" "\n"
           r"\title{\scshape Appendix B: Cross-method specificity benchmark}\author{}\date{}" "\n"
           r"\begin{document}\maketitle" "\n"
           r"\noindent This appendix accompanies Appendix A (scoring and substitution matrices). "
           r"All numbers are computed by \texttt{bench/benchmark.py}; no raw sample data is reproduced." "\n\n"
           + body + "\n" r"\end{document}" "\n")
    TEX.write_text(doc)
    print("wrote", TEX)


if __name__ == "__main__":
    main()
