#!/usr/bin/env python3
"""Leave-one-out-by-epitope evaluation of VDJAM and region-aware (NDN-weighted) scoring.

VDJAM is derived from same-antigen substitution pairs, so testing on the same epitopes is circular.
For each held-out epitope E* we re-derive the NDN substitution matrix from **all other** epitopes,
then run VDJdb-vs-VDJdb retrieval on E*: a single unit-cost search (scope 2, no indels) fixes the
candidate set, and each (query, hit) pair is scored four ways and ranked. Micro PR-AUC for
"hit shares E*" measures each scoring's specificity resolution:

  unit          : Hamming (substitution count)
  BLOSUM62      : sum of per-substitution dissimilarity
  VDJAM (flat)  : sum of NDN-VDJAM dissimilarity, all positions equal
  VDJAM region  : NDN-VDJAM dissimilarity weighted by (1 - germline retention) per position
                  -> germline V/J flanks discounted, NDN core emphasized

    python bench/loo_vdjam.py --chain TRB --top 8 --max-queries 600
"""
from __future__ import annotations

import argparse
import os
import statistics as st

import polars as pl
from Bio.Align import substitution_matrices
from seqtree import Index, SearchParams

import _bench
from metrics import pr_auc, roc_auc, pr_auc_balanced  # noqa: E402
from vdjmatch import db
from vdjmatch.match import regions
from gen_vdjam import hamming1_events, position_background, estimate_matrix, AA  # noqa: E402


def dissim_from_scores(scores: dict[tuple[str, str], float]) -> dict[tuple[str, str], float]:
    """similarity log-odds -> non-negative dissimilarity (most-similar pair = 0)."""
    off = {k: v for k, v in scores.items() if k[0] != k[1]}
    hi = max(off.values())
    return {k: hi - v for k, v in off.items()}


def named_dissim(name: str) -> dict[tuple[str, str], float]:
    """similarity matrix (BLOSUM62 / PAM250 / ...) -> non-negative off-diagonal dissimilarity."""
    b = substitution_matrices.load(name)
    vals = {(a, c): float(b[a, c]) for a in AA for c in AA if a != c}
    hi = max(vals.values())
    return {k: hi - v for k, v in vals.items()}


# standard genetic code (DNA codons -> 1-letter AA); '*' = stop
_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def codon_dissim() -> dict[tuple[str, str], float]:
    """VDJAMr: a substitution matrix from the genetic code alone -- how mutationally accessible one
    amino acid is from another. M(a,b) = P(a single random nucleotide substitution of a uniformly
    random codon of a yields a codon of b); symmetrised, then converted to off-diagonal dissimilarity
    (codon-adjacent pairs = cheap). This is the recombination/mutation null: it scores the substitution
    propensity built into the code, with no chemistry or selection."""
    from collections import defaultdict
    codons = defaultdict(list)
    for c, a in _CODE.items():
        if a != "*":
            codons[a].append(c)
    nts = "ACGT"
    m = {}
    for a in AA:
        for b in AA:
            tot = hit = 0
            for c in codons[a]:
                for i in range(3):
                    for nt in nts:
                        if nt == c[i]:
                            continue
                        tot += 1
                        hit += _CODE.get(c[:i] + nt + c[i + 1:]) == b
            m[(a, b)] = hit / tot if tot else 0.0
    sim = {(a, b): (m[(a, b)] + m[(b, a)]) / 2 for a in AA for b in AA}
    off = {k: v for k, v in sim.items() if k[0] != k[1]}
    hi = max(off.values())
    return {k: hi - v for k, v in off.items()}


def structural_dissim(inc=None) -> dict | None:
    """seqtree's TeXshade structural similarity matrix (sidechain volume + hydropathy) -> dissim.
    Parsed from the seqtree source .inc (24-symbol order ARNDCQEGHILKMFPSTWYVBZX*); None if absent.
    Path from arg or $STRUCTURAL_INC (e.g. <seqtree>/src/structural.inc)."""
    import re
    from pathlib import Path
    inc = inc or os.environ.get("STRUCTURAL_INC", "")
    if not inc or not Path(inc).exists():
        return None
    order = "ARNDCQEGHILKMFPSTWYVBZX*"
    txt = Path(inc).read_text().split("kStructural[24 * 24] = {", 1)[-1]
    nums = [int(x) for x in re.findall(r"-?\d+", txt.split("};")[0])][:24 * 24]
    sim = {(order[i], order[j]): nums[i * 24 + j] for i in range(24) for j in range(24)}
    vals = {(a, c): sim[(a, c)] for a in AA for c in AA if a != c}
    hi = max(vals.values())
    return {k: hi - v for k, v in vals.items()}


def score_pairs(qseqs, qv, qj, cand, ref_cdr3, ref_epi, true_epi, chain, ret, dis, wmode=None):
    """(label, -score) pairs over all candidate hits; score = sum of (weighted) dissimilarity.
    wmode: None (flat) | 'region' (germline-retention) | 'possig' (positional significance)."""
    out = []
    for i, hl in enumerate(cand):
        if wmode == "region":
            w = regions.position_weights(len(qseqs[i]), qv[i], qj[i], chain, ret)
        elif wmode == "possig":
            w = regions.significance_weights(len(qseqs[i]))
        else:
            w = None
        for ri in hl:
            r = ref_cdr3[ri]
            if len(r) != len(qseqs[i]):
                continue  # subs-only candidates are same length
            if r == qseqs[i]:
                continue  # exact self / identical clone -- not an informative retrieval
            sc = 0.0
            for p, (a, b) in enumerate(zip(qseqs[i], r)):
                if a != b:
                    sc += (w[p] if w else 1.0) * dis.get((a, b), 0.0)
            out.append((1 if ref_epi[ri] == true_epi else 0, -sc))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc", default=None,
                    help="VDJdb export TSV (default: $VDJDB_SAMPLE or the HF-pinned release)")
    ap.add_argument("--species", default="HomoSapiens")
    ap.add_argument("--chain", default="TRB")
    ap.add_argument("--min-epi", type=int, default=30)
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--max-queries", type=int, default=600)
    ap.add_argument("--subs", type=int, default=2, help="candidate search substitution budget")
    args = ap.parse_args()

    vdj = db.load(_bench.source(args.pmhc), species=args.species).filter(pl.col("gene") == args.chain)
    uc = _bench.long_list(vdj, cap=3000, min_n=args.min_epi)  # composition-controlled clonotypes
    bg = position_background(uc["cdr3"].unique().to_list())
    ret = regions.load_retention()
    blo = named_dissim("BLOSUM62")
    pam = named_dissim("PAM250")
    struc = structural_dissim() or blo  # fall back to BLOSUM if the seqtree source isn't present
    codon = codon_dissim()              # VDJAMr: genetic-code / mutation-accessibility null
    unit = {(a, c): 1.0 for a in AA for c in AA if a != c}

    refs = uc.group_by("cdr3").agg(pl.col("epitope").first())
    ref_cdr3 = refs["cdr3"].to_list()
    ref_epi = refs["epitope"].to_list()
    index = Index.build(ref_cdr3, "aa")
    cand_params = SearchParams(max_subs=args.subs, max_total_edits=args.subs, engine="seqtm")

    sizes = (uc.group_by("epitope").agg(pl.col("cdr3").n_unique().alias("n"))
               .filter(pl.col("n") >= args.min_epi).sort("n", descending=True))
    held = sizes["epitope"].to_list()[:args.top]
    print(f"species={args.species} chain={args.chain}; held-out epitopes={len(held)}; subs={args.subs}")
    print(f"{'epitope':13}{'n':>5}{'unit':>8}{'BLOSUM':>8}{'PAM250':>8}{'struct':>8}{'VDJAM':>8}"
          f"{'VDJAMr':>8}{'VDJAM_reg':>10}{'B+possig':>10}")
    acc = {k: [] for k in ("unit", "blosum", "pam", "struct", "vdjam", "vdjamr", "region", "bps")}
    for epi in held:
        ev = []
        for e in sizes["epitope"].to_list():
            if e == epi:
                continue
            cds = uc.filter(pl.col("epitope") == e)["cdr3"].unique().to_list()
            ev += [x for x in hamming1_events(cds, args.min_epi) if 4 <= x[2] < x[3] - 6]  # NDN
        vdis = dissim_from_scores(estimate_matrix(ev, bg))

        q = uc.filter(pl.col("epitope") == epi).unique("cdr3").head(args.max_queries)
        qs, qv, qj = q["cdr3"].to_list(), q["v"].to_list(), q["j"].to_list()
        cand = [[h.ref_id for h in hl] for hl in index.search_batch(qs, cand_params, 0)]
        sp = lambda dis, w: pr_auc_balanced(score_pairs(qs, qv, qj, cand, ref_cdr3, ref_epi, epi,  # noqa: E731
                                                        args.chain, ret, dis, w))
        row = {"unit": sp(unit, None), "blosum": sp(blo, None), "pam": sp(pam, None),
               "struct": sp(struc, None), "vdjam": sp(vdis, None), "vdjamr": sp(codon, None),
               "region": sp(vdis, "region"), "bps": sp(blo, "possig")}
        for k in acc:
            acc[k].append(row[k])
        print(f"{epi:13}{len(qs):>5}{row['unit']:>8.3f}{row['blosum']:>8.3f}{row['pam']:>8.3f}"
              f"{row['struct']:>8.3f}{row['vdjam']:>8.3f}{row['vdjamr']:>8.3f}{row['region']:>10.3f}"
              f"{row['bps']:>10.3f}")
    mean = {k: st.mean(acc[k]) for k in acc}
    print(f"\nmean balanced PR-AUC:  unit {mean['unit']:.3f}  BLOSUM {mean['blosum']:.3f}  PAM250 {mean['pam']:.3f}  "
          f"struct {mean['struct']:.3f}  VDJAM {mean['vdjam']:.3f}  VDJAMr {mean['vdjamr']:.3f}  "
          f"VDJAM-region {mean['region']:.3f}  BLOSUM+possig {mean['bps']:.3f}")
    print(f"VDJAMr vs BLOSUM: {st.mean(acc['vdjamr']) - st.mean(acc['blosum']):+.3f}  | "
          f"VDJAMr>BLOSUM in {sum(1 for a,b in zip(acc['vdjamr'],acc['blosum']) if a>b)}/{len(held)}")
    print(f"region vs flat VDJAM: {st.mean(acc['region']) - st.mean(acc['vdjam']):+.3f}  | "
          f"region vs BLOSUM: {st.mean(acc['region']) - st.mean(acc['blosum']):+.3f}  | "
          f"region>BLOSUM in {sum(1 for a,b in zip(acc['region'],acc['blosum']) if a>b)}/{len(held)}")
    print(f"BLOSUM+possig vs BLOSUM: {st.mean(acc['bps']) - st.mean(acc['blosum']):+.3f}  | "
          f"B+possig>BLOSUM in {sum(1 for a,b in zip(acc['bps'],acc['blosum']) if a>b)}/{len(held)}")


if __name__ == "__main__":
    main()
