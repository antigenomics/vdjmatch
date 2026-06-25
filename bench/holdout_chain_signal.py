"""Predict, from molecular biology, which chain carries detectable specificity signal.

Per (epitope, chain) we measure leakage-safe REPERTOIRE CONVERGENCE from the vdjdb2026 reference alone
(no test labels, no search): apex Shannon entropy (the loop tip makes all peptide contacts; germline ends
never do), dominant-V fraction, and top-motif coverage. A convergent repertoire (low apex entropy / dominant
V) means binders agree at the contact positions -> a sequence-detectable chain. A diffuse repertoire
(high entropy) means recognition is degenerate physicochemical or carried by the OTHER chain -> little
signal. We show these predict the observed per-chain NED detectability (held-out ROC) and tie them to the
PDB contact footprint.

    .venv/bin/python bench/holdout_chain_signal.py
"""
from __future__ import annotations

import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import _bench                                                        # noqa: E402
import holdout_eval as HE                                            # noqa: E402
from benchmark import A02, release, vgene                            # noqa: E402

# observed per-chain NED detectability (held-out ROC, baseline step 0) -- the thing we want to predict
NED_ROC = {("NLV", "TRB"): 0.690, ("LLW", "TRB"): 0.511, ("LLL", "TRB"): 0.622, ("ELA", "TRB"): 0.522,
           ("YLQ", "TRB"): 0.954, ("GLC", "TRB"): 0.860,
           ("NLV", "TRA"): 0.674, ("LLL", "TRA"): 0.710, ("GLC", "TRA"): 0.876, ("YLQ", "TRA"): 0.889}
# PDB contact footprint (from tcr-pmhc-structures memory): does this chain's CDR3 apex contact the peptide,
# and is the contact motif-encodable (specific residues) or degenerate (hydrophobic ridge)?
STRUCT = {"NLV": "3gsn/5d2l/5d2n; both chains, central bulge", "YLQ": "7n1f/7n6e/7pbe/7rtr; β apex motif",
          "GLC": "3o4l; hydrophobic ridge (degenerate)", "LLL": "7q9a; hydrophobic ridge (degenerate)",
          "LLW": "no public TCR-pMHC structure", "ELA": "MART-1; α-convergent (TRAV12-2), β diverse"}


def apex_entropy(cdr3s, nbins=5, edge=3):
    """mean per-relative-apex-position Shannon entropy (bits) over the loop tip; low = convergent."""
    bins = [Counter() for _ in range(nbins)]
    for s in cdr3s:
        ap = s[edge:len(s) - edge]
        if len(ap) < 2:
            continue
        for i, c in enumerate(ap):
            bins[min(nbins - 1, i * nbins // len(ap))][c] += 1
    ents = []
    for b in bins:
        tot = sum(b.values())
        if tot:
            ents.append(-sum((n / tot) * math.log2(n / tot) for n in b.values()))
    return float(np.mean(ents)) if ents else float("nan")


def top_motif_coverage(cdr3s, k=3, edge=3, top=10):
    """fraction of binders containing one of the top-`top` apex k-mers; high = convergent."""
    km = Counter()
    for s in cdr3s:
        for i in range(edge, len(s) - k + 1 - edge):
            km[s[i:i + k]] += 1
    topk = {m for m, _ in km.most_common(top)}
    if not topk:
        return float("nan")
    cov = sum(any(s[i:i + k] in topk for i in range(edge, len(s) - k + 1 - edge)) for s in cdr3s)
    return cov / len(cdr3s)


def chain_table(locus):
    rdf = _bench.valid_cdr3(release("vdjdb2026").filter(pl.col("mhc_a").str.contains(A02)
                                                        & (pl.col("gene") == locus)))
    rows = []
    for sh, e in HE.EPI.items():
        sub = rdf.filter(pl.col("epitope") == e).unique("cdr3", maintain_order=True)
        if sub.height < 20:
            continue
        cd = sub["cdr3"].to_list()
        vc = Counter(vgene(x) for x in sub["v"])
        nv = sum(vc.values())
        domV = max(vc.values()) / nv
        vent = -sum((n / nv) * math.log2(n / nv) for n in vc.values())
        rows.append({"epi": sh, "locus": locus, "n": sub.height, "apexH": apex_entropy(cd),
                     "domV": domV, "Vent": vent, "motifcov": top_motif_coverage(cd),
                     "NED": NED_ROC.get((sh, locus), float("nan"))})
    return rows


def spearman(xs, ys):
    xr = np.argsort(np.argsort(xs)).astype(float)
    yr = np.argsort(np.argsort(ys)).astype(float)
    xr -= xr.mean(); yr -= yr.mean()
    d = math.sqrt((xr * xr).sum() * (yr * yr).sum())
    return float((xr * yr).sum() / d) if d else float("nan")


def main():
    rows = chain_table("TRB") + chain_table("TRA")
    print(f"\n=== per-(epitope,chain) convergence vs NED detectability ===")
    print(f"{'epi':5}{'chain':6}{'n':>6}{'apexH':>7}{'domV':>7}{'Vent':>6}{'motifcov':>9}{'NED_ROC':>9}")
    for r in sorted(rows, key=lambda r: -r["NED"] if r["NED"] == r["NED"] else 0):
        print(f"{r['epi']:5}{r['locus']:6}{r['n']:>6}{r['apexH']:>7.2f}{r['domV']:>7.2f}"
              f"{r['Vent']:>6.2f}{r['motifcov']:>9.2f}{r['NED']:>9.3f}")
    valid = [r for r in rows if r["NED"] == r["NED"]]
    y = [r["NED"] for r in valid]
    print("\nSpearman(feature, NED_ROC) across the {} cells:".format(len(valid)))
    for f, sign in [("apexH", "neg"), ("domV", "pos"), ("Vent", "neg"), ("motifcov", "pos")]:
        print(f"  {f:9} ({sign} expected): rho = {spearman([r[f] for r in valid], y):+.3f}")
    print("\nStructural footprint (PDB; from tcr-pmhc-structures):")
    for sh in HE.EPI:
        print(f"  {sh}: {STRUCT.get(sh, '?')}")


if __name__ == "__main__":
    main()
