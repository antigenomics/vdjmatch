"""Real-world annotation evaluation harness for the held-out datasets (extreme-optimization).

Each held-out query is searched ONCE against the candidate reference (full vdjdb2026 A*02, or the >=2-ref
shortlist), its own exact match removed; the per-query hits + control neighbour costs are cached. A
pluggable SCORER then maps (hits, control, ref-stats) -> {epitope: (sim, significant)} so many scoring
formulations are evaluated without re-searching. Per epitope E: argmax assignment, confusion, ROC/PR;
correct/incorrect cases are dumped for study.

    from holdout_eval import build, evaluate, baseline_scorer
    cache = build("TRB")                       # search once (cached to disk)
    evaluate(cache, baseline_scorer)           # apply a scorer -> metrics
"""
from __future__ import annotations

import bisect
import math
import pickle
import sys
from collections import defaultdict
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import _bench                                                        # noqa: E402
from benchmark import (A02, PSSM_SCALE, SOFTV_BETA, _pssm_targets, ref_index,   # noqa: E402
                       release, shortlist, vgene)
from metrics import pr_auc_balanced, roc_auc                         # noqa: E402
from vdjmatch.evalue import background, first_hit                    # noqa: E402
from vdjmatch.match import vgene as _vg                              # noqa: E402

HD = Path.home() / "vcs/manuscripts/2026-vdjmatch/hold_out_data"
OUT = Path.home() / "vcs/manuscripts/2026-vdjmatch/extreme-optimization"
CACHE = OUT / "search_cache"
CACHE.mkdir(parents=True, exist_ok=True)
EPI = {"NLV": "NLVPMVATV", "LLW": "LLWNGPMAV", "LLL": "LLLGIGILV",
       "ELA": "ELAGIGILTV", "YLQ": "YLQPRTFLL", "GLC": "GLCTLVAML"}
A2E = {v: k for k, v in EPI.items()}
DATASETS = {"TRB": ["NLV_TRB", "LLW_TRB", "LLL_TRB", "ELA_TRB", "GLC_TRB", "YLQ_TRB"],
            "TRA": ["NLV_TRA", "LLL_TRA", "GLC_TRA", "YLQ_TRA"]}


def query_set(locus):
    """All held-out queries for a locus: (cdr3, v, true_epi_aa_or_None, dataset). Dedup by cdr3."""
    seen = {}
    for name in DATASETS[locus]:
        f = HD / f"{name}.tsv"
        if not f.exists():
            continue
        epi = EPI[name.split("_")[0]]
        d = pl.read_csv(f, separator="\t", infer_schema_length=0)
        for c, v, lab in zip(d["cdr3"], d["v"], d["label"]):
            seen.setdefault(c, (vgene(v), epi if lab == "1" else None, name))
    return [(c, *t) for c, t in seen.items()]


def build(locus, ref="full", rebuild=False):
    """Search the query set once vs the candidate reference; cache (per-query hits + control costs)."""
    cf = CACHE / f"cache_{locus}_{ref}.pkl"
    if cf.exists() and not rebuild:
        return pickle.loads(cf.read_bytes())
    rdf = release("vdjdb2026").filter(pl.col("mhc_a").str.contains(A02))
    if ref == "shortlist":
        rdf = shortlist(rdf, min_refs=2)
    tgt, ref_epi, ref_v, n_epi, _, _ = ref_index(rdf, locus)
    ctrl = background(locus)
    M, Ntot = len(ctrl), len(ref_epi)
    pool = query_set(locus)
    qs = [r[0] for r in pool]
    params = first_hit.scope(5, 2, 2)
    th, cc = first_hit.scan(tgt, ref_epi, ctrl, qs, target_v=ref_v, params=params, exclude_exact=True)
    pt = _pssm_targets(tgt, ref_epi, ref_v, qs, 5)
    recs = []
    for (c, vq, tr, ds), t, ctr in zip(pool, th, cc):
        hits = [(ps, ed, he, hv) for ps, ed, he, hv in pt[c]]
        first = [(cost, e) for cost, e, *_ in t]
        recs.append({"cdr3": c, "v": vq, "true": tr, "dataset": ds,
                     "hits": hits, "first": first, "ctrl": sorted(ctr)})
    cache = {"locus": locus, "ref": ref, "M": M, "Ntot": Ntot, "n_epi": dict(n_epi), "recs": recs}
    cf.write_bytes(pickle.dumps(cache))
    print(f"[build {locus}/{ref}] {len(recs)} queries, {Ntot} ref CDR3, {len(n_epi)} candidate epitopes",
          file=sys.stderr)
    return cache


# -------- scorers: (rec, M, n_epi) -> {epitope: (sim, significant)} ----------------------------------
def baseline_scorer(rec, M, n_epi, alpha=1e-3):
    """size-invariant control-calibrated similarity for argmax + first-hit V-agnostic significance S."""
    cs = rec["ctrl"]
    Ntot = sum(n_epi.values())
    t1 = [(c, e) for c, e in rec["first"] if c <= 1]
    significant = first_hit.pvalue(t1, cs, Ntot, M)["p_enrichment"] < alpha
    sim = defaultdict(float)
    for ps, ed, he, hv in rec["hits"]:
        w = 1.0 if _vg.gene_family(rec["v"]) == _vg.gene_family(hv) else SOFTV_BETA * _vg.vsim(rec["v"], hv)
        if w <= 0:
            continue
        nc = bisect.bisect_right(cs, ed) + 1
        sim[he] += w * math.exp(-ps / PSSM_SCALE) / nc
    return {e: (s, significant) for e, s in sim.items()}


def evaluate(cache, scorer, dump=False, **kw):
    M, n_epi, recs = cache["M"], cache["n_epi"], cache["recs"]
    per = [scorer(r, M, n_epi, **kw) for r in recs]
    assigns, sims = [], []
    for sc in per:
        cand = {e: s for e, (s, sg) in sc.items() if sg}
        a = max(cand, key=cand.get) if cand else None
        assigns.append(a)
        sims.append({e: s for e, (s, sg) in sc.items()})
    tested = sorted({r["true"] for r in recs if r["true"]}, key=lambda e: A2E.get(e, e))
    print(f"\n=== {cache['locus']} / {cache['ref']} ===  n={len(recs)}")
    print(f"{'epi':5}{'n+':>5}{'TP':>5}{'FN':>5}{'FP':>4}{'prec':>7}{'rec':>7}{'F1':>7}{'ROC':>7}{'balPR':>7}")
    rows = []
    for e in tested:
        sh = A2E[e]
        pairs = [(int(r["true"] == e), sims[i].get(e, 0.0)) for i, r in enumerate(recs)]
        tp = sum(1 for i, r in enumerate(recs) if r["true"] == e and assigns[i] == e)
        fn = sum(1 for i, r in enumerate(recs) if r["true"] == e and assigns[i] != e)
        fp = sum(1 for i, r in enumerate(recs) if r["true"] != e and assigns[i] == e)
        npos = tp + fn
        prec = tp / (tp + fp) if tp + fp else float("nan")
        rec_ = tp / npos if npos else float("nan")
        f1 = 2 * prec * rec_ / (prec + rec_) if (prec + rec_) and prec == prec and rec_ == rec_ else float("nan")
        roc, bp = roc_auc(pairs), pr_auc_balanced(pairs)
        rows.append((sh, npos, roc, bp, prec, rec_, f1))
        print(f"{sh:5}{npos:>5}{tp:>5}{fn:>5}{fp:>4}{prec:>7.3f}{rec_:>7.3f}{f1:>7.3f}{roc:>7.3f}{bp:>7.3f}")
    mroc = sum(r[2] for r in rows) / len(rows)
    print(f"  mean ROC {mroc:.3f}  mean balPR {sum(r[3] for r in rows)/len(rows):.3f}")
    return rows


if __name__ == "__main__":
    for loc in ("TRB", "TRA"):
        evaluate(build(loc, "full"), baseline_scorer)
