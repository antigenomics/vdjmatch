"""Internal clustering-quality metrics for the vdjmatch TCR clustering, alongside purity/retention.

External (label-based): purity, retention (from cluster_trim). Internal (geometry of the clustering in the
clustering's OWN distance = CDR3 edit distance): Hopkins score HS (cluster tendency: ~0.5 random, ->1 strongly
clustered), CTS (normalised cluster-tendency = 2|HS-0.5|), Silhouette (precomputed edit distance; cohesion vs
separation of the vdjmatch clusters), Calinski-Harabasz (between/within dispersion ratio on a classical-MDS
embedding of the edit distance). Silhouette/CH use the vdjmatch cluster labels over non-singleton clonotypes;
HS on the MDS embedding. Computed per (reference, chain, trim), sub-sampled (seed 0) to <=SAMPLE.

    .venv/bin/python bench/cluster_metrics.py [TRA|TRB]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path.home() / "vcs/manuscripts/2026-vdjmatch/benchmarks/scripts"))
import benchmark as B                                                # noqa: E402
import cluster_results as CR                                         # noqa: E402
import cluster_trim as CT                                            # noqa: E402
from rapidfuzz.distance import Levenshtein                           # noqa: E402
from rapidfuzz.process import cdist                                  # noqa: E402
from sklearn.metrics import calinski_harabasz_score, silhouette_score   # noqa: E402

SAMPLE = 1500


def edist(cdr3s):
    """n x n CDR3 edit-distance matrix (rapidfuzz, C)."""
    return cdist(cdr3s, cdr3s, scorer=Levenshtein.distance, dtype=np.float64)


def cmds(D, ndim=20):
    """classical MDS (PCoA): Euclidean coords approximating the distance matrix D."""
    n = len(D)
    J = np.eye(n) - 1.0 / n
    Bm = -0.5 * J @ (D ** 2) @ J
    w, V = np.linalg.eigh(Bm)
    o = np.argsort(w)[::-1][:ndim]
    return V[:, o] * np.sqrt(np.clip(w[o], 0, None))


def hopkins(X, seed=0, frac=0.1):
    """Hopkins statistic on a Euclidean embedding: ~0.5 random, ->1 clustered."""
    from sklearn.neighbors import NearestNeighbors
    rng = np.random.default_rng(seed)
    n, d = X.shape
    m = max(20, int(frac * n))
    nn = NearestNeighbors(n_neighbors=2).fit(X)
    real = nn.kneighbors(X[rng.choice(n, m, replace=False)])[0][:, 1]
    rand = nn.kneighbors(rng.uniform(X.min(0), X.max(0), (m, d)))[0][:, 0]
    s = float(rand.sum() + real.sum())
    return float(rand.sum() / s) if s else 0.5


def metrics(cdr3, v, epi, trim, seed=0):
    lab = np.array(CT.cluster_labels(cdr3, v, epi, trim))
    n = len(cdr3)
    rng = np.random.default_rng(seed)
    sel = rng.choice(n, SAMPLE, replace=False) if n > SAMPLE else np.arange(n)
    D = edist([cdr3[i] for i in sel])
    X = cmds(D)
    hs = hopkins(X, seed)
    ls = lab[sel]
    keep = ls >= 0                                                       # non-singleton clonotypes
    sil = ch = float("nan")
    if keep.sum() >= 3 and len(set(ls[keep])) >= 2:
        sil = float(silhouette_score(D[np.ix_(keep, keep)], ls[keep], metric="precomputed"))
        ch = float(calinski_harabasz_score(X[keep], ls[keep]))
    pur, ret, nc, _ = CT.cluster_pr(cdr3, v, epi, trim)
    return dict(n=n, purity=pur, retention=ret, HS=hs, CTS=2 * abs(hs - 0.5), silhouette=sil, CH=ch)


def main(which):
    d = B.release("vdjdb2026"); sl = B.shortlist(d)
    print(f"{'ref':10}{'chain':6}{'trim':6}{'n':>7}{'purity':>8}{'retent':>8}{'HS':>7}{'CTS':>7}{'silh':>7}{'CH':>9}")
    for ref_name, df in (("shortlist", sl), ("full", d)):
        cdr3, epi, v = CR.single_clonotypes(df, which)
        for trim in (False, True):
            m = metrics(cdr3, v, epi, trim)
            print(f"{ref_name:10}{which:6}{'apex' if trim else 'none':6}{m['n']:>7}{m['purity']:>8.3f}"
                  f"{m['retention']:>8.3f}{m['HS']:>7.3f}{m['CTS']:>7.3f}{m['silhouette']:>7.3f}{m['CH']:>9.1f}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "TRB")
