"""Real post-selection control repertoire (isalgo/airr_control) as the prior/null background --- the
TCRNET/mirpy approach, instead of OLGA (pre-selection) or other-epitope binders.

The HF dataset ships VDJtools .aa tables (count, freq, cdr3nt, cdr3aa, v, d, j, ...). The .aa aggregation
is already productive: 0/18M rows carry a stop/out-of-frame marker. We still apply `_bench.valid_cdr3`
(C-anchor … F/W) for safety, and RANDOM-sample (seed 42) rather than take the count-ordered head, so the
composition background is not biased toward expanded clones.

    from hardcase_control import control_cdr3
    cd = control_cdr3("TRB")     # ~100k unique productive post-selection CDR3 (cached)
"""
from __future__ import annotations

import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _bench                                                        # noqa: E402
import benchmark as B                                               # noqa: E402

CACHE = Path(__file__).resolve().parent / "_hardcase"
CACHE.mkdir(exist_ok=True)
_MEM: dict = {}


def control_frame(locus: str = "TRB", n: int = 100000, seed: int = 42) -> pl.DataFrame:
    """Unique productive post-selection (cdr3, v) from isalgo/airr_control, random-sampled to ``n``."""
    cf = CACHE / f"control_{locus}_{n}_{seed}.parquet"
    if cf.exists():
        return pl.read_parquet(cf)
    from huggingface_hub import hf_hub_download
    f = hf_hub_download("isalgo/airr_control", repo_type="dataset",
                        filename=f"human.{locus.lower()}.aa.vdjtools.tsv.gz")
    d = pl.read_csv(f, separator="\t", columns=["cdr3aa", "v"], infer_schema_length=0)
    raw = d.height
    d = (d.rename({"cdr3aa": "cdr3"}).pipe(_bench.valid_cdr3)         # productive: C…F/W, standard AA
         .filter(~pl.col("cdr3").str.contains(r"[*_~]"))
         .unique("cdr3", maintain_order=True))
    d = d.with_columns(v=pl.col("v").map_elements(B.vgene, return_dtype=pl.Utf8))
    if d.height > n:
        d = d.sort("cdr3").sample(n, seed=seed)                       # deterministic unbiased sample
    print(f"[control {locus}] {raw:,} raw -> {d.height:,} productive unique (sampled, seed {seed})",
          file=sys.stderr)
    d.write_parquet(cf)
    return d


def control_cdr3(locus: str = "TRB", n: int = 100000) -> list[str]:
    if (locus, n) not in _MEM:
        _MEM[(locus, n)] = control_frame(locus, n)["cdr3"].to_list()
    return _MEM[(locus, n)]


if __name__ == "__main__":
    for loc in ("TRB", "TRA"):
        d = control_frame(loc)
        L = d["cdr3"].str.len_chars()
        print(f"{loc}: {d.height:,} unique productive control CDR3 | "
              f"length median={L.median():.0f} mean={L.mean():.1f} "
              f"| top V {d['v'].value_counts(sort=True)['v'][0]}")
