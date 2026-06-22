"""High-level annotation API.

One ergonomic entry point over the VDJdb reference + ``seqtree`` search engine, from the simplest
(``list[CDR3] -> hits``) to the fully general (``polars df -> df + annotation columns``, ids/labels
preserved), single- or paired-chain, against any VDJdb version or a custom reference.

    import vdjmatch
    ann = vdjmatch.Annotator.latest()                       # pinned HF release (or .version(tag))
    ann.hits(["CASSIRSSYEQYF", "CASSLAPGATNEKLFF"])          # list  -> long per-hit polars frame
    ann.annotate(df, cdr3="junction_aa", locus="locus")     # df    -> df + vdjmatch_* columns
    ann.annotate_paired(cell_df, cdr3a="cdr3_alpha_aa", cdr3b="cdr3_beta_aa")   # paired alpha+beta

    vdjmatch.annotate(["CASSIRSSYEQYF", ...])                # module-level shortcut (cached default ref)
"""
from __future__ import annotations

import functools

import polars as pl

from . import db
from .match.engine import VdjdbIndex
from .match.scope import search_params

DEFAULT_SCOPE = "1,0,0,1"          # Hamming-1: the signal:noise optimum (see appendix)
_LOCUS = {"TRA": "TRA", "TRB": "TRB", "A": "TRA", "B": "TRB",
          "alpha": "TRA", "beta": "TRB", "TCRA": "TRA", "TCRB": "TRB"}


def _norm_locus(x) -> str | None:
    return _LOCUS.get(str(x)) or _LOCUS.get(str(x).upper())


def _summarise(hits: pl.DataFrame, prefix: str) -> pl.DataFrame:
    """Long per-hit frame -> one row per query CDR3: top epitope (distance-weighted vote), its score,
    total hit count, and the matched MHC of the top epitope."""
    cols = ["query_cdr3", f"{prefix}epitope", f"{prefix}mhc_class", f"{prefix}score", f"{prefix}n_hits"]
    if hits.height == 0:
        return pl.DataFrame(schema={c: (pl.Utf8 if "epitope" in c or "mhc" in c or c == "query_cdr3"
                                        else pl.Float64) for c in cols})
    w = hits.with_columns((1.0 / (1 + pl.col("n_subs"))).alias("_w"))
    per = (w.group_by(["query_cdr3", "epitope", "mhc_class"])
             .agg(pl.col("_w").sum().alias("_ws"))
             .sort(["query_cdr3", "_ws"], descending=[False, True])
             .unique(subset="query_cdr3", keep="first", maintain_order=True))
    nh = hits.group_by("query_cdr3").len().rename({"len": f"{prefix}n_hits"})
    return (per.join(nh, on="query_cdr3")
               .select("query_cdr3", pl.col("epitope").alias(f"{prefix}epitope"),
                       pl.col("mhc_class").alias(f"{prefix}mhc_class"),
                       pl.col("_ws").alias(f"{prefix}score"), f"{prefix}n_hits"))


class Annotator:
    """A reusable, indexed VDJdb reference. Build once, annotate many query sets."""

    def __init__(self, index: VdjdbIndex):
        self._index = index

    # ---- construction (any VDJdb version / custom) ----
    @classmethod
    def latest(cls, *, species: str | None = "HomoSapiens", source: str = "hf") -> "Annotator":
        """Pinned HF benchmark release (``source="hf"``) or the latest GitHub release."""
        path = db.fetch_hf() if source == "hf" else db.fetch_latest(asset="default")
        return cls.from_frame(db.load(path), species=species)

    @classmethod
    def version(cls, tag: str, *, species: str | None = "HomoSapiens") -> "Annotator":
        return cls.from_frame(db.load(db.fetch_hf(tag=tag)), species=species)

    @classmethod
    def from_path(cls, path, *, species: str | None = None) -> "Annotator":
        return cls.from_frame(db.load(path), species=species)

    @classmethod
    def from_frame(cls, vdj: pl.DataFrame, *, species: str | None = None) -> "Annotator":
        """Custom VDJdb-like frame (normalized columns: ``gene,cdr3,v,j,epitope,mhc_*,...``)."""
        return cls(VdjdbIndex.build(vdj, species=species))

    @property
    def loci(self) -> list[str]:
        return self._index.genes

    # ---- queries ----
    def hits(self, cdr3s, *, locus: str = "TRB", scope: str = DEFAULT_SCOPE,
             match_v: bool = False, match_j: bool = False, **kw) -> pl.DataFrame:
        """``list[CDR3]`` (or one string) -> long per-hit frame (query/db CDR3, epitope, score, edits)."""
        if isinstance(cdr3s, str):
            cdr3s = [cdr3s]
        q = pl.DataFrame({"cdr3": list(cdr3s)}).with_columns(
            v=pl.lit(None, pl.Utf8), j=pl.lit(None, pl.Utf8), locus=pl.lit(locus), count=pl.lit(1))
        sp = scope if not isinstance(scope, str) else search_params(scope)
        return self._index.annotate(q, sp, gene=_norm_locus(locus), match_v=match_v, match_j=match_j, **kw)

    def annotate(self, data, *, cdr3: str = "cdr3", v: str | None = None, j: str | None = None,
                 locus: str | None = None, scope: str = DEFAULT_SCOPE, prefix: str = "vdjmatch_",
                 match_v: bool = False, match_j: bool = False) -> pl.DataFrame:
        """Annotate a ``list[CDR3]`` or a polars frame, returning the input with appended
        ``{prefix}epitope/mhc_class/score/n_hits`` columns. All input rows/columns/ids/labels are kept
        (we annotate unique CDR3s per locus and join back). ``locus`` may name a column (per-row locus)
        or be a single locus string (default TRB)."""
        if not isinstance(data, pl.DataFrame):
            data = pl.DataFrame({cdr3: list(data)})
        tmp = None
        if locus is not None and locus in data.columns:
            loc_col = locus
        else:
            loc_col = tmp = "_locus"
            data = data.with_columns(pl.lit(locus or "TRB").alias(tmp))
        sp = scope if not isinstance(scope, str) else search_params(scope)
        sub_summaries = []
        for lv in data[loc_col].unique().to_list():
            gene = _norm_locus(lv)
            if gene not in self.loci:
                continue
            part = data.filter(pl.col(loc_col) == lv)
            sel = {cdr3: "cdr3"} | ({v: "v"} if v else {}) | ({j: "j"} if j else {})
            uq = part.select(list(sel)).unique().rename(sel)
            for c in ("v", "j"):
                if c not in uq.columns:
                    uq = uq.with_columns(pl.lit(None, pl.Utf8).alias(c))
            uq = uq.with_columns(count=pl.lit(1), locus=pl.lit(gene))
            h = self._index.annotate(uq, sp, gene=gene, match_v=match_v, match_j=match_j)
            s = _summarise(h, prefix).rename({"query_cdr3": cdr3}).with_columns(pl.lit(lv).alias(loc_col))
            sub_summaries.append(s)
        out = data
        if sub_summaries:
            out = out.join(pl.concat(sub_summaries), on=[cdr3, loc_col], how="left")
        return out.drop(tmp) if tmp else out

    def annotate_paired(self, data: pl.DataFrame, *, cdr3a: str = "cdr3_alpha_aa",
                        cdr3b: str = "cdr3_beta_aa", scope: str = DEFAULT_SCOPE,
                        prefix: str = "vdjmatch_") -> pl.DataFrame:
        """Paired alpha+beta: annotate each chain, then call the epitope supported by **both** chains
        (intersection; score = sum of the two chain scores). Rows/ids/labels preserved. (The
        control-calibrated paired *E-value* lives in ``match.PairedVdjdbIndex.annotate_pairs``.)"""
        a = self.annotate(data, cdr3=cdr3a, locus="TRA", scope=scope, prefix="_a_")
        ab = self.annotate(a, cdr3=cdr3b, locus="TRB", scope=scope, prefix="_b_")
        agree = (pl.when(pl.col("_a_epitope") == pl.col("_b_epitope"))
                   .then(pl.col("_a_epitope")).otherwise(None))
        return ab.with_columns(
            agree.alias(f"{prefix}epitope"),
            (pl.col("_a_score").fill_null(0) + pl.col("_b_score").fill_null(0)).alias(f"{prefix}score"),
            (pl.col("_a_n_hits").fill_null(0) + pl.col("_b_n_hits").fill_null(0)).alias(f"{prefix}n_hits"),
        ).drop([c for c in ab.columns if c.startswith("_a_") or c.startswith("_b_")])


@functools.cache
def _default() -> Annotator:
    return Annotator.latest()


def annotate(data, **kw) -> pl.DataFrame:
    """Module-level shortcut: annotate against the default (latest pinned) reference. See
    :meth:`Annotator.annotate`."""
    return _default().annotate(data, **kw)
