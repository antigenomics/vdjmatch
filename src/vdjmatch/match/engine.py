"""Matching engine: build seqtree indices over VDJdb (partitioned by gene) and annotate query
clonotypes by fuzzy CDR3 search.

A VDJdb CDR3 can carry several epitope records; we index *unique* CDR3s per gene and expand each
hit back to all of its records with a vectorized polars join (no per-hit Python dicts on the hot
path). Alignment/CIGAR is computed lazily (only for the returned hits, optional).
"""
from __future__ import annotations

import polars as pl
from seqtree import Index, SearchParams

from . import cigar

_DB_COLS = ["cdr3", "v", "j", "epitope", "mhc_a", "mhc_b", "mhc_class",
            "antigen_gene", "antigen_species", "vdjdb_score", "complex_id"]
_HIT_SCHEMA = [("qi", pl.UInt32), ("ref_id", pl.UInt32), ("score", pl.Int32),
               ("n_subs", pl.UInt16), ("n_ins", pl.UInt16), ("n_dels", pl.UInt16)]


def _genefam(col: pl.Expr) -> pl.Expr:
    """Normalize a V/J gene call to its gene name (drop allele ``*NN`` and IMGT decorations)."""
    return col.str.replace(r"\*.*$", "").str.replace(r"/.*$", "")


class VdjdbIndex:
    """Searchable VDJdb partitioned by gene (TRA/TRB). Build once, reuse across samples."""

    def __init__(self, by_gene: dict[str, tuple[Index, pl.DataFrame, pl.DataFrame]]):
        self._by_gene = by_gene  # gene -> (Index, unique_cdr3[ref_id], db_records)

    @classmethod
    def build(cls, df: pl.DataFrame, species: str | None = None) -> "VdjdbIndex":
        if species is not None:
            df = df.filter(pl.col("species") == species)
        df = df.select([c for c in (["gene", *_DB_COLS]) if c in df.columns])
        by_gene = {}
        for gene in df["gene"].unique().to_list() if "gene" in df.columns else [None]:
            part = df.filter(pl.col("gene") == gene) if gene is not None else df
            uc = part.select("cdr3").unique(maintain_order=True).with_row_index("ref_id")
            idx = Index.build(uc["cdr3"].to_list(), alphabet="aa")
            by_gene[gene] = (idx, uc, part)
        return cls(by_gene)

    @property
    def genes(self) -> list[str]:
        return list(self._by_gene)

    def index_for(self, gene: str) -> Index | None:
        g = self._by_gene.get(gene)
        return g[0] if g else None

    def annotate(self, queries: pl.DataFrame, params: SearchParams, *, gene: str,
                 threads: int = 0, match_v: bool = False, match_j: bool = False,
                 align: bool = False, region_aware: bool = False,
                 progress: bool = False, chunk: int = 2000) -> pl.DataFrame:
        """Annotate single-gene query clonotypes; returns a long per-hit frame.

        ``queries`` schema: ``cdr3, v, j, locus, count`` (one gene). Output: query_*, db_*
        (epitope/mhc/...), n_subs/n_ins/n_dels, score, and (if ``align``) cigar/match.

        ``progress`` shows a tqdm bar over the query batch; ``search_batch`` has no callback, so the
        queries are chunked (``chunk`` per step). When ``progress`` is False the whole batch is
        searched in a single call (no behaviour change).
        """
        g = self._by_gene.get(gene)
        if g is None or queries.height == 0:
            return _empty_hits()
        idx, uc, db = g
        q = queries.with_row_index("qi")
        cdr3s = q["cdr3"].to_list()
        if progress:
            from .._util import chunked, progress as _progress
            res = []
            nchunks = -(-len(cdr3s) // chunk) if chunk > 0 else 1
            for c in _progress(chunked(cdr3s, chunk), total=nchunks,
                               desc=f"{gene} search", enable=True):
                res.extend(idx.search_batch(c, params, threads))
        else:
            res = idx.search_batch(cdr3s, params, threads)

        flat = [(qi, h.ref_id, h.score, h.n_subs, h.n_ins, h.n_dels)
                for qi, hl in enumerate(res) for h in hl]
        if not flat:
            return _empty_hits()
        hits = pl.DataFrame(flat, schema=_HIT_SCHEMA, orient="row")
        hits = hits.join(uc, on="ref_id")                                    # +db cdr3
        hits = hits.join(db, on="cdr3")                                      # expand records
        hits = hits.join(q.rename({c: f"q_{c}" for c in q.columns if c != "qi"}), on="qi")

        if match_v:
            hits = hits.filter(_genefam(pl.col("v")) == _genefam(pl.col("q_v")))
        if match_j:
            hits = hits.filter(_genefam(pl.col("j")) == _genefam(pl.col("q_j")))
        if hits.height == 0:
            return _empty_hits()

        out = hits.select(
            "ref_id",
            pl.col("q_cdr3").alias("query_cdr3"), pl.col("q_v").alias("query_v"),
            pl.col("q_j").alias("query_j"), pl.col("q_count").alias("count"),
            pl.col("cdr3").alias("db_cdr3"), pl.col("v").alias("db_v"), pl.col("j").alias("db_j"),
            "epitope", "mhc_a", "mhc_class", "antigen_species", "vdjdb_score", "complex_id",
            "score", "n_subs", "n_ins", "n_dels",
        )
        if align or region_aware:
            ret = pen = None
            if region_aware:
                from . import regions
                ret, pen = regions.load_retention(), regions.vdjam_penalties()
            cig, ml, rs = [], [], []
            qv, qj = out["query_v"].to_list(), out["query_j"].to_list()
            for k, (ref, qc) in enumerate(zip(out["ref_id"], out["query_cdr3"])):
                aln = idx.align(int(ref), qc, params)
                cig.append(cigar.to_cigar(aln.ops))
                ml.append(cigar.match_line(aln.ops))
                if region_aware:
                    rs.append(regions.aligned_score(aln.aligned_query, aln.aligned_ref, aln.ops,
                                                    qv[k], qj[k], gene, ret, pen))
            cols = [pl.Series("cigar", cig), pl.Series("match", ml)]
            if region_aware:
                cols.append(pl.Series("region_score", rs))
            out = out.with_columns(cols)
        return out.drop("ref_id")


def _empty_hits() -> pl.DataFrame:
    str_cols = ["query_cdr3", "query_v", "query_j", "db_cdr3", "db_v", "db_j", "epitope",
                "mhc_a", "mhc_class", "antigen_species"]
    schema = {c: pl.Utf8 for c in str_cols}
    schema |= {"count": pl.Int64, "vdjdb_score": pl.Int64, "complex_id": pl.Int64,
               "score": pl.Int32, "n_subs": pl.UInt16, "n_ins": pl.UInt16, "n_dels": pl.UInt16}
    return pl.DataFrame(schema=schema)
