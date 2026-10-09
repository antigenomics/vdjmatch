"""One native batch per locus, with identity-preserving reference expansion."""

from __future__ import annotations

import polars as pl
from seqtree import Index, SearchParams

from ..db import schema
from . import cigar

_HIT_SCHEMA = [
    ("qi", pl.UInt32),
    ("ref_id", pl.UInt32),
    ("score", pl.Int32),
    ("n_subs", pl.UInt16),
    ("n_ins", pl.UInt16),
    ("n_dels", pl.UInt16),
]


def _genefam(col: pl.Expr) -> pl.Expr:
    return col.str.replace(r"\*.*$", "").str.replace(r"/.*$", "")


class VdjdbIndex:
    """Explicitly owned reference index; unique junction search keys, independent records retained."""

    def __init__(self, by_gene, reference=None):
        self._by_gene = by_gene
        self.reference = reference

    @classmethod
    def build(cls, df: pl.DataFrame, species: str | None = None) -> "VdjdbIndex":
        df = schema.normalize(df).with_columns(
            pl.col("vdjdb_score").cast(pl.Int64, strict=False).fill_null(0)
        )
        if species is not None:
            if df.height and df["species"].null_count() == df.height:
                raise ValueError(
                    "reference has no species metadata; use species=None explicitly"
                )
            df = df.filter(pl.col("species") == species)
        all_records = df
        df = df.filter(pl.col("reference_valid"))
        # Observation IDs stay distinct from sequence search keys, including custom frames.
        if "record_id" not in df.columns:
            df = df.with_row_index("record_id").with_columns(
                pl.col("record_id").cast(pl.String)
            )
        by_gene = {}
        for gene in sorted(df["gene"].drop_nulls().unique().to_list()):
            part = df.filter(pl.col("gene") == gene)
            uc = (
                part.select("cdr3").unique(maintain_order=True).with_row_index("ref_id")
            )
            by_gene[gene] = (Index.build(uc["cdr3"].to_list(), alphabet="aa"), uc, part)
        return cls(by_gene, all_records)

    @property
    def genes(self) -> list[str]:
        return list(self._by_gene)

    def index_for(self, gene: str) -> Index | None:
        g = self._by_gene.get(gene)
        return g[0] if g else None

    def records_for(self, gene: str) -> pl.DataFrame:
        """Read-only reference observations for candidate coverage/provenance."""
        g = self._by_gene.get(gene)
        return (
            g[2].clone()
            if g
            else schema.normalize(
                pl.DataFrame(
                    schema={"gene": pl.String, "cdr3": pl.String, "epitope": pl.String}
                )
            )
        )

    def empty_hits(self, query_id_type=pl.UInt32, *, align=False, region_aware=False):
        """Typed empty evidence with the same reference metadata columns as ordinary hits."""
        out = _empty_hits(query_id_type)
        if self._by_gene:
            records = next(iter(self._by_gene.values()))[2]
            mapping = {"cdr3": "db_cdr3", "v": "db_v", "j": "db_j"}
            for c, dtype in records.schema.items():
                name = mapping.get(c, c)
                if name in {
                    "query_id",
                    "query_cdr3",
                    "query_v",
                    "query_j",
                    "query_locus",
                    "count",
                    "score",
                    "n_subs",
                    "n_ins",
                    "n_dels",
                    "qi",
                    "ref_id",
                }:
                    name = "db_metadata_" + c
                out = out.with_columns(pl.lit(None, dtype=dtype).alias(name))
        if align or region_aware:
            out = out.with_columns(
                pl.lit(None, dtype=pl.String).alias("cigar"),
                pl.lit(None, dtype=pl.String).alias("match"),
            )
        if region_aware:
            out = out.with_columns(pl.lit(None, dtype=pl.Float64).alias("region_score"))
        return out

    def annotate(
        self,
        queries: pl.DataFrame,
        params: SearchParams,
        *,
        gene: str,
        threads: int = 1,
        match_v: bool = False,
        match_j: bool = False,
        align: bool = False,
        region_aware: bool = False,
        progress: bool = False,
        chunk: int = 2000,
    ) -> pl.DataFrame:
        """Long hit evidence, preserving query_id/locus and all reference metadata.

        Search is one native batch even with progress enabled. ``threads`` controls seqtree;
        Polars performs sequential post-search stages in its existing process-wide pool.
        ``chunk`` is retained for signature compatibility; it does not split native search.
        """
        if isinstance(threads, bool) or not isinstance(threads, int) or threads < 0:
            raise ValueError(
                "threads must be a nonnegative integer (0 selects native automatic)"
            )
        q = queries
        if "query_id" not in q.columns:
            q = q.with_row_index("query_id")
        if q["query_id"].null_count() or q["query_id"].n_unique() != q.height:
            raise ValueError("query_id must be unique and non-null within a batch")
        for c, dtype, value in [
            ("v", pl.String, None),
            ("j", pl.String, None),
            ("count", pl.Int64, 1),
            ("locus", pl.String, gene),
            ("sequence_id", pl.String, None),
            ("pair_id", pl.String, None),
        ]:
            if c not in q.columns:
                q = q.with_columns(pl.lit(value, dtype=dtype).alias(c))
        g = self._by_gene.get(gene)
        if g is None or q.height == 0:
            return self.empty_hits(
                q.schema["query_id"], align=align, region_aware=region_aware
            )
        idx, uc, records = g
        q = q.filter(
            pl.col("cdr3").is_not_null()
            & pl.col("cdr3").str.contains(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
        ).with_row_index("qi")
        if q.height == 0:
            return self.empty_hits(
                q.schema["query_id"], align=align, region_aware=region_aware
            )
        if progress:
            import sys

            print(
                f"{gene}: searching {q.height:,} queries in one native batch",
                file=sys.stderr,
            )
        res = idx.search_batch(q["cdr3"].to_list(), params, threads)
        if progress:
            print(f"{gene}: native search complete", file=sys.stderr)
        flat = [
            (qi, h.ref_id, h.score, h.n_subs, h.n_ins, h.n_dels)
            for qi, hl in enumerate(res)
            for h in hl
        ]
        if not flat:
            return self.empty_hits(
                q.schema["query_id"], align=align, region_aware=region_aware
            )
        hits = pl.DataFrame(flat, schema=_HIT_SCHEMA, orient="row").join(
            uc, on="ref_id"
        )
        hits = hits.join(
            q.select(
                "qi",
                "query_id",
                pl.col("cdr3").alias("query_cdr3"),
                pl.col("v").alias("query_v"),
                pl.col("j").alias("query_j"),
                pl.col("locus").alias("query_locus"),
                "count",
                *[
                    pl.col(c).alias("query_" + c)
                    for c in ("sequence_id", "pair_id")
                    if c in q.columns
                ],
            ),
            on="qi",
            validate="m:1",
        )
        # Align sequence pairs before many-record expansion. No repeated work per observation.
        if align or region_aware:
            pairs = hits.select("ref_id", "query_cdr3", "query_v", "query_j").unique(
                maintain_order=True
            )
            ret = pen = None
            if region_aware:
                from . import regions

                ret, pen = regions.load_retention(), regions.vdjam_penalties()
            cs, ml, rs = [], [], []
            for ref, qc, qv, qj in pairs.iter_rows():
                aln = idx.align(int(ref), qc, params)
                cs.append(cigar.to_cigar(aln.ops))
                ml.append(cigar.match_line(aln.ops))
                if region_aware:
                    rs.append(
                        regions.aligned_score(
                            aln.aligned_query,
                            aln.aligned_ref,
                            aln.ops,
                            qv,
                            qj,
                            gene,
                            ret,
                            pen,
                        )
                    )
            pairs = pairs.with_columns(pl.Series("cigar", cs), pl.Series("match", ml))
            if region_aware:
                pairs = pairs.with_columns(pl.Series("region_score", rs))
            hits = hits.join(
                pairs,
                on=["ref_id", "query_cdr3", "query_v", "query_j"],
                nulls_equal=True,
                validate="m:1",
            )
        shared = (set(hits.columns) & set(records.columns)) - {"cdr3"}
        records = records.rename({c: f"db_metadata_{c}" for c in shared})
        hits = hits.join(records, on="cdr3").rename(
            {"cdr3": "db_cdr3", "v": "db_v", "j": "db_j"}
        )
        if match_v:
            hits = hits.filter(_genefam(pl.col("db_v")) == _genefam(pl.col("query_v")))
        if match_j:
            hits = hits.filter(_genefam(pl.col("db_j")) == _genefam(pl.col("query_j")))
        out = hits.drop("qi", "ref_id").sort(["query_id", "score", "db_cdr3"])
        template = self.empty_hits(
            q.schema["query_id"], align=align, region_aware=region_aware
        )
        return out.select(
            *[pl.col(c).cast(dtype) for c, dtype in template.schema.items()],
            *[pl.col(c) for c in out.columns if c not in template.columns],
        )


def _empty_hits(query_id_type=pl.UInt32) -> pl.DataFrame:
    strings = [
        "query_cdr3",
        "query_v",
        "query_j",
        "query_locus",
        "query_sequence_id",
        "query_pair_id",
        "db_cdr3",
        "db_v",
        "db_j",
        "epitope",
        "mhc_a",
        "mhc_b",
        "mhc_class",
        "antigen_gene",
        "antigen_species",
        "reference_id",
        "record_id",
        "gene",
        "species",
    ]
    return pl.DataFrame(
        schema={
            **{c: pl.String for c in strings},
            "query_id": query_id_type,
            "count": pl.Int64,
            "vdjdb_score": pl.Int64,
            "complex_id": pl.Int64,
            "score": pl.Int32,
            "n_subs": pl.UInt16,
            "n_ins": pl.UInt16,
            "n_dels": pl.UInt16,
        }
    )
