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
    return col.str.replace(r"\*.*$", "")


def _search_unique_junctions(index, junctions, params, threads):
    """Search each distinct junction once; expand read-only hit lists to query order.

    Query V/J calls, counts and pairing affect later evidence, not sequence search.
    Shared inner hit lists must not be mutated by callers. Nothing persists across calls.
    """
    keys = list(dict.fromkeys(junctions))
    results = index.search_batch(keys, params, threads)
    by_junction = dict(zip(keys, results))
    return [by_junction[junction] for junction in junctions]


class VdjdbIndex:
    """Explicitly owned reference index; unique junction search keys, independent records retained."""

    def __init__(self, by_gene, reference=None):
        self._by_gene = by_gene
        self.reference = reference

    @classmethod
    def build(cls, df: pl.DataFrame, species: str | None = None, *, reference_population='curated-benchmark') -> "VdjdbIndex":
        """Build the declared fixed population; all-valid bypasses source curation only."""
        df = schema.normalize(df).with_columns(
            pl.col("vdjdb_score").cast(pl.Int64, strict=False).fill_null(0)
        )
        if species is not None:
            if df.height and df["species"].null_count() == df.height:
                raise ValueError(
                    "reference has no species metadata; use species=None explicitly"
                )
            df = df.filter(pl.col("species") == species)
        if reference_population not in ('curated-benchmark', 'all-valid'):
            raise ValueError('unknown reference population')
        if reference_population == 'all-valid':
            df = df.with_columns(pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False).alias('reference_valid'))
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
        exclude_exact: bool = False,
        align: bool = False,
        region_aware: bool = False,
        progress: bool = False,
        chunk: int = 2000,
        consumer=None,
        _reference_records=None,
        _candidate_only=False,
    ) -> pl.DataFrame:
        """Long hit evidence, preserving query_id/locus and all reference metadata.

        Search is one native batch of distinct junctions even with progress enabled. ``threads`` controls seqtree;
        Polars performs sequential post-search stages in its existing process-wide pool.
        With ``consumer``, one native population batch emits owned primitive buffers
        for at most ``chunk`` complete distinct junctions per block. Python Hit
        objects and the full hit table are never retained. Duplicate junctions
        expand back to every original query identity in complete-query blocks
        of at most ``chunk`` original identities.
        ``exclude_exact`` removes zero-edit hits before alignment and observation expansion.
        """
        if isinstance(threads, bool) or not isinstance(threads, int) or threads < 0:
            raise ValueError(
                "threads must be a nonnegative integer (0 selects native automatic)"
            )
        if consumer is not None and (type(chunk) is not int or chunk < 1):
            raise ValueError('positive integer output chunk required')
        if _candidate_only and consumer is None:
            raise ValueError('candidate projection requires a consumer')
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
        if consumer is not None and not callable(getattr(idx, 'search_batch_consume', None)):
            raise RuntimeError('bounded hit consumption requires seqtree Index.search_batch_consume; rebuild/install the current native source')
        if _reference_records is not None:
            if not records.select('record_id','cdr3').equals(_reference_records.select('record_id','cdr3')):
                raise ValueError('gene annotation overrides must retain every reference identity and junction')
            records = _reference_records
        if _candidate_only:
            # Keep every source row, including duplicates; only unused metadata
            # columns are omitted before the many-observation hit expansion.
            records = records.select('cdr3', 'v', 'j', 'epitope', 'mhc_a', 'mhc_b',
                'mhc_class', 'reference_id', 'vdjdb_score', 'antigen_species')
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
        query_junctions = q["cdr3"].to_list()
        reference_junctions = uc["cdr3"].to_list() if exclude_exact else None
        res = (_search_unique_junctions(idx, query_junctions, params, threads)
            if consumer is None else None)
        if progress and consumer is None:
            print(f"{gene}: native search complete", file=sys.stderr)
        def expanded(flat):
            if isinstance(flat, pl.DataFrame):
                hits = flat
            elif flat:
                hits = pl.DataFrame(flat, schema=_HIT_SCHEMA, orient="row")
            else:
                return self.empty_hits(
                    q.schema["query_id"], align=align, region_aware=region_aware
                )
            hits = hits.join(uc, on="ref_id")
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
            if exclude_exact and consumer is not None:
                hits = hits.filter(pl.col('query_cdr3') != pl.col('cdr3'))
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
            record_rows = records.rename({c: f"db_metadata_{c}" for c in shared})
            hits = hits.join(record_rows, on="cdr3").rename(
                {"cdr3": "db_cdr3", "v": "db_v", "j": "db_j"}
            )
            if match_v:
                hits = hits.filter(_genefam(pl.col("db_v")) == _genefam(pl.col("query_v")))
            if match_j:
                hits = hits.filter(_genefam(pl.col("db_j")) == _genefam(pl.col("query_j")))
            out = hits.drop("qi", "ref_id").sort(["query_id", "score", "db_cdr3"])
            if _candidate_only:
                return out
            template = self.empty_hits(
                q.schema["query_id"], align=align, region_aware=region_aware
            )
            return out.select(
                *[pl.col(c).cast(dtype) for c, dtype in template.schema.items()],
                *[pl.col(c) for c in out.columns if c not in template.columns],
            )
        if consumer is not None:
            import numpy as np
            unique = q.select('cdr3').unique(maintain_order=True).with_row_index('_native_qi')
            mapping = q.select('qi', 'cdr3').join(unique, on='cdr3', validate='m:1').select('_native_qi', 'qi').sort('_native_qi', 'qi')
            offsets = np.r_[0, np.cumsum(np.bincount(mapping['_native_qi'].to_numpy(), minlength=unique.height))]
            primitive_schema = [('_native_qi', pl.UInt32), *_HIT_SCHEMA[1:]]
            def emit(start, stop, buffer):
                hits = pl.DataFrame(np.asarray(buffer), schema=primitive_schema, orient='row')
                members = mapping.slice(int(offsets[start]), int(offsets[stop]-offsets[start]))
                for part in members.iter_slices(chunk):
                    selected = hits.join(part, on='_native_qi', validate='m:m').drop('_native_qi')
                    consumer(expanded(selected))
            idx.search_batch_consume(unique['cdr3'].to_list(), params, emit,
                threads=threads, chunk_size=chunk, max_hit_bytes=64*1024**2)
            if progress:
                print(f"{gene}: native search complete", file=sys.stderr)
            return self.empty_hits(q.schema['query_id'], align=align, region_aware=region_aware)
        flat = [
            (qi, h.ref_id, h.score, h.n_subs, h.n_ins, h.n_dels)
            for qi, hl in enumerate(res)
            for h in hl
            if not exclude_exact or query_junctions[qi] != reference_junctions[h.ref_id]
        ]
        return expanded(flat)



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
