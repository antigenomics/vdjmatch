"""Background control repertoires for E-value calibration (thin wrapper over seqtree)."""

from __future__ import annotations

from seqtree import Index
from seqtree.control import load_control

# query locus -> seqtree control name (bundled human_trb_aa; others via HuggingFace)
_CONTROL = {
    "human:TRB": "human_trb_aa",
    "human:TRA": "human_tra_aa",
    "mouse:TRB": "mouse_trb_aa",
    "mouse:TRA": "mouse_tra_aa",
}


def _organism(species: str) -> str:
    aliases = {
        "human": "human",
        "homosapiens": "human",
        "homo sapiens": "human",
        "mouse": "mouse",
        "musmusculus": "mouse",
        "mus musculus": "mouse",
    }
    organism = aliases.get(str(species).strip().lower())
    if organism is None:
        raise ValueError(f"no control for species={species!r}")
    return organism


def background(
    locus: str = "TRB",
    species: str = "human",
    size: int | None = None,
    cache_dir: str | None = None,
) -> Index:
    """Load a deduplicated background repertoire ``Index`` for the given locus/species.

    Bundled: human TRB. Others (human TRA, mouse TRA/TRB) download via ``seqtree[control]``.
    """
    organism = _organism(species)
    name = _CONTROL.get(f"{organism}:{locus}")
    if name is None:
        raise ValueError(f"no control for locus={locus!r} species={species!r}")
    return load_control(name, size=size, cache_dir=cache_dir)


def raw_background(locus="TRB", species="human", path=None, *, labelled=False):
    """Build a fresh control index from a raw table or bundled human TRB junctions.

    Returns the index and raw-source provenance; never reads/writes a persisted index.
    Supplied tables explicitly declare junctions of the requested species and locus.
    """
    import gzip
    import hashlib
    import importlib.resources
    from pathlib import Path
    from ..db.cache import sha256

    import polars as pl
    from seqtree import Index

    species = _organism(species)
    if locus not in {"TRA", "TRB"}:
        raise ValueError("raw controls support TRA or TRB only")
    if path:
        path = Path(path)
        suffix = path.with_suffix("").suffix if path.suffix == ".gz" else path.suffix
        separator = "," if suffix == ".csv" else "\t"
        headers = pl.read_csv(path, separator=separator, n_rows=0).columns
        column = next(
            (c for c in ("junction_aa", "cdr3_aa", "cdr3", "cdr3aa") if c in headers),
            None,
        )
        if column is None:
            raise ValueError("control requires a junction_aa or legacy junction column")
        if column != "junction_aa" and (
            any(c in headers for c in ("v_call", "j_call", "sequence_id"))
        ):
            raise ValueError(
                "AIRR control requires junction_aa; bare CDR3 cannot be calibrated"
            )
        genes = {kind: next((c for c in aliases if c in headers), None) for kind, aliases in
            [('v', ('v_call', 'v_gene', 'v.segm', 'v')), ('j', ('j_call', 'j_gene', 'j.segm', 'j'))]}
        frame = pl.read_csv(
            path, separator=separator, columns=list(dict.fromkeys([column, *[c for c in genes.values() if c], *[c for c in ('species','organism','locus') if c in headers]])) if labelled else [column], infer_schema_length=0
        )
        frame = frame.with_columns(pl.col(column).str.strip_chars().str.to_uppercase())
        valid = frame.filter(pl.col(column).str.contains(r"^[ACDEFGHIKLMNPQRSTVWY]+$"))
        sequences = valid[column].unique().sort().to_list()
        provenance = {
            "kind": "explicit_raw_table",
            "path": str(path.resolve()),
            "sha256": sha256(path),
            "source_rows": frame.height,
            "productive_rows": valid.height,
            "excluded_rows": frame.height - valid.height,
            "declared_species": species,
            "declared_locus": locus,
            "sequence_column": column,
        }
    else:
        if _organism(species) != "human" or locus != "TRB":
            raise ValueError(
                "this control is not bundled; supply --control with the matching raw species/locus repertoire"
            )
        asset = importlib.resources.files("seqtree").joinpath(
            "data/control_human_trb_aa.txt.gz"
        )
        payload = asset.read_bytes()
        raw_sequences = gzip.decompress(payload).decode().splitlines()
        sequences = pl.Series(raw_sequences).unique().sort().to_list()
        provenance = {
            "kind": "bundled_raw",
            "asset": "seqtree/data/control_human_trb_aa.txt.gz",
            "sha256": hashlib.sha256(payload).hexdigest(),
            "source_rows": len(raw_sequences),
        }
    if not sequences:
        raise ValueError("control has no valid junctions for the requested locus")
    provenance["unique_junctions"] = len(sequences)
    if labelled:
        if not path:
            raise ValueError('gene-labelled calibration requires an explicit raw table')
        labelled_frame = valid.select(pl.col(column).alias('cdr3'),
            *[pl.col(c).alias(kind) if c else pl.lit(None,dtype=pl.String).alias(kind) for kind,c in genes.items()], *[c for c in ('species','organism','locus') if c in valid.columns])
        labelled_frame = control_population(labelled_frame,species=species,locus=locus)
        provenance['population_rows']=labelled_frame.height
        labelled_frame = impute_gene_calls(labelled_frame, species=species, locus=locus)
        provenance.update(observation_unit='distinct full junction within each allele-pooled V/J predicate; abundance ignored',
            gene_columns=genes, gene_imputation='label-independent observed mode; bundled model mode if no supplied calls',
            v_imputed_rows=int(labelled_frame['v_imputed'].sum()), j_imputed_rows=int(labelled_frame['j_imputed'].sum()))
        return labelled_frame, provenance
    return Index.build(sequences, "aa"), provenance


def raw_gene_background(path, locus="TRB", species="human", *, deduplicate=True, include_j=False):
    """Fresh eligible junction/V-allele keys for total-distance control counts.

    Uses supplied V annotations; never chooses one conflicting V call or weights
    rows by abundance. Distinct V assignments for the same junction are distinct
    keys when deduplicating; generated draw multiplicities can be retained.
    Present species/locus columns are filtered; absent columns require the
    caller's explicit declaration. With ``include_j``, returned annotation rows
    retain conflicting J calls for boundary consensus; population counts still
    use junction/V keys, and callers collapse annotation rows before counting.
    """
    from pathlib import Path
    import polars as pl
    from ..db.cache import sha256
    from ..match.geometry import load_v_loops, resolve_v_alleles

    organism = _organism(species)
    if locus not in {"TRA", "TRB"}:
        raise ValueError("gene controls require TRA or TRB")
    path = Path(path)
    suffix = path.with_suffix("").suffix if path.suffix == ".gz" else path.suffix
    separator = "," if suffix == ".csv" else "\t"
    headers = pl.read_csv(path, separator=separator, n_rows=0).columns
    sequence = next((c for c in ("junction_aa", "cdr3_aa", "cdr3", "cdr3aa") if c in headers), None)
    v = next((c for c in ("v_call", "v_gene", "v.segm", "v") if c in headers), None)
    if sequence is None or v is None:
        raise ValueError("total-distance controls require full junction and V columns")
    if sequence != "junction_aa" and any(c in headers for c in ("v_call", "j_call", "sequence_id")):
        raise ValueError("AIRR control requires junction_aa; bare CDR3 cannot be calibrated")
    j = next((c for c in ("j_call", "j_gene", "j.segm", "j") if c in headers), None) if include_j else None
    species_column = next((c for c in ("species", "organism") if c in headers), None)
    locus_column = next((c for c in ("locus", "chain") if c in headers), None)
    raw = pl.read_csv(path, separator=separator,
                      columns=[sequence, v, *[c for c in (species_column, locus_column, j) if c]],
                      infer_schema_length=0)
    selected = raw
    if species_column:
        aliases = {}
        for value in raw[species_column].drop_nulls().unique():
            try:
                aliases[value] = _organism(value)
            except ValueError:
                aliases[value] = None
        selected = selected.filter(pl.col(species_column).replace_strict(
            aliases, default=None, return_dtype=pl.String).eq(organism).fill_null(False))
    if locus_column:
        selected = selected.filter(pl.col(locus_column).str.strip_chars().str.to_uppercase()
            .replace({'ALPHA':'TRA','BETA':'TRB'}).eq(locus).fill_null(False))
    frame = selected.select(pl.col(sequence).str.strip_chars().str.to_uppercase().alias("cdr3"),
                       pl.col(v).alias("v"),
                       *([pl.col(j).alias("j") if j else pl.lit(None,dtype=pl.String).alias("j")] if include_j else []))
    frame = frame.with_columns(resolve_v_alleles(frame["v"], load_v_loops(organism), locus=locus).alias("_allele"))
    valid = frame.filter(pl.col("cdr3").str.contains(r"^[ACDEFGHIKLMNPQRSTVWY]{8,}$").fill_null(False))
    eligible = valid.filter(pl.col("_allele").is_not_null())
    keys = eligible.select("cdr3", "_allele", *(["j"] if include_j else []))
    unique_keys = keys.select("cdr3", "_allele").unique()
    keys = (keys.unique() if deduplicate else keys).sort("cdr3", "_allele")
    if not keys.height:
        raise ValueError("control has no eligible junction/V-allele keys")
    return keys, {
        "kind": "explicit_raw_gene_table", "path": str(path.resolve()), "sha256": sha256(path),
        "declared_species": organism, "declared_locus": locus,
        "sequence_column": sequence, "v_column": v, "source_rows": raw.height,
        "species_column": species_column, "locus_column": locus_column,
        "excluded_metadata_rows": raw.height-selected.height,
        "excluded_sequence_rows": frame.height-valid.height,
        "excluded_gene_rows": valid.height-eligible.height,
        "eligible_rows": eligible.height, "unique_keys": unique_keys.height,
        "population_size": unique_keys.height if deduplicate else keys.height, "deduplicated": deduplicate,
        "observation_unit": ("distinct full junction + resolved V allele; abundance ignored" if deduplicate else
                             "eligible generated receptor draw; repeated junction/V keys retained"),
        **({"j_column":j,"j_annotation_policy":"retained for boundary consensus; conflicts never select a first call"} if include_j else {}),
        "v_annotation_policy": "source-supplied calls; ambiguous/unmodeled excluded, no first-call choice",
    }


def impute_gene_calls(frame, *, species, locus, kinds=('v', 'j'), reference=None):
    """Retain source rows/calls and flag label-independent modal gene imputation."""
    import polars as pl
    from .mock_reference import _modal_gene

    for kind in kinds:
        if kind not in frame.columns:
            frame=frame.with_columns(pl.lit(None,dtype=pl.String).alias(kind))
        genes=pl.col(kind).cast(pl.String).str.to_uppercase().str.split(',').list.eval(pl.element().str.split('*').list.first().str.strip_chars()).list.unique()
        value=genes.list.first()
        valid=((genes.list.len()==1)&value.str.contains('^'+locus+kind.upper()+r'[A-Z0-9/.-]+$')).fill_null(False)
        prior=pl.col(kind+'_imputed') if kind+'_imputed' in frame.columns else pl.lit(False)
        original=pl.col('original_'+kind) if 'original_'+kind in frame.columns else pl.col(kind)
        frame=frame.with_columns(original.alias('original_'+kind), (prior|~valid).alias(kind+'_imputed'))
        if frame[kind+'_imputed'].any():
            population=frame if reference is None else reference
            pool=population.select(genes.alias('_genes')).filter(pl.col('_genes').list.len()==1).select(pl.col('_genes').list.first().alias(kind)).filter(pl.col(kind).str.contains('^'+locus+kind.upper()+r'[A-Z0-9/.-]+$').fill_null(False))
            replacement=_modal_gene(pool,kind,species=species,locus=locus,kind=kind)
            frame=frame.with_columns(pl.when(valid).then(value).otherwise(pl.lit(replacement)).alias(kind))
        else:
            frame=frame.with_columns(value.alias(kind))
    return frame


def conditioned_counts(queries, controls, params, *, match_v=False, match_j=False, threads=1, exclude_exact=False):
    """One native histogram batch per fixed gene-conditioned control population.

    Count distinct full junctions within each allele-pooled V/J predicate. Model
    sizes remain fixed under exact hit exclusion; read abundance is ignored.
    Returns query/distance CDF, full-scope counts and per-query model sizes.
    """
    import numpy as np
    import polars as pl
    from seqtree import Index

    if not params.max_total_edits and any((params.max_subs,params.max_ins,params.max_dels)):
        raise ValueError('gene-conditioned controls require an explicit total-edit cap')
    kinds=[k for k,flag in [('v',match_v),('j',match_j)] if flag]
    if not kinds or not {'cdr3',*kinds} <= set(controls.columns):
        raise ValueError('gene-conditioned controls require junction and requested gene columns')
    parts=[];sizes=[]
    controls=controls.select('cdr3',*kinds).unique()
    for key,q in queries.partition_by(kinds,as_dict=True,maintain_order=True).items():
        population=controls.filter(pl.all_horizontal([pl.col(k).eq(value) for k,value in zip(kinds,key)]))['cdr3'].unique().to_list()
        sequences=q['cdr3'].unique(maintain_order=True).to_list()
        if population and params.max_total_edits:
            histogram=np.asarray(Index.build(population,'aa').edit_histogram_batch(sequences,params,threads,False),dtype=np.uint64)
        else:
            histogram=np.zeros((len(sequences),params.max_total_edits+1),dtype=np.uint64)
            present=set(population)
            histogram[:,0]=[sequence in present for sequence in sequences]
        if exclude_exact:histogram[:,0]=0
        cumulative=histogram.cumsum(axis=1,dtype=np.uint64)
        parts.append(pl.DataFrame({'cdr3':np.repeat(sequences,cumulative.shape[1]),'edits':np.tile(np.arange(cumulative.shape[1],dtype=np.uint32),len(sequences)),'_nc':cumulative.ravel()}).join(q.select('query_id','cdr3'),on='cdr3',validate='m:m').drop('cdr3'))
        sizes.append(q.select('query_id').with_columns(pl.lit(len(population),dtype=pl.UInt64).alias('_control_size')))
    cdf=pl.concat(parts)
    totals=cdf.filter(pl.col('edits')==params.max_total_edits).select('query_id',pl.col('_nc').alias('n_control'))
    return cdf,totals,pl.concat(sizes)


def control_population(frame, *, species, locus):
    """Select supplied population metadata before imputing genuinely missing calls."""
    import polars as pl
    organism=_organism(species)
    for column in ('species','organism'):
        if column in frame.columns:
            aliases={value:_organism(value) for value in frame[column].drop_nulls().unique() if str(value).strip()}
            frame=frame.filter(pl.col(column).replace_strict(aliases,default=None,return_dtype=pl.String).fill_null(organism)==organism)
    gene_cols=[c for c in ('v','j') if c in frame.columns]
    alpha=pl.any_horizontal([pl.col(c).cast(pl.String).str.to_uppercase().str.contains(r'TRA[VJ]').fill_null(False) for c in gene_cols]) if gene_cols else pl.lit(False)
    beta=pl.any_horizontal([pl.col(c).cast(pl.String).str.to_uppercase().str.contains(r'TRB[VJ]').fill_null(False) for c in gene_cols]) if gene_cols else pl.lit(False)
    # Contradictory segment calls cannot infer a locus. Keep the explicitly
    # supplied population; impute_gene_calls retains/flags the original calls.
    inferred=pl.when(alpha&beta).then(None).when(alpha).then(pl.lit('TRA')).when(beta).then(pl.lit('TRB')).otherwise(None)
    declared=pl.col('locus').cast(pl.String).str.strip_chars().str.to_uppercase().replace({'ALPHA':'TRA','BETA':'TRB'}) if 'locus' in frame.columns else pl.lit(None,dtype=pl.String)
    declared=pl.when(declared=='').then(None).otherwise(declared)
    if frame.filter(inferred.is_not_null()&declared.is_not_null()&(inferred!=declared)).height:
        raise ValueError('control V/J locus conflicts with explicit locus metadata')
    return frame.filter(pl.coalesce(declared,inferred,pl.lit(locus))==locus)
