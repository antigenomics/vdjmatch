"""Reproduce accepted legacy matrix and epitope-withheld retrieval tables."""
from pathlib import Path


def register(subparsers):
    p = subparsers.add_parser('scoring-audit', help='derive origin-aware matrices or reproduce legacy matrix/LOO tables')
    p.add_argument('--task', choices=['matrix', 'loo', 'origin-matrix'], required=True)
    p.add_argument('--pmhc', required=True, help='raw VDJdb export; legacy cdr3 denotes anchored junction')
    p.add_argument('--species', choices=['HomoSapiens', 'MusMusculus'], default='HomoSapiens')
    p.add_argument('--locus', choices=['TRA', 'TRB'], required=True)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--seed', type=int, default=0, help='legacy composition-cap seed (LOO only)')
    p.add_argument('--min-epi', type=int, default=30, help='legacy unique-junction eligibility; not the full-VDJdb source-record gate')
    p.add_argument('--min-records', type=int, default=30, help='origin-matrix: valid source-record eligibility before deduplication')
    p.add_argument('--origin-neighbours', action='store_true', help='origin-matrix: derive exact flank-conditioned N-region neighbour models and source-held-out predictions')
    p.add_argument('--origin-rates', action='store_true', help='origin-matrix: all-site N-region escape exposures, retaining zero-neighbour sites; requires vdjtools>=5.8.0')
    p.add_argument('--origin-mixed', action='store_true', help='origin-matrix: full local G/P/N/template posterior and exchangeable block matrices; requires vdjtools>=5.9.0')
    p.add_argument('--vlen', type=int, default=4, help='matrix task V prefix length')
    p.add_argument('--jlen', type=int, default=6, help='matrix task J suffix length')
    p.add_argument('--legacy', help='explicit legacy similarity TSV; default packaged vdjam.txt')
    p.add_argument('--retention', help='explicit V/J retention model TSV; default packaged human model')
    p.add_argument('--significance', help='explicit positional significance model TSV; default packaged model')
    p.add_argument('--structural-inc', help='explicit seqtree structural.inc; omit to skip structural row, never alias BLOSUM')
    p.add_argument('--top', type=int, default=8, help='LOO held-out epitopes by descending unique-junction count')
    p.add_argument('--max-queries', type=int, default=600, help='LOO lexicographically first unique queries per epitope')
    p.add_argument('--subs', type=int, default=2, help='LOO substitution-only candidate radius')
    p.add_argument('--output-prefix', required=True)
    p.set_defaults(func=main)
    return p


def main(a):
    if a.task == 'origin-matrix':
        return origin_matrix(a)
    if a.origin_neighbours or a.origin_rates or a.origin_mixed:
        raise ValueError('--origin-neighbours/--origin-rates/--origin-mixed require --task origin-matrix')
    import json
    import math
    import time
    from importlib import metadata, resources
    import polars as pl
    from .. import db
    from ..match import regions, scoring_audit as legacy
    from ..db.cache import sha256
    from ..runner.annotation import _validate_output_paths

    for name in ('threads', 'min_epi', 'top', 'max_queries'):
        if getattr(a, name) < 1:
            raise ValueError(f'{name.replace("_", "-")} must be positive')
    if min(a.vlen, a.jlen, a.subs) < 0:
        raise ValueError('region lengths and substitution radius must be nonnegative')
    if a.task == 'loo' and a.species == 'MusMusculus' and not a.retention:
        raise ValueError('mouse LOO requires an explicit species-matched --retention model')
    start = time.perf_counter()
    root = resources.files('vdjmatch.resources')
    paths = {'legacy': Path(a.legacy or root / 'vdjam.txt')} if a.task == 'matrix' else {}
    if a.task == 'loo':
        paths.update(retention=Path(a.retention or root / 'trimming' / 'human_vj_retention.tsv'),
                     significance=Path(a.significance or root / 'trimming' / 'position_significance.tsv'))
        if a.structural_inc:
            paths['structural'] = Path(a.structural_inc)
    prefix = _validate_output_paths(a.output_prefix, [a.pmhc,*paths.values(),*[path for path in
        (a.legacy,a.retention,a.significance,a.structural_inc) if path]])
    models = {name: {'path': str(path.resolve()), 'sha256': sha256(path)} for name, path in paths.items()}
    raw = db.load(a.pmhc, species=a.species, gene=a.locus)
    # The legacy estimators use standard-residue sanitation, not modern curation.
    valid = legacy.valid_cdr3(raw)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    if a.task == 'matrix':
        matrices, correlations = legacy.matrix_audit(raw, min_epi=a.min_epi, vlen=a.vlen,
            jlen=a.jlen, legacy=legacy.legacy_vdjam_offdiag(paths['legacy']), threads=a.threads)
        tables = {'matrix': matrices, 'correlations': correlations}
    else:
        structural = None
        if a.structural_inc:
            structural = legacy.structural_dissim(a.structural_inc)
            if structural is None or len(structural) != len(legacy.AA)*(len(legacy.AA)-1):
                raise ValueError('structural-inc does not supply the complete standard-AA structural matrix')
        retention = regions.load_retention(paths['retention'])
        significance = regions.load_significance(paths['significance'])
        if not retention or not all(significance):
            raise ValueError('retention and significance models must be nonempty')
        if any(not math.isfinite(p) or not 0 <= p <= 1
               for profile in [*retention.values(), *significance] for p in profile):
            raise ValueError('retention/significance probabilities must be finite and in [0,1]')
        tasks, summary, pairs, queries, refs = legacy.loo_audit(raw, locus=a.locus, retention=retention,
            significance=significance, structural=structural, min_epi=a.min_epi, top=a.top,
            max_queries=a.max_queries, subs=a.subs, seed=a.seed, threads=a.threads)
        tables = {'tasks': tasks, 'summary': summary, 'pairs': pairs, 'queries': queries,
                  'references': refs.with_row_index('reference_key_id')}
    for name, table in tables.items():
        table.with_columns(pl.lit(a.species).alias('species'), pl.lit(a.locus).alias('locus')).write_csv(
            str(prefix)+'.'+name+'.tsv', separator='\t')
    manifest = {
        'reference': db.provenance(a.pmhc), 'models': models,
        'options': vars(a) | {'func': 'scoring-audit'},
        'versions': {name: metadata.version(name) for name in ('vdjmatch', 'seqtree', 'polars', 'numpy')},
        'source_chain_records': raw.height, 'valid_aa_chain_records': valid.height,
        'invalid_aa_chain_records': raw.height-valid.height,
        'contract': 'legacy accepted evidence; not current annotation or full-VDJdb benchmark eligibility',
        'sequence_convention': 'producer-normalized anchored junction; no fabricated anchors',
        'matrix_estimator': 'same-epitope unordered Hamming-1 pairs; position background; alpha=1 log2 odds',
        'loo_contract': 'spike-study records removed; unique(cdr3,v,j,epitope,mhc_class); cap=3000; seeded canonical input; first epitope per junction; exact sequence pairs excluded; NDN training positions [4,L-6)',
        'metric': 'legacy tied-score-group trapezoidal balanced PR-AUC at reference prevalence 0.5; not average precision',
        'structural_status': 'explicit input' if a.structural_inc else 'omitted; no BLOSUM fallback',
        'output_rows': {name: table.height for name, table in tables.items()},
        'seconds': time.perf_counter()-start,
    }
    Path(str(prefix)+'.provenance.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return 0


def origin_matrix(a):
    """Fresh exact-origin block derivation; no legacy estimator/model inputs."""
    import json
    import time
    import sys
    from importlib import metadata, resources
    import polars as pl
    from arda.cdr3fix import load_anchors, resolve_allele
    from vdjtools.model import load_bundled, n_region_posterior_batch
    from vdjtools.model.native import pack, _gene_idx
    import vdjtools._core as native
    from .. import db
    from ..db.cache import sha256
    from ..match.scoring_audit import origin_block_matrices, origin_neighbour_tables, origin_rate_tables, origin_substitution_pairs, mixed_origin_block_matrices
    from ..runner.annotation import _validate_output_paths
    from .epitope import pool_epitopes, write_tables

    if a.threads < 1 or a.min_records < 1:
        raise ValueError("threads and min-records must be positive")
    if any(
        getattr(a, k, None)
        for k in ["legacy", "retention", "significance", "structural_inc"]
    ):
        raise ValueError("origin-matrix does not accept legacy scoring models")
    start = time.perf_counter()
    organism = {"HomoSapiens": "human", "MusMusculus": "mouse"}[a.species]
    model_root = (
        resources.files("vdjtools.model")
        / "_bundled"
        / "isalgo"
        / (a.locus if organism == "human" else "mouse_" + a.locus)
    )
    model_files = sorted(Path(p) for p in model_root.iterdir() if p.is_file())
    prefix = _validate_output_paths(a.output_prefix, [a.pmhc, *model_files])
    raw = db.load(a.pmhc, species=a.species)
    pooled, cohort = pool_epitopes(raw, min_records=a.min_records)
    source = pooled.filter(
        (pl.col("gene") == a.locus)
        & pl.col("reference_valid").fill_null(False)
        & pl.col("cdr3").str.contains("^[ACDEFGHIKLMNPQRSTVWY]+$")
    ).select(
        "epitope", "cdr3", pl.col("v").alias("source_v"), pl.col("j").alias("source_j"),
        (pl.col("reference_id") if "reference_id" in pooled.columns else pl.lit(None, dtype=pl.String)).alias("reference_id")
    )
    all_labels = raw.filter(
        (pl.col("gene") == a.locus)
        & pl.col("reference_valid").fill_null(False)
        & pl.col("cdr3").str.contains("^[ACDEFGHIKLMNPQRSTVWY]+$")
    ).select(
        "epitope", "cdr3", pl.col("v").alias("source_v"), pl.col("j").alias("source_j")
    )
    model = load_bundled(a.locus, source="isalgo", organism=organism)
    _, vi, ji = pack(model)
    anchors = None
    call_rows = []
    for column, kind, index in [("source_v", "V", vi), ("source_j", "J", ji)]:
        inverse = {i: name for name, i in index.items()}
        mapping = {}
        for call in all_labels[column].unique().sort(nulls_last=True).to_list():
            resolved = None
            status = "unresolved"
            ambiguous = bool(call and (any(c in call for c in ",;| ") or "/TR" in call))
            if ambiguous:
                status = "ambiguous_call"
            if call and not ambiguous:
                try:
                    resolved = inverse[_gene_idx(index, call, kind)]
                    status = "model_resolved"
                except KeyError:
                    if anchors is None:
                        anchors = load_anchors(organism)
                    candidate = resolve_allele(call, kind, anchors, organism=organism)
                    if candidate:
                        try:
                            resolved = inverse[_gene_idx(index, candidate, kind)]
                            status = "arda_resolved_absent_call"
                        except KeyError:
                            pass
            mapping[call] = resolved
            call_rows.append(
                dict(segment=kind, call=call, resolved=resolved, status=status)
            )
        source = source.with_columns(
            pl.col(column)
            .replace_strict(mapping, return_dtype=pl.String)
            .alias(kind.lower())
        )
        all_labels = all_labels.with_columns(
            pl.col(column)
            .replace_strict(mapping, return_dtype=pl.String)
            .alias(kind.lower())
        )
    ready = source.filter(pl.col("v").is_not_null() & pl.col("j").is_not_null())
    queries = (
        ready.select("cdr3", "v", "j")
        .unique()
        .sort("cdr3", "v", "j")
        .with_row_index("query_id")
        .with_columns(pl.col("query_id").cast(pl.Int64))
    )
    origin_function=n_region_posterior_batch
    classes=None
    emission_rows=[]
    if a.origin_mixed:
        from vdjtools.model import codon_origin_posterior_batch
        origin_function=codon_origin_posterior_batch
        for name,table in model.tables.items():
            if name.endswith(('_dinucl','_first_nt')):
                minimum=float(table['p'].min())
                if not minimum>0:
                    raise ValueError('mixed block eligibility needs strictly positive N emissions; use exact alternative masses for a zero-emission model')
                emission_rows.append(dict(table=name,probabilities=table.height,minimum_probability=minimum))
        if not emission_rows:
            raise ValueError('mixed block eligibility requires declared N emission tables')
        stage_start=time.perf_counter()
        print(json.dumps(dict(stage='local_origin_posterior',queries=queries.height,native_threads=a.threads)),file=sys.stderr,flush=True)
        classes=origin_function(model,queries['cdr3'].to_list(),v=queries['v'].to_list(),
            j=queries['j'].to_list(),threads=a.threads).rename({'row':'query_id'}).with_columns(
            (pl.col('codon_amino_acids').str.len_chars()>1).fill_null(False).alias('exchangeable'))
        print(json.dumps(dict(stage='local_origin_posterior_complete',classes=classes.height,
            seconds=time.perf_counter()-stage_start)),file=sys.stderr,flush=True)
        origins=classes.group_by('query_id','position').agg(
            pl.when(pl.col('probability').is_not_null().any())
            .then(pl.col('probability').filter(pl.col('origin')=='NNN').sum())
            .otherwise(None).alias('n_probability')).sort('query_id','position')
    else:
        origins = n_region_posterior_batch(
            model,
            queries["cdr3"].to_list(),
            v=queries["v"].to_list(),
            j=queries["j"].to_list(),
            threads=a.threads,
        ).rename({"row": "query_id"})
    membership = (
        ready.group_by("epitope", "cdr3", "v", "j")
        .len(name="source_records")
        .join(queries, on=["cdr3", "v", "j"], validate="m:1")
        .sort("epitope", "query_id")
    )
    known_labels = (
        all_labels.select("cdr3", "v", "j", "epitope")
        .unique()
        .join(queries, on=["cdr3", "v", "j"], how="inner", validate="m:1")
        .select("query_id", "epitope")
    )
    matrix, blocks, support, epitope_pairs, contrast = origin_block_matrices(
        membership, origins, known_labels
    )
    resolved = membership.join(support, on="query_id", validate="m:1").with_columns(
        pl.when(pl.col("origin_supported"))
        .then(pl.lit("origin_supported"))
        .otherwise(pl.lit("zero_sequence_likelihood"))
        .alias("status"),
        (pl.col("cdr3").str.len_chars() < 8).alias("legacy_length_lt8"),
    )
    accounting = resolved.group_by("epitope", "status", "legacy_length_lt8").agg(
        pl.col("source_records").sum(), pl.len().alias("clonotypes")
    )
    missing = (
        source.filter(pl.col("v").is_null() | pl.col("j").is_null())
        .with_columns((pl.col("cdr3").str.len_chars() < 8).alias("legacy_length_lt8"))
        .group_by("epitope", "legacy_length_lt8")
        .agg(
            pl.len().alias("source_records"),
            pl.struct("cdr3", "source_v", "source_j").n_unique().alias("clonotypes"),
        )
        .with_columns(pl.lit("unresolved_call").alias("status"))
    )
    length_accounting = pl.concat(
        [accounting, missing.select(accounting.columns)], how="vertical_relaxed"
    ).sort("epitope", "status", "legacy_length_lt8")
    accounting = (
        length_accounting.group_by("epitope", "status")
        .agg(pl.col("source_records").sum(), pl.col("clonotypes").sum())
        .sort("epitope", "status")
    )
    if int(accounting["source_records"].sum()) != source.height:
        raise RuntimeError("origin accounting lost source records")
    sites = (
        origins.join(queries, on="query_id", validate="m:1")
        .join(membership.select("query_id", "epitope"), on="query_id", validate="m:m")
        .filter(pl.col("n_probability").is_not_null())
        .with_columns(
            pl.col("cdr3").str.len_chars().alias("length"),
            pl.col("cdr3").str.slice(pl.col("position"), 1).alias("aa"),
        )
    )
    exposures = sites.group_by("epitope", "length", "position", "aa").agg(
        pl.col("n_probability").sum().alias("n_site_weight"),
        pl.len().alias("clonotype_sites"),
    )
    coordinate = ["epitope", "length", "position"]
    exposures = (
        exposures.with_columns(
            pl.col("n_site_weight").sum().over(coordinate).alias("n_coordinate_weight")
        )
        .with_columns(
            pl.when(pl.col("n_coordinate_weight") > 0)
            .then(pl.col("n_site_weight") / pl.col("n_coordinate_weight"))
            .otherwise(None)
            .alias("aa_probability_given_N")
        )
        .sort(*coordinate, "aa")
    )
    tables = {
        "matrix": matrix,
        "epitope-pairs": epitope_pairs,
        "exposures": exposures,
        "blocks": blocks,
        "support": accounting,
        "length-band-support": length_accounting,
        "calls": pl.DataFrame(
            call_rows,
            schema={
                "segment": pl.String,
                "call": pl.String,
                "resolved": pl.String,
                "status": pl.String,
            },
        ).sort("segment", "call", nulls_last=True),
    }
    tables.update(contrast)
    if classes is not None:
        stage_start=time.perf_counter()
        print(json.dumps(dict(stage='mixed_origin_block_counts',classes=classes.height)),file=sys.stderr,flush=True)
        tables.update(mixed_origin_block_matrices(membership,classes))
        print(json.dumps(dict(stage='mixed_origin_block_counts_complete',seconds=time.perf_counter()-stage_start)),file=sys.stderr,flush=True)
        tables['mixed-emissions']=pl.DataFrame(emission_rows)
        local=classes.join(queries.select('query_id','cdr3'),on='query_id',validate='m:1').join(
            membership.select('query_id','epitope'),on='query_id',validate='m:m').with_columns(
            pl.col('cdr3').str.len_chars().alias('length'),
            pl.col('cdr3').str.slice(pl.col('position'),1).alias('aa'))
        tables['mixed-exposures']=local.group_by('epitope','length','position','aa','origin','exchangeable').agg(
            pl.when(pl.col('probability').is_not_null().any()).then(pl.col('probability').sum())
            .otherwise(None).alias('posterior_site_weight'),
            pl.col('query_id').n_unique().alias('clonotype_sites'),
            pl.col('probability').is_null().sum().alias('unknown_origin_sites')).sort(
                'epitope','length','position','aa','origin','exchangeable')

    pair_data = origin_substitution_pairs(queries, origins, threads=a.threads) if a.origin_neighbours or a.origin_rates else None
    if a.origin_rates:
        tables.update(origin_rate_tables(model, queries, membership, origins, threads=a.threads, pair_data=pair_data))
    if a.origin_neighbours:
        source_membership = ready.select("epitope", "cdr3", "v", "j", "reference_id").join(
            queries, on=["cdr3", "v", "j"], validate="m:1").select("query_id", "epitope", "reference_id").unique()
        tables.update(origin_neighbour_tables(model, queries, membership, origins, known_labels,
            source_membership=source_membership, threads=a.threads, pair_data=pair_data))
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for name, table in tables.items():
        table.with_columns(
            pl.lit(a.species).alias("species"), pl.lit(a.locus).alias("locus")
        ).write_csv(str(prefix) + "." + name + ".tsv", separator="\t")
    write_tables(prefix, cohort)
    manifest = dict(
        reference=db.provenance(a.pmhc),
        options=vars(a) | {"func": "scoring-audit"},
        models={p.name: sha256(p) for p in model_files},
        versions={
            name: metadata.version(name)
            for name in [
                "vdjmatch",
                "vdjtools",
                "arda-mapper",
                "seqtree",
                "polars",
                "numpy",
            ]
        },
        source_records=source.height,
        origin_queries=queries.height,
        seconds=time.perf_counter() - start,
        source_sha256={
            name: sha256(Path(path))
            for name, path in {
                "cli": __file__,
                "counter": origin_block_matrices.__code__.co_filename,
                "cohort": pool_epitopes.__code__.co_filename,
                "origin_wrapper": origin_function.__code__.co_filename,
                "native": native.__file__,
                "call_resolution": resolve_allele.__code__.co_filename,
            }.items()
        },
        sequence_convention="anchor-inclusive junction; zero-based AA codons",
        estimator="pure-N posterior products; same-epitope/V/J/length blocks; unique junction/V/J unit clusters",
        recognition_background="same V/J/length/position/N opportunity; complete known curated label sets disjoint; unmatched positive mass retained separately",
        recognition_scores=[
            "unordered count density ratio",
            "mean directional log odds; single-stratum independent composition cancels",
        ],
        neighbour_estimator=(dict(
            opportunity="20 exact AA alternatives conditional on unchanged flanks, V/J and pure N; products consumed in memory",
            composition="shared epitope-residue log preferences; zero-row-sum pair-effect gauge; composition-only baseline refitted",
            prior="one total neutral pseudo-event per pair/preference block; no reference-matrix prior",
            validation="largest recorded source per epitope held out; every held-out junction excluded from training",
            completion="reversible compatibility walk with matched N-site stationary exposure at fixed0.70bits information",
            pair_units="distinct junction/V/J pair incidences per eligible epitope; N-posterior products are fractional weights",
        ) if a.origin_neighbours else None),
        rate_estimator=(dict(
            exposure="all source sites; source N posterior times exact19-alternative N mass times distinct model-supported partner count",
            units="directed endpoint incidences and expected opportunities; no-neighbour and zero/unknown N sites retained",
            population="epitope/V/J/junction length; constant model mass normalizer belongs to nuisance intercept",
            interpretation="working convergent-neighbour intensity, not an ancestral mutation probability",
        ) if a.origin_rates else None),
        mixed_estimator=(dict(
            posterior="complete local G/P/N/template partition",
            eligibility="genetic-code alternatives with immutable letters fixed; all declared N emissions strictly positive",
            channels="same exchangeable; same zero alternative; different local class",
            counts="ordered aligned pairs between unique junction/V/J unit clusters; mandatory identities excluded from the exchangeable core",
            NNN="marginal of the same complete posterior; no second origin batch",
            limitation="local class coarsens full histories; model posterior is not observed nucleotide truth",
        ) if a.origin_mixed else None),
        known_label_source_records=all_labels.height,
        known_label_query_memberships=known_labels.height,
        length_policy="all producer-declared junction lengths retained; no anchor motif screen or fabrication; legacy_length_lt8 is an audit band only",
        weighting=["pair-mass pooled", "equal contributing epitope"],
        score="unrounded log2(q_ab/(p_a*p_b)); distance=max(0,S_aa+S_bb-2*S_ab), bits",
        missing="zero-likelihood and unresolved calls retained in support; unsupported cells null; no pseudocount",
        output_rows={k: t.height for k, t in tables.items()},
    )
    Path(str(prefix) + ".provenance.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    return 0
