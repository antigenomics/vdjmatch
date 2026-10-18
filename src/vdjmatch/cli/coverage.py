"""Occurrence-weighted clonotype/Leiden coverage, with matched control populations."""
from pathlib import Path


def register(subparsers):
    p = subparsers.add_parser('reference-coverage', help=__doc__)
    p.add_argument('--vdjdb', required=True)
    p.add_argument('--species', default='human')
    p.add_argument('--locus', choices=['TRA', 'TRB'], required=True)
    p.add_argument('--epitope', action='append')
    p.add_argument('--targets-from-sample')
    p.add_argument('--reference-population', choices=['curated-benchmark', 'all-valid'], default='curated-benchmark')
    add_graph_options(p)
    p.add_argument('--control', help='raw joint junction/V/J recombination-event control')
    p.add_argument('--control-model', help='prepared immutable control directory from prepare-control')
    p.add_argument('--mock-replicates', type=int, default=0)
    p.add_argument('--pgen-features', help='validated PREFIX.features.parquet Hamming-1 bank for control matching')
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--output-prefix', required=True)
    from .epitope import add_options
    add_options(p)
    p.set_defaults(func=main)


def add_graph_options(p, *, combined=False):
    p.add_argument('--geometry', choices=['levenshtein', 'sw', *(['both'] if combined else [])],
                   default='both' if combined else None, required=not combined)
    p.add_argument('--scope', default='1,1,1,1')
    p.add_argument('--sw-radius', type=int, help='explicit inclusive radius in summed raw BLOSUM62 self-score-distance units')
    p.add_argument('--sw-radius-knee', action='store_true', help='derive a global within-pMHC nearest-neighbour SW Kneedle radius')
    p.add_argument('--cdr-loops', choices=['junction', 'cdr12', 'cdr12plus'], default='cdr12plus')
    p.add_argument('--gap-open', type=int, default=11)
    p.add_argument('--gap-extend', type=int, default=1)
    p.add_argument('--leiden-resolution', type=float, default=1.)
    p.add_argument('--nboot', type=int, default=200, help='vdjtools augmented-assemblage bootstrap replicates')
    p.add_argument('--pair-budget', type=int, default=2_000_000)
    p.add_argument('--membership-budget-bytes', type=int, default=2*1024**3,
        help='estimated raw reference/mock membership byte ceiling (default: 2 GiB; maximum: 8 GiB)')


def main(a, *, inputs=None):
    from tempfile import TemporaryDirectory
    with TemporaryDirectory(prefix='vdjmatch-edges-') as scratch:
        return _main(a, inputs=inputs, scratch=Path(scratch))


def _main(a, *, inputs, scratch):
    import hashlib
    import json
    import time
    from importlib.metadata import version
    import polars as pl
    import numpy as np
    from .. import db
    from ..aggregate.reference_coverage import reference_graph, leiden_taxa, abundance_coverage, sw_knee
    from ..db.cache import sha256
    from ..evalue.control import _organism
    from ..evalue.mock_reference import control_events, compute_pgen, iter_mock_references, EVENT_SAMPLER, PGEN_BINNING, _pgen_model_identity, _prepare_mock_control, load_pgen_features
    from ..runner.annotation import _validate_output_paths
    from .structure import reference_groups

    start = time.monotonic()
    validate_options(a)
    membership_budget = getattr(a, 'membership_budget_bytes', 2*1024**3)
    protected = [a.vdjdb] + ([a.control] if a.control else []) + ([a.targets_from_sample] if a.targets_from_sample else [])
    if a.pgen_features:
        protected += [a.pgen_features, a.pgen_features.removesuffix('.features.parquet')+'.manifest.json']
    prefix = _validate_output_paths(a.output_prefix, protected)
    organism = _organism(a.species)
    species = {'human': 'HomoSapiens', 'mouse': 'MusMusculus'}[organism]
    from .epitope import load_reference
    if inputs is None:
        raw, epitope_info = load_reference(a, species)
    else:
        raw, epitope_info = inputs['raw'], inputs['epitope_info']
    valid = pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)
    if a.reference_population == 'curated-benchmark':
        valid &= pl.col('reference_valid')
    groups = reference_groups(raw, a)
    populations, records = [], []
    def stable_seed(*identity):
        return int.from_bytes(hashlib.sha256(json.dumps([a.seed, species, a.locus, *identity],
            sort_keys=True).encode()).digest()[:4], 'little') % 2**31
    for i, (task, source) in enumerate(groups):
        selected = source.filter(valid).select('cdr3', 'v', 'j').with_columns(pl.lit(i).alias('group_id'))
        records.append(selected)
        populations.append(dict(group_id=i, species=species, locus=a.locus, **task,
            population='reference', replicate=-1, source_records=source.height,
            retained_records=selected.height, excluded_records=source.height-selected.height))
    if not records or not sum(part.height for part in records):
        raise ValueError('no valid reference occurrences in the requested population')
    reference = pl.concat(records)
    estimated_membership_rows = reference.height * (1+a.mock_replicates)
    estimated_membership_bytes = estimated_membership_rows * 1024
    if estimated_membership_bytes > membership_budget:
        raise ValueError(f'raw reference/mock membership estimate {estimated_membership_bytes} bytes exceeds explicit budget {membership_budget} bytes')
    loops = () if a.cdr_loops == 'junction' else ('cdr1', 'cdr2') if a.cdr_loops == 'cdr12' else ('cdr1', 'cdr2', 'cdr25')
    knee_curve = None
    if a.geometry == 'sw' and a.sw_radius_knee:
        a.sw_radius, knee_curve = sw_knee(reference, species=organism, locus=a.locus,
            cdr_loops=loops, gap_open=a.gap_open, gap_extend=a.gap_extend,
            threads=a.threads, pair_budget=a.pair_budget)
    control_info, feature_identity = None, None
    if a.control:
        from .combined import _raw_control
        if inputs is None:
            if getattr(a, 'control_model', None):
                from .control_model import load_control_model
                events, _, _, control_info = load_control_model(a.control_model, a.control,
                    species=organism, locus=a.locus, modes=())
            else:
                control_raw, control_info = _raw_control(a.control, retain_missing_genes=True)
                events = control_events(control_raw, species=organism, locus=a.locus)
                del control_raw
            if a.pgen_features:
                pgen, feature_identity = load_pgen_features(a.pgen_features,
                    reference.rename({'cdr3': 'junction_aa'}), events, species=organism, locus=a.locus)
            else:
                pgen = compute_pgen(reference.rename({'cdr3': 'junction_aa'}), events,
                    species=organism, locus=a.locus, threads=a.threads)
            prepared_control = _prepare_mock_control(reference.rename({'cdr3': 'junction_aa'}), events,
                pgen=pgen, species=organism, locus=a.locus, threads=a.threads)
        else:
            prepared_control = inputs['prepared_control']
            feature_identity = inputs['pgen_features']
            events, pgen, control_info = inputs['events'], inputs['pgen'], dict(inputs['control_info'])
        control_info.update(sha256=sha256(Path(a.control)), events=int(events['event_weight'].sum()),
            event_keys=events.height, sampler=EVENT_SAMPLER, replicates=a.mock_replicates, pgen_binning=PGEN_BINNING,
            pgen_features=feature_identity,
            pgen_model=feature_identity['model'] if feature_identity else _pgen_model_identity(organism, a.locus),
            pgen_source_sha256=sha256(Path(compute_pgen.__code__.co_filename)))
        for i, (task, _) in enumerate(groups):
            selected = reference.filter(pl.col('group_id') == i).rename({'cdr3': 'junction_aa'})
            if not selected.height:
                continue
            for mock in iter_mock_references(selected, events, replicates=a.mock_replicates,
                    seed=stable_seed(task, 'control'), pgen=pgen, species=organism, locus=a.locus, threads=a.threads,
                    prepared_control=prepared_control):
                replicate = int(mock['replicate'][0])
                group = len(populations)
                records.append(mock.select(pl.col('junction_aa').alias('cdr3'), 'v', 'j',
                    pl.lit(group).alias('group_id')))
                populations.append(dict(group_id=group, species=species, locus=a.locus, **task,
                    population='matched_control', replicate=replicate, source_records=selected.height,
                    retained_records=mock.height, excluded_records=0))
    def edge_sink(part):
        for (group,), frame in part.partition_by('group_id', as_dict=True).items():
            path = scratch/f'{group}.tsv'
            exists = path.exists()
            with path.open('ab') as stream:
                frame.drop('group_id').write_csv(stream, separator='\t', include_header=not exists)

    actual_membership_rows = sum(part.height for part in records)
    nodes, edges = reference_graph(pl.concat(records), geometry=a.geometry, radius=a.sw_radius,
        species=organism, locus=a.locus, threads=a.threads, scope=a.scope,
        gap_open=a.gap_open, gap_extend=a.gap_extend, pair_budget=a.pair_budget,
        cdr_loops=loops, edge_sink=edge_sink)
    graph_seconds = time.monotonic()-start
    print(f'reference-coverage {a.geometry} graph prepared: {graph_seconds:.3f}s', flush=True)
    del records, reference
    metadata = pl.DataFrame(populations, infer_schema_length=None)
    summaries, memberships = [], []
    prefix.parent.mkdir(parents=True, exist_ok=True)
    curve_path = Path(str(prefix)+'.rarefaction.tsv')
    written = False
    edges_path = Path(str(prefix)+'.edges.tsv')
    edges.write_csv(edges_path, separator='\t')
    if not nodes['group_id'].is_sorted():
        raise ValueError('coverage graph nodes must be sorted by group_id')
    group_ids = nodes['group_id'].to_numpy()
    group_offsets = np.searchsorted(group_ids, np.arange(len(populations)+1, dtype=group_ids.dtype))
    for context in populations:
        begin, end = group_offsets[context['group_id']:context['group_id']+2]
        selected = nodes.slice(int(begin), int(end-begin))
        if not selected.height:
            summaries.append(context | dict(taxon_resolution='unavailable', coverage_status='no_valid_occurrences'))
            continue
        identity = {k: v for k, v in context.items() if k not in
            ('group_id', 'source_records', 'retained_records', 'excluded_records')}
        modeled = selected.filter(pl.col('geometry_available'))
        unavailable = selected.filter(~pl.col('geometry_available')).with_columns(pl.lit(None, dtype=pl.Int64).alias('community'))
        edge_path = scratch/f"{context['group_id']}.tsv"
        local_edges = (pl.read_csv(edge_path, separator='\t', schema=dict(edges.schema)).sort('a', 'b')
            if edge_path.exists() else edges)
        with edges_path.open('ab') as stream:
            local_edges.write_csv(stream, separator='\t', include_header=False)
        if edge_path.exists():
            edge_path.unlink()
        membership = leiden_taxa(modeled, local_edges, resolution=a.leiden_resolution,
            seed=stable_seed(identity, 'leiden'))
        del local_edges
        memberships.append(pl.concat([membership, unavailable]))
        community_counts = membership.group_by('community').agg(pl.col('occurrences').sum()).sort('community')['occurrences']
        availability = dict(geometry_available_records=int(modeled['occurrences'].sum()),
            geometry_unavailable_records=int(unavailable['occurrences'].sum()),
            geometry_available_clonotypes=modeled.height, geometry_unavailable_clonotypes=unavailable.height)
        for label, counts in (('exact_clonotype', selected['occurrences']), ('leiden_community', community_counts)):
            if not len(counts):
                summaries.append(context | availability | dict(taxon_resolution=label, coverage_status='geometry_unavailable'))
                continue
            summary, curve = abundance_coverage(counts.to_numpy(), nboot=a.nboot,
                seed=stable_seed(identity, label), threads=a.threads)
            summaries.append(context | availability | dict(taxon_resolution=label, **summary))
            curve = curve.with_columns(*[pl.lit(value).alias(key) for key, value in
                (context | dict(taxon_resolution=label)).items()])
            with curve_path.open('ab' if written else 'wb') as stream:
                curve.write_csv(stream, separator='\t', include_header=not written)
            written = True
    from .epitope import write_tables
    tables = [curve_path, edges_path, *write_tables(prefix, epitope_info)]
    if knee_curve is not None:
        knee_path = Path(str(prefix)+'.knee.tsv')
        knee_curve.join(metadata, on='group_id', validate='m:1').write_csv(knee_path, separator='\t')
        tables.append(knee_path)
    for suffix, table in [('.coverage.tsv', pl.DataFrame(summaries, infer_schema_length=None)),
            ('.populations.tsv', metadata), ('.communities.tsv', pl.concat(memberships).join(metadata, on='group_id', validate='m:1'))]:
        path = Path(str(prefix)+suffix)
        table.write_csv(path, separator='\t')
        tables.append(path)
    from vdjtools.stats import rarefaction_batch
    import vdjtools._core as stats_native
    import seqtree._core as search_native
    from ..match.geometry import load_v_loops
    manifest = dict(reference=db.provenance(a.vdjdb), reference_population=a.reference_population, epitope=a.epitope,
        **epitope_info,
        shared_invocation_inputs=inputs is not None,
        targets_sha256=sha256(Path(a.targets_from_sample)) if a.targets_from_sample else None,
        species=species, locus=a.locus, geometry=a.geometry, scope=a.scope if a.geometry == 'levenshtein' else None,
        sw_radius=a.sw_radius, sw_radius_selection=('reference-only within-full-pMHC nearest-other-clonotype distances pooled over selected species/chain; convex Kneedle S=1 strongest confirmed interior knee; no smoothing' if a.sw_radius_knee else 'explicit'), cdr_loops=a.cdr_loops if a.geometry == 'sw' else None,
        gap_open=a.gap_open if a.geometry == 'sw' else None, gap_extend=a.gap_extend if a.geometry == 'sw' else None,
        v_model_sha256=sha256(Path(load_v_loops.__code__.co_filename).parents[1]/'resources/vloops'/f'{organism}_v_loops.tsv'),
        sw_distance='sum_loop [SW(s,s)+SW(t,t)-2*SW(s,t)]; raw BLOSUM62; full junction including anchors; selected V loops stripped of alignment dots; no trim/weights/kernel',
        clonotype_key='full junction + supplied normalized V call + supplied normalized J call; alleles retained; missing calls remain explicit',
        v_model_policy='existing source-aware allele resolver, then flagged fixed first complete co-locus V-loop template for unresolved geometry; raw junction/V/J keys and occurrence counts retained; lexical template order; no target labels',
        leiden=dict(objective='RBConfiguration modularity', resolution=a.leiden_resolution, iterations='until stable', edges='unweighted; undirected; isolates retained'),
        sampling_unit='VDJdb source-record occurrence; community abundance sums member-clonotype occurrences; control draws use matched event slots',
        estimator='vdjtools iNEXT finite-sample Chao1; f2>0: S+(n-1)/n*f1^2/(2*f2); f2=0: S+(n-1)/n*f1*(f1-1)/2',
        coverage='vdjtools sample_coverage; unseen mass=1-C; unseen-type fraction=(estimated-observed)/estimated',
        uncertainty='vdjtools augmented-assemblage bootstrap; richness normal 95% CI, coverage percentile 95% envelope; partitions fixed at full observed graph; no assertion of unseen-family exhaustion',
        bootstrap_probability_api='vdjtools.stats.bootstrap_abundances; exact compressed unseen-taxon law owned by vdjtools',
        control=control_info, nboot=a.nboot, seed=a.seed, seed_policy='SHA256(base seed,species,locus,task/population,resolution); first 4 bytes little-endian modulo 2**31',
        membership_budget_bytes=membership_budget, estimated_membership_rows=estimated_membership_rows,
        estimated_membership_bytes=estimated_membership_bytes, actual_membership_rows=actual_membership_rows,
        membership_unit='raw reference/mock source-record rows before graph-key deduplication',
        threads=a.threads, pair_budget=a.pair_budget, search_batches=1+int(knee_curve is not None),
        source_sha256=sha256(Path(__file__)), reducer_sha256=sha256(Path(reference_graph.__code__.co_filename)),
        stats_source_sha256=sha256(Path(rarefaction_batch.__code__.co_filename)),
        stats_native_sha256=sha256(Path(stats_native.__file__)), search_native_sha256=sha256(Path(search_native.__file__)),
        software={name:version(name) for name in ('vdjmatch','seqtree','vdjtools','polars','numpy','igraph','leidenalg')},
        edge_storage='bounded native chunks; temporary per-population files; one population loaded for Leiden; final edges sorted by node IDs',
        graph_preparation_seconds=graph_seconds, tables={path.name:sha256(path) for path in tables}, seconds=time.monotonic()-start)
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return 0


def validate_options(a):
    """Reject invalid/unsupported graph requests before any raw computation."""
    membership_budget = getattr(a, 'membership_budget_bytes', 2*1024**3)
    if not 1 <= membership_budget <= 8*1024**3:
        raise ValueError('membership budget must be positive and at most 8 GiB')
    if a.threads < 1 or a.nboot < 2 or a.seed < 0 or a.mock_replicates < 0 or a.mock_replicates == 1:
        raise ValueError('positive threads, nboot>=2, seed>=0 and mock-replicates zero or >=2 required')
    if getattr(a, 'control_model', None) and not a.control:
        raise ValueError('--control-model requires its original --control input for validation')
    if bool(a.control) != bool(a.mock_replicates):
        raise ValueError('--control and --mock-replicates>=2 must be supplied together')
    if a.pgen_features and (not a.control or not a.mock_replicates):
        raise ValueError('--pgen-features requires --control and --mock-replicates>=2')
    if a.geometry == 'sw' and a.sw_radius is None and not a.sw_radius_knee:
        raise ValueError('full SW requires --sw-radius or --sw-radius-knee')
    if a.sw_radius is not None and a.sw_radius_knee:
        raise ValueError('select explicit radius or reference Kneedle, not both')
    if a.geometry == 'levenshtein' and (a.sw_radius is not None or a.sw_radius_knee):
        raise ValueError('--sw-radius applies only to SW geometry')
    if a.geometry == 'sw' and a.scope != '1,1,1,1':
        raise ValueError('--scope applies only to Levenshtein geometry')
    try:
        import igraph
        import leidenalg
    except ImportError as exc:
        raise ValueError('Leiden requires vdjmatch[diagnostics]') from exc
    import math
    if a.pair_budget < 1 or a.gap_open < 0 or a.gap_extend < 0 or not math.isfinite(a.leiden_resolution) or a.leiden_resolution <= 0:
        raise ValueError('positive pair budget/resolution and nonnegative gap costs required')
    if a.sw_radius is not None and not 0 <= a.sw_radius < 2**28:
        raise ValueError('SW radius must be nonnegative and below 2**28')
    if a.geometry == 'levenshtein':
        from ..match.scope import parse_scope
        _, insertions, deletions, _ = parse_scope(a.scope)
        if insertions != deletions:
            raise ValueError('undirected edit graphs require symmetric insertion/deletion budgets')
