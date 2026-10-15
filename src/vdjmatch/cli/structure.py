"""Reference-only graph, junction-length and V-usage composition."""
from pathlib import Path


def register(subparsers):
    p = subparsers.add_parser('reference-structure', help=__doc__)
    p.add_argument('--vdjdb', required=True)
    p.add_argument('--species', default='HomoSapiens')
    p.add_argument('--locus', choices=['TRA','TRB'], required=True)
    p.add_argument('--epitope', action='append')
    p.add_argument('--targets-from-sample', help='compatible pMHC unions declared in an AIRR table')
    p.add_argument('--reference-population', choices=['curated-benchmark','all-valid'], default='curated-benchmark')
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--output-prefix', required=True)
    p.set_defaults(func=main)


def main(a):
    import json
    import time
    import numpy as np
    import polars as pl
    from importlib.metadata import version
    import seqtree._core as native_core
    from .. import db
    from ..aggregate.candidates import PMHC
    from ..aggregate.structure import composition, graph_profiles
    from ..cluster import overlap
    from ..db.cache import sha256
    from ..evalue.control import _organism
    from ..match.geometry import load_v_loops, resolve_v_alleles
    from ..runner.annotation import _validate_output_paths

    start = time.monotonic()
    if a.threads < 1:
        raise ValueError('positive threads required')
    protected = [a.vdjdb]+([a.targets_from_sample] if a.targets_from_sample else [])
    prefix = _validate_output_paths(a.output_prefix, protected)
    organism = _organism(a.species)
    species = {'human':'HomoSapiens','mouse':'MusMusculus'}[organism]
    raw = db.load(a.vdjdb, species=species, gene=a.locus, epitope=a.epitope)
    if a.reference_population == 'all-valid':
        raw = raw.with_columns(pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False).alias('reference_valid'))
    model = load_v_loops(organism)
    raw = raw.with_columns(resolve_v_alleles(raw['v'],model,locus=a.locus).alias('_allele'))
    available = pl.col('reference_valid') & pl.col('_allele').is_not_null() & (pl.col('cdr3').str.len_chars() >= 8)
    if a.targets_from_sample:
        from ..io.airr import _read_table
        from .historical import assayed_targets
        targets = _read_table(a.targets_from_sample)
        if a.epitope:
            targets = targets.filter(pl.col('epitope').is_in(a.epitope))
        if 'query_id' not in targets.columns:
            targets = targets.with_row_index('query_id')
        groups = [(task, selected) for task, _, selected in
            assayed_targets(targets.select('query_id'),targets,raw,a.species)]
    else:
        groups = [(dict(zip(PMHC,key)), part) for key,part in raw.group_by(PMHC,maintain_order=True)]
    # One batched native search for the whole reference population, then subset its
    # edges by declared target membership; no per-epitope native search loop.
    mixture_available = pl.col('reference_valid') & pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)
    selected_nodes = [source.filter(mixture_available).select('cdr3') for _,source in groups]
    nodes = (pl.concat(selected_nodes).unique(maintain_order=True) if selected_nodes else raw.head(0).select('cdr3')).with_row_index('_node')
    pairs = overlap(nodes['cdr3'].to_list(),scope='1,0,0,1',threads=a.threads) if nodes.height else None
    edges = pairs.filter(pl.col('a_idx')<pl.col('b_idx')).select('a_idx','b_idx') if pairs is not None else None
    memberships = []
    for group_id, (_, source) in enumerate(groups):
        memberships.append(source.filter(available).select('cdr3').unique(maintain_order=True)
            .with_row_index('_local').join(nodes,on='cdr3',validate='1:1')
            .sort('_local')
            .with_columns(pl.lit(group_id).alias('_group')))
    membership = pl.concat(memberships) if memberships else None
    grouped_edges = {}
    if membership is not None and edges is not None and edges.height:
        left = membership.select(pl.col('_node').cast(pl.Int64).alias('a_idx'),'_group',pl.col('_local').alias('a'))
        right = membership.select(pl.col('_node').cast(pl.Int64).alias('b_idx'),'_group',pl.col('_local').alias('b'))
        local = edges.join(left,on='a_idx',how='inner',validate='m:m').join(
            right,on=['b_idx','_group'],how='inner',validate='m:1')
        grouped_edges = {key[0]: frame.select('a','b') for key, frame in local.group_by('_group')}
    from ..aggregate.reference_combination import graph_mixture_weights
    mixture_parts = [source.filter(mixture_available).select('cdr3', 'v').with_columns(
        pl.lit(i, dtype=pl.UInt32).alias('group_id')) for i, (_, source) in enumerate(groups)]
    mixture_records = pl.concat(mixture_parts) if mixture_parts else None
    if mixture_records is not None and mixture_records.height:
        edge_names = edges.join(nodes.select(pl.col('_node').cast(pl.Int64).alias('a_idx'),
            pl.col('cdr3').alias('cdr3_a')), on='a_idx').join(nodes.select(
                pl.col('_node').cast(pl.Int64).alias('b_idx'), pl.col('cdr3').alias('cdr3_b')), on='b_idx').select('cdr3_a', 'cdr3_b')
        mixture_full, mixture_punctures = graph_mixture_weights(mixture_records, threads=a.threads, hamming_edges=edge_names)
    else:
        mixture_full = mixture_punctures = None
    # Degree distributions use every valid full junction, including those without
    # a resolved V model. Reuse the same native edges as the scoring mixture.
    degree_membership = pl.concat([
        source.select('cdr3').unique(maintain_order=True).with_row_index('_local')
        .join(nodes, on='cdr3', validate='1:1').with_columns(pl.lit(i).alias('_group'))
        for i, source in enumerate(selected_nodes)]) if selected_nodes else None
    degree_edges = {}
    if degree_membership is not None and edges is not None and edges.height:
        left = degree_membership.select(pl.col('_node').cast(pl.Int64).alias('a_idx'),
            '_group', pl.col('_local').alias('a'))
        right = degree_membership.select(pl.col('_node').cast(pl.Int64).alias('b_idx'),
            '_group', pl.col('_local').alias('b'))
        joined = edges.join(left, on='a_idx', validate='m:m').join(
            right, on=['b_idx', '_group'], validate='m:1')
        degree_edges = {key[0]: part.select('a', 'b') for key, part in joined.group_by('_group')}
    summaries, punctures, degree_rows = [], [], []
    for group_id, (task, source) in enumerate(groups):
        retained = source.filter(available)
        junctions = memberships[group_id]
        local_edges = grouped_edges[group_id].iter_rows() if group_id in grouped_edges else []
        full, removed = graph_profiles(junctions.height,local_edges)
        v_full, v_removed = composition(retained,junctions)
        context = dict(_mixture_group=group_id,species=species,locus=a.locus,**task,source_observations=source.height,
                       unavailable_records=source.height-retained.height)
        n = selected_nodes[group_id].unique().height
        endpoints = degree_edges[group_id].to_numpy().ravel() if group_id in degree_edges else np.array([], dtype=np.int64)
        degrees = np.bincount(endpoints, minlength=n)
        values, counts = np.unique(degrees, return_counts=True)
        edge_count = len(endpoints) // 2
        assert int(degrees.sum()) == 2 * edge_count
        degree_context = {k: v for k, v in context.items() if k != '_mixture_group'}
        degree_context['modeled_unavailable_records'] = degree_context.pop('unavailable_records')
        degree_context['invalid_reference_records'] = source.height - source.filter(mixture_available).height
        degree_context.update(graph_population='all_valid_junctions', graph_scope='full_reference',
                              junctions=n, graph_edges=edge_count)
        degree_rows.extend(degree_context | dict(degree=int(d), junction_count=int(c))
                           for d, c in zip(values, counts))
        if not n:
            degree_rows.append(degree_context | dict(degree=None, junction_count=0))
        degree_summary = dict(hamming_junctions=n, hamming_graph_edges=edge_count,
            hamming_mean_degree=float(degrees.mean()) if n else None,
            hamming_degree_variance=float(degrees.var()) if n else None,
            hamming_max_degree=int(degrees.max()) if n else None,
            hamming_isolated_fraction=float(np.count_nonzero(degrees == 0) / n) if n else None)
        degree_summary.update({f'hamming_{name}_degree': float(np.quantile(degrees, q)) if n else None
                               for name, q in [('q10', .1), ('median', .5), ('q90', .9)]})
        summaries.append(context | full | v_full | degree_summary | dict(structure_status='available' if full['junctions'] else 'empty_reference'))
        for sequence, graph, usage in zip(junctions['cdr3'],removed,v_removed):
            punctures.append(context | dict(excluded_junction=sequence, structure_status='available' if graph['junctions'] else 'empty_after_puncture') | graph | usage)
    prefix.parent.mkdir(parents=True,exist_ok=True)
    summary_path = Path(str(prefix)+'.structures.tsv')
    puncture_path = Path(str(prefix)+'.punctures.tsv')
    degree_path = Path(str(prefix)+'.degrees.tsv')
    if not summaries:
        raise ValueError('no declared reference targets')
    summary = pl.DataFrame(summaries, infer_schema_length=None).with_columns(pl.col('_mixture_group').cast(pl.UInt32))
    puncture_table = (pl.DataFrame(punctures, infer_schema_length=None).with_columns(
        pl.col('_mixture_group').cast(pl.UInt32)) if punctures else
        summary.head(0).with_columns(pl.lit(None, dtype=pl.String).alias('excluded_junction')))
    mixture_new_punctures = 0
    if mixture_full is not None:
        weights = {'original_weight', 'gapblock_weight', 'stem4_weight'}
        def names(frame):
            return frame.rename({c: '_mixture_group' if c == 'group_id' else
                c if c in weights or c == 'excluded_junction' else 'mixture_' + c for c in frame.columns})
        summary = summary.join(names(mixture_full), on='_mixture_group', how='left', validate='1:1', maintain_order='left')
        context_columns = ['species', 'locus', *PMHC, 'source_observations', 'unavailable_records']
        puncture_table = puncture_table.with_columns(pl.lit(True).alias('_modeled_profile')).join(
            names(mixture_punctures), on=['_mixture_group', 'excluded_junction'], how='full', coalesce=True,
            validate='1:1').join(summary.select('_mixture_group', *context_columns), on='_mixture_group',
                how='left', suffix='_context', validate='m:1').with_columns(*[
                    pl.coalesce(c, c + '_context').alias(c) for c in context_columns],
                    pl.col('_modeled_profile').is_null().alias('mixture_only_puncture')).drop(
                        '_modeled_profile', *[c + '_context' for c in context_columns])
        mixture_new_punctures = puncture_table['mixture_only_puncture'].sum()
    else:
        summary = summary.with_columns(*[pl.lit(None, dtype=pl.Float64).alias(c)
            for c in ('original_weight', 'gapblock_weight', 'stem4_weight')],
            pl.lit('empty_reference').alias('mixture_weight_status'))
        puncture_table = puncture_table.with_columns(*[pl.lit(None, dtype=pl.Float64).alias(c)
            for c in ('original_weight', 'gapblock_weight', 'stem4_weight')],
            pl.lit('empty_reference').alias('mixture_weight_status'),
            pl.lit(False).alias('mixture_only_puncture'))
    summary.drop('_mixture_group').write_csv(summary_path, separator='\t')
    puncture_table.drop('_mixture_group').write_csv(puncture_path, separator='\t')
    pl.DataFrame(degree_rows, schema_overrides={'degree': pl.Int64, 'junction_count': pl.Int64},
        infer_schema_length=None).write_csv(degree_path, separator='\t')
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(dict(reference=db.provenance(a.vdjdb),
        reference_population=a.reference_population, species=species,locus=a.locus,threads=a.threads,
        epitope=a.epitope, targets_sha256=sha256(Path(a.targets_from_sample)) if a.targets_from_sample else None,
        software={name:version(name) for name in ('vdjmatch','seqtree','polars','numpy')},
        native_sha256=sha256(Path(native_core.__file__)),
        v_model_sha256=sha256(Path(load_v_loops.__code__.co_filename).parents[1]/'resources/vloops'/f'{organism}_v_loops.tsv'),
        loader_sha256=sha256(Path(db.load.__code__.co_filename)),
        search_source_sha256=sha256(Path(overlap.__code__.co_filename)),
        native_batches=1 if nodes.height else 0, searched_unique_junctions=nodes.height,
        global_hamming_edges=edges.height if edges is not None else 0, targets=len(summaries),
        puncture_rows=puncture_table.height, mixture_only_punctures=mixture_new_punctures,
        mixture_weights='gap=1-sum_L p_L^2; stem4=sum_LV p_LV^2-sum_LVC p_LVC^2; original=remainder',
        mixture_population='all valid full AA junction/V-gene keys; unknown or ambiguous V is explicit category; source records counted separately',
        modeled_population='existing resolved V and junction length >=8; unchanged modeled feature columns',
        mixture_source_sha256=sha256(Path(graph_mixture_weights.__code__.co_filename)),
        graph_scope='full-junction Hamming<=1; no indels',
        graph_node_unit='unique full junction', v_usage_unit='unique full junction/resolved V allele keys; gene-level categories',
        degree_population='all valid full junctions, including unresolved V; full reference before query puncture',
        degree_quantiles='linear interpolation of the unique-junction degree distribution',
        entropy_units='natural logarithms', length_units='amino-acid residues, conserved anchors included',
        exclusion='puncture removes the complete full junction, all V annotations and incident edges; articulation fragments recomputed exactly',
        interpretation='reference-only morphology and fixed graph-mixture heuristic; no label fit, recognition probability or P/E-value',
        source_sha256=sha256(Path(__file__)), reducer_sha256=sha256(Path(graph_profiles.__code__.co_filename)),
        tables={str(p.name):sha256(p) for p in (summary_path,puncture_path,degree_path)}, seconds=time.monotonic()-start),indent=2)+'\n')
    return 0
