"""Reference-only graph, junction-length and V/J-usage composition."""
from pathlib import Path


def register(subparsers):
    p = subparsers.add_parser('reference-structure', help=__doc__)
    p.add_argument('--vdjdb', required=True)
    p.add_argument('--species', default='HomoSapiens')
    p.add_argument('--locus', choices=['TRA','TRB'], required=True)
    p.add_argument('--epitope', action='append')
    p.add_argument('--targets-from-sample', help='compatible pMHC unions declared in an AIRR table')
    p.add_argument('--reference-population', choices=['curated-benchmark','all-valid'], default='curated-benchmark')
    p.add_argument('--control', help='raw rearrangement control for descriptive V/J event-frequency ratios')
    p.add_argument('--motif-mode', choices=['exact3', 'imgt3', 'wild4', 'imgt_wild4'],
                   default='imgt_wild4', help='control-relative motif projection (default: imgt_wild4)')
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
    from ..aggregate.motif_information import motif_population, motif_selection_information
    from ..aggregate.germline_information import germline_projection_map, iter_germline_information
    from ..aggregate.structure import composition, graph_profiles, reference_gene_usage, control_gene_usage, compare_gene_usage, gene_selection_information
    from ..cluster import overlap
    from ..db.cache import sha256
    from ..evalue.control import _organism
    from ..match.geometry import load_v_loops, resolve_v_alleles
    from ..match.motifs import masked_motif_presence_counts
    from ..evalue.epitope_rules import _tokens
    from ..runner.annotation import _validate_output_paths

    start = time.monotonic()
    if a.threads < 1:
        raise ValueError('positive threads required')
    if a.motif_mode != 'imgt_wild4' and not a.control:
        raise ValueError('--motif-mode requires --control')
    protected = [a.vdjdb]+([a.targets_from_sample] if a.targets_from_sample else [])+([a.control] if a.control else [])
    prefix = _validate_output_paths(a.output_prefix, protected)
    organism = _organism(a.species)
    species = {'human':'HomoSapiens','mouse':'MusMusculus'}[organism]
    control_usage, control_info, control_motifs = None, None, None
    motif_path, motif_rows, motif_written = None, 0, False
    if a.control:
        from .combined import _raw_control
        from ..evalue.mock_reference import control_events
        control_raw, control_info = _raw_control(a.control)
        events = control_events(control_raw, species=organism, locus=a.locus)
        control_usage = control_gene_usage(events, locus=a.locus)
        control_motifs = motif_population(events.rename({'junction_aa': 'cdr3'}), locus=a.locus, mode=a.motif_mode)
        prefix.parent.mkdir(parents=True, exist_ok=True)
        motif_path = Path(str(prefix)+'.motifs.tsv')
        control_info.update(path=str(Path(a.control).resolve()), sha256=sha256(Path(a.control)),
            event_keys=events.height, total_events=int(events['event_weight'].sum()),
            reader_sha256=sha256(Path(_raw_control.__code__.co_filename)),
            event_normalizer_sha256=sha256(Path(control_events.__code__.co_filename)))
        del control_raw, events
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
        mixture_full, _ = graph_mixture_weights(mixture_records, threads=a.threads, hamming_edges=edge_names, include_punctures=False)
    else:
        mixture_full = None
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
    degree_nodes = {key[0]: part.sort('_local') for key, part in
                    degree_membership.group_by('_group')} if degree_membership is not None else {}
    summaries, degree_rows, length_tables, gene_tables = [], [], [], []
    for group_id, (task, source) in enumerate(groups):
        retained = source.filter(available)
        junctions = memberships[group_id]
        local_edges = grouped_edges[group_id].iter_rows() if group_id in grouped_edges else []
        full, _ = graph_profiles(junctions.height,local_edges,include_punctures=False)
        v_full, _ = composition(retained,junctions,include_punctures=False)
        context = dict(_mixture_group=group_id,species=species,locus=a.locus,**task,source_observations=source.height,
                       unavailable_records=source.height-retained.height)
        gene_context = {k: v for k, v in context.items() if k != '_mixture_group'}
        gene_tables.append(reference_gene_usage(source.filter(mixture_available), locus=a.locus).with_columns(
            *[pl.lit(value).alias(key) for key, value in gene_context.items()]))
        if control_motifs is not None:
            motifs = motif_selection_information(source.filter(mixture_available), control_motifs, locus=a.locus, mode=a.motif_mode).with_columns(
                *[pl.lit(value).alias(key) for key, value in gene_context.items()])
            with motif_path.open('ab' if motif_written else 'wb') as stream:
                motifs.write_csv(stream, separator='\t', include_header=not motif_written)
            motif_written = True
            motif_rows += motifs.height
            del motifs
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
        # Local node indices, not join order, align junction lengths with degrees.
        length_graph = degree_nodes[group_id].select(
            pl.col('cdr3').str.len_chars().alias('junction_length')).with_columns(
                pl.Series('degree', degrees)).group_by('junction_length').agg(
                    pl.len().alias('junctions'),
                    (pl.col('degree').sum() // 2).alias('graph_edges'),
                    (pl.col('degree') == 0).sum().alias('isolated_junctions'),
                    pl.col('degree').mean().alias('mean_degree'),
                    (pl.col('degree') == 0).mean().alias('isolated_fraction')) if n else None
        length_counts = source.filter(mixture_available).with_columns(
            pl.col('cdr3').str.len_chars().alias('junction_length'), available.alias('_modeled')
        ).group_by('junction_length').agg(pl.len().alias('valid_reference_records'),
                                        pl.col('_modeled').sum().alias('modeled_reference_records'))
        length_context = {k: v for k, v in degree_context.items() if k not in ('junctions', 'graph_edges')}
        if n:
            length_table = length_graph.join(length_counts, on='junction_length', validate='1:1')
        else:
            length_table = pl.DataFrame(dict(junction_length=[None], junctions=[0], graph_edges=[0],
                isolated_junctions=[0], mean_degree=[None], isolated_fraction=[None],
                valid_reference_records=[0], modeled_reference_records=[0]))
        length_tables.append(length_table.with_columns(
            *[pl.lit(value).alias(key) for key, value in length_context.items()]))
        degree_summary = dict(hamming_junctions=n, hamming_graph_edges=edge_count,
            hamming_mean_degree=float(degrees.mean()) if n else None,
            hamming_degree_variance=float(degrees.var()) if n else None,
            hamming_max_degree=int(degrees.max()) if n else None,
            hamming_isolated_fraction=float(np.count_nonzero(degrees == 0) / n) if n else None)
        degree_summary.update({f'hamming_{name}_degree': float(np.quantile(degrees, q)) if n else None
                               for name, q in [('q10', .1), ('median', .5), ('q90', .9)]})
        summaries.append(context | full | v_full | degree_summary | dict(structure_status='available' if full['junctions'] else 'empty_reference'))
    prefix.parent.mkdir(parents=True,exist_ok=True)
    summary_path = Path(str(prefix)+'.structures.tsv')
    degree_path = Path(str(prefix)+'.degrees.tsv')
    length_path = Path(str(prefix)+'.lengths.tsv')
    gene_path = Path(str(prefix)+'.genes.tsv')
    if not summaries:
        raise ValueError('no declared reference targets')
    summary = pl.DataFrame(summaries, infer_schema_length=None).with_columns(pl.col('_mixture_group').cast(pl.UInt32))
    if mixture_full is not None:
        weights = {'original_weight', 'gapblock_weight', 'stem4_weight'}
        def names(frame):
            return frame.rename({c: '_mixture_group' if c == 'group_id' else
                c if c in weights else 'mixture_' + c for c in frame.columns})
        summary = summary.join(names(mixture_full), on='_mixture_group', how='left', validate='1:1', maintain_order='left')
    else:
        summary = summary.with_columns(*[pl.lit(None, dtype=pl.Float64).alias(c)
            for c in ('original_weight', 'gapblock_weight', 'stem4_weight')],
            pl.lit('empty_reference').alias('mixture_weight_status'))
    summary.drop('_mixture_group').write_csv(summary_path, separator='\t')
    tables = [summary_path, degree_path, length_path, gene_path]
    gene_usage = pl.concat(gene_tables, how='vertical_relaxed')
    germline_info = None
    if control_usage is not None:
        gene_usage = compare_gene_usage(gene_usage, control_usage)
        gene_usage, information = gene_selection_information(gene_usage)
        information_path = Path(str(prefix)+'.information.tsv')
        information.write_csv(information_path, separator='\t')
        tables.append(information_path)
        control_path = Path(str(prefix)+'.control-genes.tsv')
        control_usage.with_columns(pl.lit(species).alias('species'), pl.lit(a.locus).alias('locus')).write_csv(
            control_path, separator='\t')
        tables.extend([control_path, motif_path])
        mapping = germline_projection_map(model, locus=a.locus)
        mapping_path = Path(str(prefix)+'.germline-map.tsv')
        factor_path = Path(str(prefix)+'.germline.tsv')
        information_path = Path(str(prefix)+'.germline-information.tsv')
        mapping.write_csv(mapping_path, separator='\t')
        factor_rows, information_rows = 0, 0
        for number, (factors, information) in enumerate(iter_germline_information(gene_usage, control_usage, mapping)):
            for path, frame in ((factor_path, factors), (information_path, information)):
                with path.open('ab' if number else 'wb') as stream:
                    frame.write_csv(stream, separator='\t', include_header=number == 0)
            factor_rows += factors.height
            information_rows += information.height
        tables.extend([mapping_path, factor_path, information_path])
        germline_info = dict(factor_rows=factor_rows, information_rows=information_rows,
            sampling_unit='original distinct full-junction/V-gene keys; no deduplication after residue projection',
            projections='one-based aligned cdr1/cdr2/cdr25 columns; full AA and 11 IMGT classes separately',
            allele_policy='all bundled alleles of a gene must agree per column; no default allele or inferred allele frequencies',
            missing_policy='missing gene/model/loop or allele-dependent column keeps all mass and unavailable comparison; no smoothing',
            gap_policy='dot is aligned absence within a populated loop; an entirely empty loop template is unavailable',
            interpretation='consensus bits/KL describe gene-consensus categories; residue_bits_lower allows every compatible or unresolved control allele; correlated positions are not additive',
            reducer_sha256=sha256(Path(iter_germline_information.__code__.co_filename)))
    gene_usage.write_csv(gene_path, separator='\t')
    pl.DataFrame(degree_rows, schema_overrides={'degree': pl.Int64, 'junction_count': pl.Int64},
        infer_schema_length=None).write_csv(degree_path, separator='\t')
    pl.concat(length_tables, how='vertical_relaxed').sort('species', 'locus', *PMHC, 'junction_length').write_csv(
        length_path, separator='\t')
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
        mixture_weights='gap=1-sum_L p_L^2; stem4=sum_LV p_LV^2-sum_LVC p_LVC^2; original=remainder',
        mixture_population='all valid full AA junction/V-gene keys; unknown or ambiguous V is explicit category; source records counted separately',
        modeled_population='existing resolved V and junction length >=8; unchanged modeled feature columns',
        mixture_source_sha256=sha256(Path(graph_mixture_weights.__code__.co_filename)),
        graph_scope='full-junction Hamming<=1; no indels',
        graph_node_unit='unique full junction', v_usage_unit='unique full junction/resolved V allele keys; gene-level categories',
        degree_population='all valid full junctions, including unresolved V; fixed full reference',
        degree_quantiles='linear interpolation of the unique-junction degree distribution',
        length_coverage='exact junction length; all valid reference records and unique junctions; modeled records counted separately; no cross-length Hamming edges',
        gene_usage='all valid reference rows; V/J alleles pooled; unresolved calls retained; full and exact-length strata; record and junction/gene-key denominators explicit',
        control=control_info,
        gene_control_comparison=('reference record/key fractions divided by control event fractions; no smoothing; unresolved genes, absent control genes and absent control lengths explicit; descriptive composition only'
            if control_info is not None else None),
        motif_information=(dict(rows=motif_rows, mode=a.motif_mode,
            k=4 if 'wild' in a.motif_mode else 3, masks=[0,1,2,3] if 'wild' in a.motif_mode else [],
            alphabet='11 IMGT classes' if 'imgt' in a.motif_mode else '20 amino acids', token_scope='full junction; once-per-key presence',
            sampling_unit='distinct full-junction/V-gene keys on both sides',
            scopes='within V, pooled lengths and exact length',
            statistic='presence log2 likelihood ratio and binary presence/absence KL; overlapping tokens are not independent',
            joint_projection='stratum_bits + presence_bits = joint_presence_bits; whole-population junction/V-gene totals; one token projection only',
            row_policy='reference-observed tokens only; no automatic scorer/threshold; missing control support stays unavailable',
            reference_keys_without_fourmer=int(gene_usage.filter((pl.col("segment")=="V") &
                (pl.col("length_scope")=="exact_length") & (pl.col("junction_length")<4))["junction_gene_keys"].sum()),
            reducer_sha256=sha256(Path(motif_selection_information.__code__.co_filename)),
            counter_sha256=sha256(Path(masked_motif_presence_counts.__code__.co_filename)),
            tokenizer_sha256=sha256(Path(_tokens.__code__.co_filename))) if control_info is not None else None),
        selection_information=(dict(sampling_unit='distinct full-junction/segment-gene keys in both populations; V and J separate',
            units='log2 likelihood-ratio bits per junction/gene key',
            decomposition='I(length,gene)=I(length)+I(gene|length); marginal gene information also reported',
            missing_support='all original mass retained; incomplete expectations null; no smoothing or renormalization',
            reducer_sha256=sha256(Path(gene_selection_information.__code__.co_filename))) if control_info is not None else None),
        germline_information=germline_info,
        entropy_units='natural logarithms', length_units='amino-acid residues, conserved anchors included',
        exclusion='none: model composition and weights use the fixed full reference; benchmark exclusion applies only to hit lists',
        interpretation='reference-only morphology and fixed graph-mixture heuristic; no label fit, recognition probability or P/E-value',
        source_sha256=sha256(Path(__file__)), reducer_sha256=sha256(Path(graph_profiles.__code__.co_filename)),
        tables={str(p.name):sha256(p) for p in tables}, seconds=time.monotonic()-start),indent=2)+'\n')
    return 0
