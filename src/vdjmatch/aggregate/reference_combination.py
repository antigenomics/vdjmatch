"""Reference graph mixture and matched three-channel null normalization.

The fixed mixture partitions two uniform junction/V-gene key draws. Different
lengths select the gap-block channel; same length/V but different Hamming-one
components select stem4; the remainder selects original evidence. This is an
explicit structural heuristic, not a learned recognition model.
"""
from __future__ import annotations

import numpy as np
import polars as pl

from .combination import THREE_WEIGHTS, null_normalized_evidence
from .structure import graph_profiles, composition


def _component_colour_squares(n, edges, colours, *, include_punctures=True):
    """Full and every vertex-punctured component/V square sums in two DFS passes.

    ``colours`` holds one count for each unique junction/V key. Low-link
    articulation fragments inherit subtree colour histograms. Work is
    O((vertices+edges)*V), space O(vertices*V+edges), where V is the observed
    V-category count; no vertex-specific BFS or scientific score recomputation.
    """
    neighbours = [[] for _ in range(n)]
    for a, b in edges:
        neighbours[a].append(b)
        neighbours[b].append(a)
    parent = np.full(n, -1, dtype=np.int64)
    seen = np.full(n, -1, dtype=np.int64)
    low = np.zeros(n, dtype=np.int64)
    component = np.zeros(n, dtype=np.int64)
    subtree = colours.copy()
    separated = np.zeros_like(colours)
    fragment_squares = np.zeros(n, dtype=np.int64)
    totals = []
    clock = 0
    for root in range(n):
        if seen[root] >= 0:
            continue
        seen[root] = low[root] = clock
        clock += 1
        component[root] = len(totals)
        stack = [(root, iter(neighbours[root]))]
        while stack:
            u, children = stack[-1]
            v = next(children, None)
            if v is not None:
                if v == parent[u]:
                    continue
                if seen[v] < 0:
                    parent[v] = u
                    component[v] = component[u]
                    seen[v] = low[v] = clock
                    clock += 1
                    stack.append((v, iter(neighbours[v])))
                else:
                    low[u] = min(low[u], seen[v])
            else:
                stack.pop()
                p = parent[u]
                if p >= 0:
                    subtree[p] += subtree[u]
                    low[p] = min(low[p], low[u])
                    if include_punctures and low[u] >= seen[p]:
                        separated[p] += subtree[u]
                        fragment_squares[p] += subtree[u] @ subtree[u]
        totals.append(subtree[root].copy())
    total = np.asarray(totals, dtype=np.int64)
    component_square = np.einsum('ij,ij->i', total, total)
    full = int(component_square.sum())
    if not include_punctures:
        return full, None
    remainder = total[component] - colours - separated
    punctured = full - component_square[component] + fragment_squares + np.einsum('ij,ij->i', remainder, remainder)
    return full, punctured


def _weights(length_square, length_v_square, component_v_square, size):
    if not size:
        return dict(original_weight=None, gapblock_weight=None, stem4_weight=None,
                    weight_status='empty_reference')
    denominator = size * size
    # Exact integer partitions keep nonnegativity without fitted coefficients.
    return dict(original_weight=(length_square - length_v_square + component_v_square) / denominator,
                gapblock_weight=(denominator - length_square) / denominator,
                stem4_weight=(length_v_square - component_v_square) / denominator,
                weight_status='available')


def graph_mixture_weights(reference: pl.DataFrame, *, threads: int = 1, hamming_edges: pl.DataFrame | None = None,
                          include_punctures: bool = True):
    """Return full and exact-junction-punctured graph-mixture tables.

    Internal reference rows contain ``group_id,cdr3,v``; duplicates retain source
    observation counts. Mixture mass uses unique ``(full junction,V gene)`` keys,
    with alleles pooled. Missing/ambiguous V calls are an explicit unknown category.
    Every complete junction and all its V/source rows disappear in its puncture.
    One global Hamming-one native batch supplies all group graphs, unless a
    transient union ``hamming_edges`` table (``cdr3_a,cdr3_b``) is supplied for
    reuse across randomized models. No disk graph cache is read or written. Full and
    articulation-aware punctured graph/length/V diagnostics accompany weights.

    For key probabilities p, w_gap=1-sum_L p_L^2;
    w_stem=sum_(L,V) p_LV^2-sum_(L,V,C) p_LVC^2;
    w_original=1-w_gap-w_stem. L is length, V gene and C Hamming-one component.
    Components have a single length, making these mutually exclusive exhaustive
    pair events. The law has no label fit, task switch or adjustable coefficient.
    Recompute it for every mock model. Exact-hit exclusion never changes the
    full model. Optional punctures are retained only as legacy diagnostics.
    Empty diagnostic punctures have unavailable weights. Normal scoring uses
    include_punctures=False: fixed full-model weights and no query deletions.
    """
    from ..cluster import overlap

    if threads < 1 or not {'group_id', 'cdr3', 'v'}.issubset(reference.columns) or not reference.height:
        raise ValueError('graph mixture requires nonempty group_id/cdr3/v reference and positive threads')
    if reference['group_id'].null_count() or reference.filter(pl.col('cdr3').is_null() |
            ~pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$')).height:
        raise ValueError('graph reference identities and full junctions must be valid')
    gene = pl.col('v').cast(pl.String).str.strip_chars().str.to_uppercase()
    valid_v = gene.str.contains(r'^TR[AB]V[A-Z0-9/.-]+(?:\*[A-Z0-9]+)?$').fill_null(False)
    records = reference.with_columns(pl.when(valid_v).then(gene.str.split('*').list.first())
        .otherwise(pl.lit('__unknown_v__')).alias('_allele'))
    nodes = records.select('cdr3').unique().sort('cdr3').with_row_index('_node')
    if hamming_edges is None:
        pairs = overlap(nodes['cdr3'].to_list(), scope='1,0,0,1', threads=threads)
        edges = pairs.filter(pl.col('a_idx') < pl.col('b_idx')).select('a_idx', 'b_idx')
    else:
        if not {'cdr3_a', 'cdr3_b'}.issubset(hamming_edges.columns):
            raise ValueError('transient Hamming edges require cdr3_a/cdr3_b')
        edges = hamming_edges.join(nodes.rename({'cdr3': 'cdr3_a', '_node': 'a_idx'}), on='cdr3_a').join(
            nodes.rename({'cdr3': 'cdr3_b', '_node': 'b_idx'}), on='cdr3_b').select(
                pl.min_horizontal('a_idx', 'b_idx').cast(pl.Int64).alias('a_idx'),
                pl.max_horizontal('a_idx', 'b_idx').cast(pl.Int64).alias('b_idx')).unique().filter(pl.col('a_idx') < pl.col('b_idx'))
    full_rows, puncture_rows = [], []
    for key, part in records.partition_by('group_id', as_dict=True).items():
        junctions = part.select('cdr3').unique().sort('cdr3').with_row_index('_local').join(nodes, on='cdr3', validate='1:1').sort('_local')
        local_edges = edges.join(junctions.select(pl.col('_node').cast(pl.Int64).alias('a_idx'),
            pl.col('_local').alias('a')), on='a_idx').join(junctions.select(
                pl.col('_node').cast(pl.Int64).alias('b_idx'), pl.col('_local').alias('b')), on='b_idx').select('a', 'b').rows()
        graph, removed_graph = graph_profiles(junctions.height, local_edges, include_punctures=include_punctures)
        usage, removed_usage = composition(part, junctions, include_punctures=include_punctures)
        categories = part.select('_allele').unique().sort('_allele').with_row_index('_v')
        keys = part.select('cdr3', '_allele').unique().join(junctions.select('cdr3', '_local'), on='cdr3').join(categories, on='_allele')
        if junctions.height * categories.height * 8 * 8 > 256 * 1024**2:
            raise ValueError('colored graph subtree histograms exceed the 256 MiB allocation budget')
        colours = np.zeros((junctions.height, categories.height), dtype=np.int64)
        colours[keys['_local'].to_numpy(), keys['_v'].to_numpy()] = 1
        component_square, removed_component_square = _component_colour_squares(junctions.height, local_edges, colours, include_punctures=include_punctures)
        sizes = colours.sum(axis=1)
        lengths = junctions['cdr3'].str.len_chars().to_numpy()
        length_v = {length: colours[lengths == length].sum(axis=0) for length in np.unique(lengths)}
        length_counts = {length: int(counts.sum()) for length, counts in length_v.items()}
        length_square = sum(count * count for count in length_counts.values())
        length_v_square = sum(int(counts @ counts) for counts in length_v.values())
        size = keys.height
        context = dict(group_id=key[0], source_records=part.height, mixture_unit='unique full junction/V gene key; alleles pooled',
                       unknown_v_keys=keys.filter(pl.col('_allele') == '__unknown_v__').height)
        full_rows.append(context | graph | usage | dict(length_mass_square=length_square,
            length_v_mass_square=length_v_square, component_v_mass_square=component_square) |
            _weights(length_square, length_v_square, component_square, size))
        if not include_punctures:
            continue
        source_counts = dict(part.group_by('cdr3').len().iter_rows())
        unknown_index = categories.filter(pl.col('_allele') == '__unknown_v__')['_v'].to_list()
        for i, sequence in enumerate(junctions['cdr3']):
            length = lengths[i]
            remaining_length = length_counts[length] - sizes[i]
            a = length_square - length_counts[length] ** 2 + remaining_length ** 2
            after_v = length_v[length] - colours[i]
            b = length_v_square - int(length_v[length] @ length_v[length]) + int(after_v @ after_v)
            c = int(removed_component_square[i])
            puncture_rows.append(context | removed_graph[i] | removed_usage[i] | dict(
                excluded_junction=sequence, source_records=part.height-source_counts[sequence],
                unknown_v_keys=context['unknown_v_keys']-int(colours[i, unknown_index[0]])
                    if unknown_index else 0,
                length_mass_square=int(a), length_v_mass_square=int(b), component_v_mass_square=c) |
                _weights(int(a), int(b), c, int(size-sizes[i])))
    full = pl.DataFrame(full_rows, infer_schema_length=None).sort('group_id')
    punctures = (pl.DataFrame(puncture_rows, infer_schema_length=None).sort('group_id', 'excluded_junction')
                 if puncture_rows else full.head(0).with_columns(pl.lit(None, dtype=pl.String).alias('excluded_junction')))
    return full, punctures


def nearest_length_weights(reference: pl.DataFrame, graph_weights: pl.DataFrame, *,
                           species: str, locus: str, threads: int = 1) -> pl.DataFrame:
    """Candidate gap gate from fixed-reference cross-length nearest improvement.

    For each length, g=C/M: C resolved junction/V-allele keys have an unrestricted
    nearest distance strictly below their same-length nearest, among M keys with
    both neighbours. All distances use the gapblock channel's clipped/trim3/V-loop
    geometry. Ties give no cross-length preference. Absent comparisons use the
    pooled g; no pooled comparisons selects original alone. Null length denotes
    that pooled fallback. No query or outcome is used to estimate this table.

    Remaining mass preserves the old conditional original/stem4 ratio:
    w_original=(1-g)*o/(o+s), w_stem4=(1-g)*s/(o+s). These are descriptive gates,
    not probabilities of mechanisms. Geometry keys resolve alleles; graph keys
    pool V genes. Their declared populations are not silently equated.
    """
    from ..match.geometry import load_v_loops, resolve_v_alleles, reference_nearest_by_length

    if locus not in ('TRA', 'TRB') or threads < 1:
        raise ValueError('nearest-length weights require TRA/TRB and positive threads')
    required = {'group_id', *THREE_WEIGHTS}
    if (not required.issubset(graph_weights.columns) or graph_weights['group_id'].null_count() or
            graph_weights['group_id'].is_duplicated().any()):
        raise ValueError('one complete graph weight row required per reference group')
    weights = graph_weights.select(THREE_WEIGHTS).to_numpy()
    if (not np.isfinite(weights).all() or np.any(weights < 0) or
            not np.allclose(weights.sum(axis=1), 1., rtol=0., atol=1e-12) or
            np.any(weights[:, 0]+weights[:, 2] <= 0)):
        raise ValueError('graph weights require a finite convex sum and positive original/stem4 mass')
    groups = reference.select('group_id').unique().sort('group_id').with_row_index('_partition')
    if (groups.join(graph_weights.select('group_id'), on='group_id', how='anti').height or
            graph_weights.join(groups, on='group_id', how='anti').height):
        raise ValueError('graph weights and reference populations must match')
    model = load_v_loops(species)
    frame = reference.select('group_id', pl.col('cdr3').alias('seq0'),
        resolve_v_alleles(reference['v'], model, locus=locus).alias('v0')).filter(
            pl.col('v0').is_not_null() & pl.col('seq0').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]{8,}$'))
    frame = frame.unique().join(groups, on='group_id', validate='m:1')
    nearest, same = reference_nearest_by_length(frame, partition_ids=frame['_partition'].to_list(),
        species=species, threads=threads, ctrim=3, matrix='clipped', stem_depth=0)
    comparable = np.isfinite(nearest) & np.isfinite(same)
    measured = frame.with_columns(pl.col('seq0').str.len_chars().cast(pl.Int64).alias('junction_length'),
        pl.Series('_comparable', comparable), pl.Series('_closer', comparable & (nearest < same)))
    counts = [pl.len().alias('resolved_reference_keys'), pl.col('_comparable').sum().alias('comparison_keys'),
              pl.col('_closer').sum().alias('cross_length_closer_keys')]
    local = measured.group_by('group_id', 'junction_length').agg(*counts)
    pooled = measured.group_by('group_id').agg(*counts)
    scopes = reference.select('group_id', pl.col('cdr3').str.len_chars().cast(pl.Int64).alias('junction_length')).unique()
    scopes = pl.concat([scopes, groups.select('group_id', pl.lit(None, dtype=pl.Int64).alias('junction_length'))])
    values = pl.concat([local, pooled.with_columns(pl.lit(None, dtype=pl.Int64).alias('junction_length'))], how='diagonal_relaxed')
    result = scopes.join(values, on=['group_id', 'junction_length'], how='left', nulls_equal=True, validate='1:1').join(
        pooled.select('group_id', pl.col('comparison_keys').alias('pooled_comparison_keys'),
            pl.col('cross_length_closer_keys').alias('pooled_closer_keys')), on='group_id', how='left', validate='m:1')
    count_names = ['resolved_reference_keys', 'comparison_keys', 'cross_length_closer_keys',
                   'pooled_comparison_keys', 'pooled_closer_keys']
    result = result.with_columns(pl.col(count_names).fill_null(0)).join(
        graph_weights.select('group_id', *[pl.col(c).alias('legacy_'+c) for c in THREE_WEIGHTS]),
        on='group_id', validate='m:1')
    have_local = pl.col('junction_length').is_not_null() & (pl.col('comparison_keys') > 0)
    have_any = pl.col('pooled_comparison_keys') > 0
    result = result.with_columns(pl.when(have_local).then(pl.col('cross_length_closer_keys')/pl.col('comparison_keys'))
        .when(have_any).then(pl.col('pooled_closer_keys')/pl.col('pooled_comparison_keys')).otherwise(0.).alias('gapblock_weight'),
        pl.when(have_local).then(pl.lit('exact_length')).when(have_any).then(pl.lit('pooled'))
        .otherwise(pl.lit('original_fallback')).alias('weight_scope'))
    remainder = pl.col('legacy_original_weight') + pl.col('legacy_stem4_weight')
    return result.with_columns(
        pl.when(have_any).then((1-pl.col('gapblock_weight'))*pl.col('legacy_original_weight')/remainder).otherwise(1.).alias('original_weight'),
        pl.when(have_any).then((1-pl.col('gapblock_weight'))*pl.col('legacy_stem4_weight')/remainder).otherwise(0.).alias('stem4_weight')
    ).sort('group_id', 'junction_length')


def combine_reference_evidence(frame: pl.DataFrame, nulls: pl.DataFrame,
                               weights: pl.DataFrame) -> pl.DataFrame:
    """Apply model-matched graph weights and the existing normalized convex sum.

    This returns a score, not a P-value. Calibrate the complete sum/axis/label
    selection with fresh reference/query nulls and recomputed mock weights.
    """
    if 'group_id' not in frame.columns or 'group_id' not in weights.columns:
        raise ValueError('combined evidence and reference weights require group_id')
    if set(THREE_WEIGHTS).intersection(frame.columns):
        raise ValueError('combined input already contains reference weights')
    selected = weights.select('group_id', *THREE_WEIGHTS)
    if selected['group_id'].null_count() or selected['group_id'].n_unique() != selected.height:
        raise ValueError('reference weights require unique non-null group_id')
    if frame.select('group_id').unique().join(selected.select('group_id'), on='group_id', how='anti').height:
        raise ValueError('missing reference weights for combined evidence group')
    combined = frame.join(selected, on='group_id', how='left', validate='m:1', maintain_order='left')
    return null_normalized_evidence(combined, nulls, three_channel=True)
