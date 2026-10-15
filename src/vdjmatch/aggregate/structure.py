"""Reference graph summaries with exact leave-one-junction-out connectivity.

A single iterative depth-first traversal gives the components remaining after
removing each vertex, including articulation vertices. No per-query graph search.
"""
from __future__ import annotations
import math
import operator


def graph_profiles(n, edges, *, include_punctures=True):
    """Return full and vertex-punctured graph counts in O(vertices + edges).

    Vertices are unique full junctions; edges join one-substitution neighbours.
    Removing a vertex removes every associated V annotation from the caller's
    composition table. This function supplies connectivity, not recognition odds.
    Normal scoring sets include_punctures=False and returns no deletion diagnostics.
    """
    if isinstance(n, bool) or not isinstance(n, int) or n < 0:
        raise ValueError('nonnegative integer vertex count required')
    unique = set()
    for a, b in edges:
        try:
            a, b = operator.index(a), operator.index(b)
        except TypeError as exc:
            raise ValueError('integer vertex IDs required') from exc
        if not (0 <= a < n and 0 <= b < n) or a == b:
            raise ValueError('edges must join distinct valid vertex IDs')
        unique.add((min(a,b), max(a,b)))
    neighbours = [[] for _ in range(n)]
    for a,b in unique:
        neighbours[a].append(b); neighbours[b].append(a)
    seen = [-1]*n; low = [0]*n; parent = [-1]*n; sizes = [1]*n
    component = [-1]*n; separated = [[] for _ in range(n)]; component_sizes = []
    clock = 0
    for root in range(n):
        if seen[root] >= 0:
            continue
        cid = len(component_sizes)
        seen[root] = low[root] = clock; clock += 1; component[root] = cid
        stack = [(root, iter(neighbours[root]))]
        while stack:
            u, iterator = stack[-1]
            v = next(iterator, None)
            if v is not None:
                if v == parent[u]:
                    continue
                if seen[v] < 0:
                    parent[v] = u; component[v] = cid
                    seen[v] = low[v] = clock; clock += 1
                    stack.append((v, iter(neighbours[v])))
                else:
                    low[u] = min(low[u], seen[v])
            else:
                stack.pop()
                p = parent[u]
                if p >= 0:
                    sizes[p] += sizes[u]; low[p] = min(low[p], low[u])
                    if low[u] >= seen[p]:
                        separated[p].append(sizes[u])
        component_sizes.append(sizes[root])
    largest = [0,0]
    for size in component_sizes:
        if size >= largest[0]:
            largest = [size,largest[0]]
        else:
            largest[1] = max(largest[1],size)
    size_log = lambda s: s*math.log(s) if s else 0.
    total_log = math.fsum(size_log(s) for s in component_sizes)
    total_square = sum(s*s for s in component_sizes)
    singletons = sum(s == 1 for s in component_sizes)

    def summary(nodes, count, isolated, largest_size, edges_count, log_sum, square_sum):
        return dict(junctions=nodes, components=count, isolated_junctions=isolated,
            graph_edges=edges_count, largest_component=largest_size,
            largest_component_fraction=largest_size/nodes if nodes else None,
            clustered_fraction=1-isolated/nodes if nodes else None,
            component_entropy=max(0.,math.log(nodes)-log_sum/nodes) if nodes else None,
            size_weighted_component_mean=square_sum/nodes if nodes else None)

    full = summary(n,len(component_sizes),singletons,largest[0] if largest else 0,
                   len(unique),total_log,total_square)
    if not include_punctures:
        return full, []
    punctured = []
    for u in range(n):
        old = component_sizes[component[u]]
        fragments = separated[u]
        remainder = old-1-sum(fragments)
        if remainder:
            fragments = [*fragments, remainder]
        unaffected = largest[0] if old != largest[0] else largest[1] if len(largest)>1 else 0
        punctured.append(summary(n-1,len(component_sizes)-1+len(fragments),
            singletons-(old == 1)+sum(s == 1 for s in fragments),
            max(unaffected,max(fragments,default=0)),len(unique)-len(neighbours[u]),
            total_log-size_log(old)+math.fsum(size_log(s) for s in fragments),
            total_square-old*old+sum(s*s for s in fragments)))
    return full, punctured


def composition(records, junctions, *, include_punctures=True):
    """Full and exact-junction-punctured length/V summaries; no query labels."""
    import math
    import numpy as np
    import polars as pl

    seqs = junctions['cdr3'].to_list()
    lengths = np.array([len(s) for s in seqs], dtype=np.float64)
    keys = records.select('cdr3','_allele').unique().with_columns(
        pl.col('_allele').str.split('*').list.first().alias('_gene'))
    counts = keys.group_by('_gene').len().sort('_gene')
    sizes = counts['len'].to_numpy().astype(float)
    total = keys.height
    flog = lambda x: x*math.log(x) if x else 0.
    v_log = math.fsum(flog(x) for x in sizes)
    length_counts = dict(zip(*np.unique(lengths, return_counts=True)))
    length_log = math.fsum(flog(x) for x in length_counts.values())
    sum_length = lengths.sum(); sum_square = (lengths*lengths).sum()
    records_per_seq = dict(records.group_by('cdr3').len().iter_rows())
    removed_keys = dict(keys.group_by('cdr3').len().iter_rows())
    # Each sequence removes its full contribution from every V category.
    v_changes = {}
    if include_punctures:
        changes = keys.group_by('cdr3','_gene').len().rename({'len':'_removed'}).join(
            counts,on='_gene',validate='m:1').sort('cdr3','_gene')
        for sequence, part in changes.group_by('cdr3',maintain_order=True):
            before = part['len'].to_list(); after = (part['len']-part['_removed']).to_list()
            v_changes[sequence[0]] = (math.fsum(flog(new)-flog(old) for old,new in zip(before,after)),
                                      sum(x==0 for x in after))

    def row(removed=None):
        n = len(seqs)-(removed is not None)
        length = len(removed) if removed is not None else 0
        k = removed_keys[removed] if removed is not None else 0
        remain = total-k
        delta, lost_genes = v_changes[removed] if removed is not None else (0.,0)
        v_entropy = max(0., math.log(remain)-(v_log+delta)/remain) if remain else None
        f_length = length_log
        if removed is not None:
            count = length_counts[float(length)]
            f_length += flog(count-1)-flog(count)
        mean = (sum_length-length)/n if n else None
        variance = max(0.,(sum_square-length*length)/n-mean*mean) if n else None
        return dict(clonotype_keys=remain, v_genes=len(sizes)-lost_genes,
            v_entropy=v_entropy, effective_v_genes=math.exp(v_entropy) if v_entropy is not None else None,
            length_mean=mean, length_sd=math.sqrt(variance) if variance is not None else None,
            length_entropy=max(0.,math.log(n)-f_length/n) if n else None,
            reference_records=records.height-(records_per_seq[removed] if removed is not None else 0),
            excluded_reference_records=records_per_seq[removed] if removed is not None else 0,
            excluded_clonotype_keys=k)
    return row(), [row(sequence) for sequence in seqs] if include_punctures else []


def _gene_expression(locus, segment):
    """One shared allele-pooled gene projection; unresolved calls stay null."""
    import polars as pl
    gene = pl.col(segment.lower()).cast(pl.String).str.strip_chars().str.to_uppercase()
    resolved = gene.str.contains('^' + locus + segment + r'[A-Z0-9/.-]+(?:\*[A-Z0-9]+)?$').fill_null(False)
    return pl.when(resolved).then(gene.str.split('*').list.first()).otherwise(None)


def reference_gene_usage(records, *, locus):
    """Count fixed-reference V/J usage overall and within each junction length.

    Alleles are pooled into unambiguous genes; unresolved calls remain a null
    category. Record fractions use all supplied records. Key fractions use
    distinct full-junction/gene keys within the segment and length scope,
    including unresolved keys. These are composition, not background enrichment.
    """
    import polars as pl

    if locus not in {'TRA', 'TRB'} or not {'cdr3', 'v', 'j'}.issubset(records.columns):
        raise ValueError('gene usage requires TRA/TRB and cdr3/v/j reference columns')
    if records['cdr3'].null_count() or not records['cdr3'].str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').all():
        raise ValueError('gene usage requires valid full amino-acid junctions')
    parts = []
    for segment in ('V', 'J'):
        frame = records.select('cdr3', pl.col('cdr3').str.len_chars().alias('junction_length'),
            _gene_expression(locus, segment).alias('gene'))
        for scope in ('all_lengths', 'exact_length'):
            selected = frame.with_columns(pl.lit(None, dtype=pl.UInt32).alias('junction_length')) if scope == 'all_lengths' else frame
            counts = selected.group_by('junction_length', 'gene').agg(
                pl.len().alias('reference_records'), pl.col('cdr3').n_unique().alias('junction_gene_keys'))
            parts.append(counts.with_columns(pl.lit(segment).alias('segment'),
                pl.lit(scope).alias('length_scope'),
                pl.when(pl.col('gene').is_not_null()).then(pl.lit('resolved')).otherwise(pl.lit('unresolved')).alias('gene_status'),
                pl.col('reference_records').sum().over('junction_length').alias('total_reference_records'),
                pl.col('junction_gene_keys').sum().over('junction_length').alias('total_junction_gene_keys')).with_columns(
                    (pl.col('reference_records') / pl.col('total_reference_records')).alias('record_fraction'),
                    (pl.col('junction_gene_keys') / pl.col('total_junction_gene_keys')).alias('key_fraction')))
    return pl.concat(parts).sort('segment', 'length_scope', 'junction_length', 'gene', nulls_last=True)


def control_gene_usage(events, *, locus):
    """V/J composition of normalized control_events, pooling gene alleles.

    Event fractions preserve recombination-event multiplicity; key fractions use
    unique full-junction/gene keys. Read abundance never supplies event weights.
    The caller selects species/locus through control_events before this reduction.
    """
    import polars as pl
    from ..io.columns import _integer_count

    if not {'junction_aa', 'v', 'j', 'event_weight'}.issubset(events.columns) or not events.height:
        raise ValueError('control gene usage requires nonempty normalized control events')
    events = events.with_columns(_integer_count('event_weight').alias('event_weight'))
    if events['event_weight'].null_count() or (events['event_weight'] <= 0).any():
        raise ValueError('control event weights must be positive integers')
    usage = reference_gene_usage(events.rename({'junction_aa': 'cdr3'}), locus=locus)
    if usage['gene'].null_count():
        raise ValueError('control events require unambiguous matching V/J genes')
    parts = []
    for segment in ('V', 'J'):
        frame = events.select(pl.col('junction_aa').str.len_chars().alias('junction_length'),
            pl.col(segment.lower()).cast(pl.String).str.strip_chars().str.to_uppercase()
                .str.split('*').list.first().alias('gene'), 'event_weight')
        for scope in ('all_lengths', 'exact_length'):
            selected = frame.with_columns(pl.lit(None, dtype=pl.UInt32).alias('junction_length')) if scope == 'all_lengths' else frame
            parts.append(selected.group_by('junction_length', 'gene').agg(
                pl.col('event_weight').sum().alias('control_events')).with_columns(
                pl.lit(segment).alias('segment'), pl.lit(scope).alias('length_scope'),
                pl.col('control_events').sum().over('junction_length').alias('total_control_events')))
    keys = ['segment', 'length_scope', 'junction_length', 'gene']
    return usage.select(*keys, pl.col('junction_gene_keys').alias('control_junction_gene_keys'),
        pl.col('total_junction_gene_keys').alias('total_control_junction_gene_keys')).join(
        pl.concat(parts), on=keys, nulls_equal=True, validate='1:1', maintain_order='left').with_columns(
        (pl.col('control_events') / pl.col('total_control_events')).alias('control_event_fraction'),
        (pl.col('control_junction_gene_keys') / pl.col('total_control_junction_gene_keys')).alias('control_key_fraction'))


def compare_gene_usage(usage, control):
    """Descriptive reference/control ratios with unchanged reference denominators.

    One species/locus control table supplies all reference pMHC populations.
    Missing reference genes, unobserved control genes and missing control lengths
    remain distinct. No pseudocount, tail, selection claim or binding confidence.
    """
    import polars as pl

    keys = ['segment', 'length_scope', 'junction_length', 'gene']
    strata = keys[:-1]
    totals = ['total_control_events', 'total_control_junction_gene_keys']
    counts = ['control_events', 'control_junction_gene_keys']
    if (not {*keys, 'record_fraction', 'key_fraction'}.issubset(usage.columns) or
            not {*keys, *totals, *counts}.issubset(control.columns) or control.select(keys).is_duplicated().any()):
        raise ValueError('gene usage comparison requires unique normalized control gene strata')
    result = usage.join(control.select(*strata, *totals).unique(), on=strata,
        how='left', nulls_equal=True, validate='m:1', maintain_order='left').join(
        control.select(*keys, *counts), on=keys, how='left', nulls_equal=True, validate='m:1', maintain_order='left')
    available = pl.col('gene').is_not_null() & pl.col('total_control_events').is_not_null()
    result = result.with_columns(*[pl.when(available).then(pl.col(c).fill_null(0)).otherwise(None).alias(c) for c in counts])
    result = result.with_columns(
        (pl.col('control_events') / pl.col('total_control_events')).alias('control_event_fraction'),
        (pl.col('control_junction_gene_keys') / pl.col('total_control_junction_gene_keys')).alias('control_key_fraction'),
        pl.when(pl.col('gene').is_null()).then(pl.lit('unresolved_reference_gene'))
          .when(pl.col('total_control_events').is_null()).then(pl.lit('no_control_length'))
          .when(pl.col('control_events') == 0).then(pl.lit('gene_not_observed'))
          .otherwise(pl.lit('available')).alias('control_status'))
    return result.with_columns(*[pl.when(pl.col('control_event_fraction') > 0).then(
        pl.col(source) / pl.col('control_event_fraction')).otherwise(None).alias(name)
        for source, name in [('record_fraction', 'record_to_control_event_ratio'),
                             ('key_fraction', 'key_to_control_event_ratio')]])


def gene_selection_information(comparison):
    """Pointwise selection bits and conditional V/J--length information.

    Consume :func:`compare_gene_usage` output with species, locus and full pMHC
    columns. Both populations use distinct full-junction/segment-gene keys;
    V and J are separate projections, not independent evidence to add.
    Return the unchanged input rows/columns plus pointwise log2 ratios, and one
    summary per context/segment/length scope. ``all_lengths`` measures marginal
    gene information; ``exact_length`` decomposes joint information into length
    and target-weighted conditional gene information. Null restrictions remain
    valid grouping keys. Zero/missing control support or unresolved reference
    genes retain their original mass and make affected information unavailable;
    no pseudocount, renormalization, score selection or confidence is introduced.
    """
    import polars as pl
    from .candidates import PMHC

    context = ['species', 'locus', *PMHC, 'segment']
    keys = [*context, 'length_scope', 'junction_length', 'gene']
    required = {*keys, 'junction_gene_keys', 'total_junction_gene_keys',
                'total_control_junction_gene_keys', 'reference_records', 'key_fraction', 'control_key_fraction'}
    if required - set(comparison.columns):
        raise ValueError('selection information requires compared V/J gene counts and full pMHC context')
    if comparison.select(keys).is_duplicated().any():
        raise ValueError('selection information requires unique context/gene/length rows')
    if (comparison.filter(~pl.col('length_scope').is_in(['all_lengths', 'exact_length']).fill_null(False)).height or
            comparison.filter((pl.col('junction_gene_keys') <= 0) |
                              (pl.col('total_junction_gene_keys') < pl.col('junction_gene_keys'))).height):
        raise ValueError('selection information requires fixed positive gene counts and declared length scopes')
    frame = comparison.with_row_index('_selection_row')
    totals = frame.filter(pl.col('length_scope') == 'all_lengths').select(*context,
        pl.col('total_junction_gene_keys').alias('reference_total'),
        pl.col('total_control_junction_gene_keys').alias('control_total')).unique()
    rows = frame.join(totals, on=context, how='left', nulls_equal=True, validate='m:1', maintain_order='left')
    if rows['reference_total'].null_count():
        raise ValueError('selection information requires an all-length marginal for every context')
    rows = rows.with_columns(
        (pl.col('total_junction_gene_keys') / pl.col('reference_total')).alias('reference_length_fraction'),
        (pl.col('total_control_junction_gene_keys') / pl.col('control_total')).alias('control_length_fraction'))
    rows = rows.with_columns(
        (pl.col('reference_length_fraction') / pl.col('control_length_fraction')).log(2).alias('length_bits'),
        pl.when((pl.col('control_key_fraction') > 0) & pl.col('gene').is_not_null()).then(
            (pl.col('key_fraction') / pl.col('control_key_fraction')).log(2)).otherwise(None).alias('gene_bits'),
        (pl.col('key_fraction') * pl.col('reference_length_fraction')).alias('joint_reference_fraction'))
    rows = rows.with_columns((pl.col('length_bits') + pl.col('gene_bits')).alias('joint_bits'))
    summary = rows.group_by(*context, 'length_scope').agg(
        pl.col('reference_records').sum(), pl.col('junction_gene_keys').sum().alias('reference_keys'),
        pl.col('control_total').first().alias('control_keys'),
        pl.col('junction_gene_keys').filter(pl.col('gene_bits').is_null()).sum().alias('unresolved_comparison_keys'),
        pl.col('junction_gene_keys').filter(pl.col('length_bits').is_null()).sum().alias('unavailable_length_keys'),
        (pl.col('joint_reference_fraction') * pl.col('length_bits')).sum().alias('length_information_bits'),
        (pl.col('joint_reference_fraction') * pl.col('gene_bits')).sum().alias('gene_information_bits'),
        (pl.col('joint_reference_fraction') * pl.col('joint_bits')).sum().alias('joint_information_bits'))
    complete_gene = pl.col('unresolved_comparison_keys') == 0
    complete_length = pl.col('unavailable_length_keys') == 0
    summary = summary.with_columns(
        pl.when(complete_length).then(pl.col('length_information_bits')).otherwise(None).alias('length_information_bits'),
        pl.when(complete_gene).then(pl.col('gene_information_bits')).otherwise(None).alias('gene_information_bits'),
        pl.when(complete_gene & complete_length).then(pl.col('joint_information_bits')).otherwise(None).alias('joint_information_bits'),
        pl.when(complete_gene & complete_length).then(pl.lit('available')).otherwise(pl.lit('incomplete_support')).alias('information_status'))
    return rows.sort('_selection_row').drop('_selection_row'), summary.sort(*context, 'length_scope', nulls_last=True)
