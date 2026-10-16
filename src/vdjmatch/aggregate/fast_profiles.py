"""Literal fixed-reference classes from V concentration and edit-1 enrichment."""
from __future__ import annotations

import math
import polars as pl

from .structure import _gene_expression

PROFILE_RECIPE = 'v1:source>=30;V-key-majority>=1/2;nonidentical-Levenshtein1-within/rest-pair-density>1;pooled-and-pMHC'
_MHC = ['mhc_a', 'mhc_b', 'mhc_class']
_CONTEXT = ['level', 'species', 'gene', 'epitope', *_MHC]
_SCHEMA = {name: pl.String for name in ['level', 'species', 'locus', 'epitope', *_MHC,
    'mhc_context_status', 'support_status', 'neighbor_status', 'max_v_gene', 'profile_class']}
_SCHEMA.update({name: pl.Int64 for name in ['source_records', 'valid_source_records',
    'excluded_source_records', 'junctions', 'background_junctions', 'exact_shared_junctions',
    'within_opportunities', 'between_opportunities', 'within_neighbor_pairs', 'between_neighbor_pairs',
    'neighbor_junctions', 'reference_v_keys', 'unresolved_v_keys']})
_SCHEMA.update({name: pl.Float64 for name in ['within_pair_density', 'between_pair_density',
    'neighbor_enrichment', 'neighbor_coverage', 'max_v_gene_share', 'unresolved_v_fraction']})
_SCHEMA.update({name: pl.Boolean for name in ['classification_eligible', 'v_majority', 'neighbor_enriched']})


def _near_lengths(histogram, keys):
    """Sum counts at lengths L-1,L,L+1 using only small length histograms."""
    return pl.concat([histogram.with_columns((pl.col('_length')+delta).alias('_length'))
        for delta in (-1, 0, 1)]).group_by(*keys, '_length').agg(pl.col('_n').sum())


def _population(rows, *, locus, threads):
    """Reduce one species/locus graph; no per-target search or full-frame scans."""
    from .. import cluster

    selected = rows.filter(pl.col('_valid'))
    nodes = selected.select('cdr3').unique().sort('cdr3').with_row_index('_node').with_columns(
        pl.col('_node').cast(pl.Int64), pl.col('cdr3').str.len_chars().cast(pl.Int64).alias('_length'))
    membership = selected.select('_target', 'epitope', 'cdr3').unique().join(nodes, on='cdr3', validate='m:1')
    # The native wrapper owns hit flattening. One whole-population batch only.
    edges = (cluster.overlap(nodes['cdr3'].to_list(), scope='1,1,1,1', threads=threads)
        .select(pl.col('a_idx').alias('_a'), pl.col('b_idx').alias('_b')) if nodes.height else
        pl.DataFrame(schema={'_a': pl.Int64, '_b': pl.Int64}))
    degrees = edges.group_by('_a').len(name='_degree').with_columns(pl.col('_degree').cast(pl.Int64))
    peptide_nodes = membership.select('epitope', '_node', '_length').unique()
    exclusive = peptide_nodes.with_columns(pl.len().over('_node').alias('_labels')).filter(
        pl.col('_labels') == 1).drop('_labels')
    exclusive_neighbors = edges.join(exclusive.select(pl.col('_node').alias('_b'), 'epitope'),
        on='_b', validate='m:1').group_by('_a', 'epitope').len(name='_exclusive_neighbors')
    node_counts = membership.join(degrees, left_on='_node', right_on='_a', how='left', validate='m:1').join(
        exclusive_neighbors, left_on=['_node', 'epitope'], right_on=['_a', 'epitope'],
        how='left', validate='m:1').join(exclusive.select('_node', 'epitope').with_columns(
            pl.lit(1).alias('_exclusive')), on=['_node', 'epitope'], how='left', validate='m:1').with_columns(
                pl.col('_degree', '_exclusive_neighbors', '_exclusive').fill_null(0))
    summary = node_counts.group_by('_target').agg(pl.len().alias('junctions'),
        (pl.col('_degree')-pl.col('_exclusive_neighbors')).sum().alias('between_neighbor_pairs'),
        (pl.len()-pl.col('_exclusive').sum()).alias('exact_shared_junctions'))
    # Join only graph endpoints sharing the same target; never expand all label pairs.
    within = edges.join(membership.select('_target', pl.col('_node').alias('_a')), on='_a', validate='m:m').join(
        membership.select('_target', pl.col('_node').alias('_b')), on=['_target', '_b'], validate='m:1')
    within_counts = within.group_by('_target').agg((pl.len()//2).alias('within_neighbor_pairs'),
        pl.col('_a').n_unique().alias('neighbor_junctions'))
    histogram = membership.group_by('_target', 'epitope', '_length').len(name='_n').with_columns(pl.col('_n').cast(pl.Int64))
    global_histogram = nodes.group_by('_length').len(name='_n').with_columns(pl.col('_n').cast(pl.Int64))
    exclusive_histogram = exclusive.group_by('epitope', '_length').len(name='_n').with_columns(pl.col('_n').cast(pl.Int64))
    opportunities = histogram.join(_near_lengths(histogram.select('_target', '_length', '_n'), ['_target']).rename(
        {'_n': '_own_near'}), on=['_target', '_length'], validate='1:1').join(
        _near_lengths(global_histogram, []).rename({'_n': '_global_near'}), on='_length', validate='m:1').join(
        _near_lengths(exclusive_histogram, ['epitope']).rename({'_n': '_exclusive_near'}),
        on=['epitope', '_length'], how='left', validate='m:1').with_columns(pl.col('_exclusive_near').fill_null(0)).group_by(
            '_target').agg(((pl.col('_n')*(pl.col('_own_near')-1)).sum()//2).alias('within_opportunities'),
                (pl.col('_n')*(pl.col('_global_near')-pl.col('_exclusive_near'))).sum().alias('_cross_including_self'))
    summary = summary.join(within_counts, on='_target', how='left', validate='1:1').join(
        opportunities, on='_target', validate='1:1').with_columns(pl.col('within_neighbor_pairs', 'neighbor_junctions').fill_null(0)).with_columns(
            (pl.col('_cross_including_self')-pl.col('exact_shared_junctions')).alias('between_opportunities'))
    # Source multiplicity never weights either the graph or V-gene concentration.
    genes = selected.select('_target', 'cdr3', _gene_expression(locus, 'V').alias('_v')).unique().group_by(
        '_target', '_v').len(name='_keys').sort('_target', '_keys', '_v', descending=[False, True, False], nulls_last=True)
    gene_counts = genes.group_by('_target').agg(pl.col('_keys').sum().alias('reference_v_keys'),
        pl.col('_keys').filter(pl.col('_v').is_null()).sum().alias('unresolved_v_keys'),
        pl.col('_v').filter(pl.col('_v').is_not_null()).first().alias('max_v_gene'),
        pl.col('_keys').filter(pl.col('_v').is_not_null()).first().fill_null(0).alias('_max_v_keys'))
    excluded_counts = dict(exclusive.group_by('epitope').len().iter_rows())
    result = summary.join(gene_counts, on='_target', validate='1:1')
    by_target = {row['_target']: row for row in result.to_dicts()}
    return {target: by_target.get(target, {}) | {'background_junctions': nodes.height-excluded_counts.get(epitope, 0)}
        for target, epitope in rows.select('_target', 'epitope').unique().iter_rows()}


def reference_profiles(reference: pl.DataFrame, *, threads=4) -> pl.DataFrame:
    """Describe pooled epitopes and full pMHC groups in a normalized reference.

    Keep source-record eligibility separate from unique-junction graph evidence.
    Background is the union of junctions annotated to any other peptide within
    the same species/locus, including junctions also shared with the target.
    Modal pooled MHC is descriptive only. No pairing, binding or confidence claim
    follows from a pair of marginal chain classes.
    """
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise ValueError('threads must be a positive integer')
    required = {'species', 'gene', 'epitope', 'cdr3', 'v', 'reference_valid', *_MHC}
    if required-set(reference.columns):
        raise ValueError('fast profiles require the normalized reference schema')
    if reference.is_empty():
        return pl.DataFrame(schema=_SCHEMA)
    from ..cli.epitope import pool_epitopes
    raw = reference.select(*required).with_columns(
        (pl.col('reference_valid').fill_null(False) &
         pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)).alias('_valid'))
    pooled, _ = pool_epitopes(raw, min_records=1)
    labels = pooled.select('species', 'epitope', *_MHC).unique()
    all_rows = pl.concat([raw.drop(_MHC).join(labels, on=['species', 'epitope'], how='left',
        nulls_equal=True, validate='m:1').with_columns(pl.lit('epitope').alias('level')),
        raw.with_columns(pl.lit('pmhc').alias('level'))], how='diagonal_relaxed')
    targets = all_rows.select(_CONTEXT).unique().sort(_CONTEXT, nulls_last=True).with_row_index('_target')
    all_rows = all_rows.join(targets, on=_CONTEXT, nulls_equal=True, validate='m:1')
    inventory = all_rows.group_by('_target').agg(pl.len().alias('source_records'),
        pl.col('_valid').sum().alias('valid_source_records')).join(targets, on='_target', validate='1:1').sort('_target')
    profiles = {}
    for (_, locus), population in all_rows.partition_by(['species', 'gene'], as_dict=True, maintain_order=True).items():
        profiles.update(_population(population, locus=locus, threads=threads))
    output = []
    for target in inventory.to_dicts():
        data = profiles.get(target['_target'], {})
        n, bg = data.get('junctions', 0), data.get('background_junctions', 0)
        m, cross = data.get('within_neighbor_pairs', 0), data.get('between_neighbor_pairs', 0)
        own_op, cross_op = data.get('within_opportunities', 0), data.get('between_opportunities', 0)
        own_rate, cross_rate = m/own_op if own_op else None, cross/cross_op if cross_op else None
        if not n or not bg or not own_op or not cross_op:
            ratio, enriched, status = None, None, 'insufficient_pair_opportunities'
        elif not m and not cross:
            ratio, enriched, status = None, None, 'no_neighbor_evidence'
        else:
            ratio = own_rate/cross_rate if cross_rate else math.inf
            enriched, status = ratio > 1, 'zero_between_neighbors' if not cross else 'observed'
        keys, missing, maximum = data.get('reference_v_keys', 0), data.get('unresolved_v_keys', 0), data.get('_max_v_keys', 0)
        v_majority = (True if 2*maximum >= keys else None if 2*(maximum+missing) >= keys else False) if keys else None
        eligible = target['valid_source_records'] >= 30
        label = ('v_majority' if v_majority else 'no_v_majority')+'__'+('neighbor_enriched' if enriched else 'not_neighbor_enriched')
        if not eligible:
            label = 'insufficient_source_records'
        elif v_majority is None or enriched is None:
            label = 'unresolved_axes'
        output.append({key: target[key] for key in ['level', 'species', 'epitope', *_MHC]} |
            dict(locus=target['gene'], source_records=target['source_records'],
                valid_source_records=target['valid_source_records'],
                excluded_source_records=target['source_records']-target['valid_source_records'],
                classification_eligible=eligible, support_status='eligible' if eligible else 'insufficient_source_records',
                mhc_context_status='modal_annotation' if target['level']=='epitope' else
                    'incomplete_mhc' if any(not target[c] for c in _MHC) else 'observed_full_pmhc',
                junctions=n, background_junctions=bg, exact_shared_junctions=data.get('exact_shared_junctions', 0),
                within_opportunities=own_op, between_opportunities=cross_op,
                within_neighbor_pairs=m, between_neighbor_pairs=cross, neighbor_junctions=data.get('neighbor_junctions', 0),
                within_pair_density=own_rate, between_pair_density=cross_rate,
                neighbor_enrichment=ratio, neighbor_coverage=data.get('neighbor_junctions', 0)/n if n else None,
                neighbor_enriched=enriched, neighbor_status=status, reference_v_keys=keys, unresolved_v_keys=missing,
                unresolved_v_fraction=missing/keys if keys else None, max_v_gene=data.get('max_v_gene'),
                max_v_gene_share=maximum/keys if keys else None, v_majority=v_majority, profile_class=label))
    return pl.DataFrame(output, schema=_SCHEMA).sort('species', 'locus', 'level', 'epitope', *_MHC, nulls_last=True)
