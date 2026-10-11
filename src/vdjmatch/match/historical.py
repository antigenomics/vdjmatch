"""Reproduce the manuscript's single-chain PSSM density component.

This is a baseline, not the full manuscript scorer: sparse-reference germline priors
and cohort rank fusion are separate calculations. Reference representatives follow
source order within each pMHC/junction, as in the original per-epitope producer.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from itertools import accumulate
import math

import polars as pl
from seqtree import Index, SearchParams
from seqtree.evalue import evalue_result

from ..aggregate.candidates import PMHC
from ..evalue import first_hit
from .regions import significance_pssm, gene_family
from .vgene import vsim


def control_histograms(control, sequences, *, threads=1, radii=None):
    """Transient cumulative unit-edit counts; no neighbour objects cross into Python."""
    if not hasattr(control, 'edit_histogram_batch'):
        raise RuntimeError('historical scoring requires seqtree edit_histogram_batch; install the repaired source')
    sequences=list(sequences)
    if radii is None:
        radii={s:5 for s in sequences}
    groups=defaultdict(list)
    for s in sequences:
        radius=radii[s]
        if isinstance(radius,bool) or not isinstance(radius,int) or not 0<=radius<=5:
            raise ValueError('control radius must be an integer from zero to five')
        if radius:groups[radius].append(s)
    out={s:[0]*6 for s in sequences}
    for radius,batch in groups.items():
        counts=control.edit_histogram_batch(batch,first_hit.scope(radius),threads,True)
        for s,row in zip(batch,counts):
            cumulative=list(accumulate(row))
            out[s]=cumulative+[cumulative[-1]]*(5-radius)
    return out


def _density_neighbours(seqs, rs, threads):
    idx = Index.build(rs, 'aa')
    by_length = defaultdict(list)
    for s in seqs:
        by_length[len(s)].append(s)
    hits = {}
    for length, batch in by_length.items():
        p = SearchParams(max_subs=5, max_ins=0, max_dels=0, max_total_edits=5, engine='seqtm')
        p.pos_matrix = significance_pssm(length)
        hits.update(zip(batch, idx.search_batch(batch, p, threads)))
    unit = idx.edit_histogram_batch(seqs, first_hit.scope(1,1,1), threads, True)
    radii = {s:max([1 if hs[1] else 0]+[h.n_subs for h in hits[s] if rs[h.ref_id]!=s])
             for s,hs in zip(seqs,unit)}
    return hits, unit, radii


def shared_control_histograms(control, batches, *, threads=1):
    """One transient control CDF per junction at its largest needed target radius.

    Reference searches are repeated during scoring so no cross-target neighbour
    graph is retained. Counts below the largest radius are exactly unchanged.
    """
    if threads < 1 or not len(control):
        raise ValueError('positive threads and a nonempty control are required')
    radii = {}
    for queries, reference in batches:
        if not queries.height or not reference.height:
            continue
        seqs = queries['cdr3'].unique(maintain_order=True).to_list()
        _, _, needed = _density_neighbours(seqs, reference.unique('cdr3', keep='first', maintain_order=True)['cdr3'].to_list(), threads)
        for s, radius in needed.items():
            radii[s] = max(radii.get(s, 0), radius)
    return control_histograms(control, radii, threads=threads, radii=radii)


def density(queries, reference, control, *, threads=1, pool_reference=False, control_counts=None,
            puncture_populations=False):
    """Return keyed PSSM density and original radius<=1 significance component.

    Query IDs and contradictory labels remain distinct. Exact junctions are excluded
    symmetrically. Neither significance nor density is posterior specificity.
    ``pool_reference`` reproduces the original task-level reference union; the caller
    selects compatible MHC records first. Heterogeneous MHC output fields are null.
    Candidate ``puncture_populations`` also removes exact matches from population
    sizes. An empty remaining population yields missing score/significance.
    """
    if threads < 1 or not len(control):
        raise ValueError("positive threads and a nonempty control are required")
    seqs = queries['cdr3'].unique(maintain_order=True).to_list()
    control_sizes = dict.fromkeys(seqs, len(control))
    if puncture_populations:
        exact = Counter(control.ref_seqs())
        control_sizes = {s:len(control)-exact[s] for s in seqs}
    rows = []
    if pool_reference:
        for column in ('epitope', 'species', 'gene'):
            if column in reference.columns and reference[column].n_unique() > 1:
                raise ValueError('pooled reference requires one '+column)
        key = tuple(reference[c][0] if reference.height and reference[c].n_unique() == 1
                    else None for c in PMHC)
        groups = [(key, reference)] if reference.height else []
    else:
        groups = reference.group_by(PMHC, maintain_order=True)
        if control_counts is None and reference.select(PMHC).unique().height > 1:
            control_counts = shared_control_histograms(control,
                ((queries, obs) for _, obs in groups), threads=threads)
    for key, obs in groups:
        # Reproduce the original per-epitope reference representative rule explicitly.
        ref = obs.unique('cdr3', keep='first', maintain_order=True)
        rs, rv = ref['cdr3'].to_list(), ref['v'].to_list()
        hits, unit, radii = _density_neighbours(seqs, rs, threads)
        # Original significance used the nearest nonexact target radius, capped at one.
        # Only target-supported radii enter the original formula. Narrower control
        # balls have the identical CDF at those radii; unsupported queries need none.
        counts = control_histograms(control,seqs,threads=threads,radii=radii) if control_counts is None else control_counts
        reference_keys = set(rs)
        reference_sizes = {s:len(rs)-int(puncture_populations and s in reference_keys) for s in seqs}
        significance = {s:(evalue_result(hs[1],counts[s][1],reference_sizes[s],control_sizes[s])['p_enrichment']
                        if hs[1] else 1.0) if reference_sizes[s] and control_sizes[s] else None
                        for s,hs in zip(seqs,unit)}
        for query in queries.iter_rows(named=True):
            s, v = query['cdr3'], query['v']
            nr, nc_population = reference_sizes[s], control_sizes[s]
            total = 0.0 if nr and nc_population else None
            for hit in hits[s] if total is not None else []:
                if rs[hit.ref_id] == s:
                    continue
                ref_v = rv[hit.ref_id]
                weight = 1.0 if gene_family(v) == gene_family(ref_v) else .25*vsim(v,ref_v)
                nc = counts[s][hit.n_subs]
                total += weight*math.exp(-hit.score/400.0)/max(nr/nc_population*nc,.01)
            rows.append({'query_id':query['query_id'], **dict(zip(PMHC,key)),
                         'score':total, 'p_enrichment':significance[s],
                         'n_reference':nr,'control_size':nc_population,
                         'estimator':'punctured-pssm-density-v1' if puncture_populations else 'historical-pssm-density-v1'})
    schema={'query_id':queries.schema['query_id'], **{c:pl.String for c in PMHC},
            'score':pl.Float64,'p_enrichment':pl.Float64,'n_reference':pl.UInt64,
            'control_size':pl.UInt64,'estimator':pl.String}
    return pl.DataFrame(rows,schema=schema)


def germline_prior(queries, reference, background, *, include_length=True, exclude_exact=False):
    """Original smoothed V/J[/junction-length] log likelihood ratio, keyed by query.

    Frames use ``cdr3,v,j``. Reference junctions use their first source-order V/J;
    background rows retain raw multiplicity. Alpha reproduction deliberately passes
    the original beta background with ``include_length=False``; no correction is
    inferred here. Unknown/missing-input disposition belongs to the caller.
    ``exclude_exact`` removes every full-junction match from both populations and
    feature counts per query. This candidate option preserves historical defaults;
    an emptied population yields missing evidence, not a prior-only assignment.
    """
    ref = reference.unique('cdr3', keep='first', maintain_order=True)
    nr, nb = ref.height, background.height
    if not nr or not nb:
        raise ValueError('germline prior requires nonempty reference and background')
    columns = ['v', 'j'] + (['length'] if include_length else [])

    def features(frame):
        return frame.with_columns(
            pl.col('v', 'j').str.replace(r'\*.*$', ''),
            pl.col('cdr3').str.len_chars().alias('length'))

    out = features(queries).select('query_id', *(['cdr3'] if exclude_exact else []), *columns)
    ref, bg = features(ref), features(background)
    out = out.with_columns(pl.lit(0.0).alias('germline_lr'),
        pl.lit(nr, dtype=pl.UInt64).alias('n_reference'),
        pl.lit(nb, dtype=pl.UInt64).alias('background_size'))
    if exclude_exact:
        for population, size in [(ref, 'n_reference'), (bg, 'background_size')]:
            removed = population.group_by('cdr3').len().rename({'len': '_removed'})
            out = (out.join(removed, on='cdr3', how='left', validate='m:1', maintain_order='left')
                .with_columns((pl.col(size) - pl.col('_removed').fill_null(0)).alias(size))
                .drop('_removed'))
    for column in columns:
        rc = ref.group_by(column).len().rename({'len': '_ref_count'})
        bc = bg.group_by(column).len().rename({'len': '_bg_count'})
        out = (out.join(rc, on=column, how='left', maintain_order='left')
               .join(bc, on=column, how='left', maintain_order='left')
               .with_columns(pl.col('_ref_count', '_bg_count').fill_null(0)))
        if exclude_exact:
            keys = ['cdr3', column]
            for population, count in [(ref, '_ref_count'), (bg, '_bg_count')]:
                removed = population.group_by(keys).len().rename({'len': '_removed'})
                out = (out.join(removed, on=keys, how='left', validate='m:1', maintain_order='left')
                    .with_columns((pl.col(count) - pl.col('_removed').fill_null(0)).alias(count))
                    .drop('_removed'))
        out = (out
               .with_columns((pl.col('germline_lr') +
                   (((pl.col('_ref_count') + .5) / (pl.col('n_reference') + 30)) /
                    ((pl.col('_bg_count') + .5) / (pl.col('background_size') + 30))).log())
                   .alias('germline_lr'))
               .drop('_ref_count', '_bg_count'))
    if exclude_exact:
        out = out.with_columns(pl.when((pl.col('n_reference') > 0) & (pl.col('background_size') > 0))
            .then(pl.col('germline_lr')).otherwise(None).alias('germline_lr'))
    return out.select('query_id', 'germline_lr', 'n_reference', 'background_size')


def _score_series(values):
    values = pl.Series(values, dtype=pl.Float64)
    if values.null_count() or not values.is_finite().all():
        raise ValueError('original fusion requires finite, nonmissing component scores')
    return values


def unpaired_score(densities, priors, reference_size):
    """Original sparse (<500 unique reference junctions) or dense unpaired score."""
    density_values, prior_values = _score_series(densities), _score_series(priors)
    if isinstance(reference_size, pl.Series):
        sizes = _score_series(reference_size)
        if len(sizes) != len(density_values) or (sizes < 1).any() or (sizes != sizes.floor()).any():
            raise ValueError('aligned components and a positive integer reference size are required')
        if len(density_values) != len(prior_values) or (density_values < 0).any():
            raise ValueError('aligned components and nonnegative density are required')
        return pl.DataFrame({'density': density_values, 'prior': prior_values, 'size': sizes}).select(
            pl.when(pl.col('size') < 500).then((pl.col('density') + 1e-6).log() + pl.col('prior'))
            .otherwise(pl.col('density')).alias('score'))['score']
    if len(density_values) != len(prior_values) or reference_size < 1:
        raise ValueError('aligned components and a positive reference size are required')
    if (density_values < 0).any():
        raise ValueError('density must be nonnegative')
    score = ((density_values + 1e-6).log() + prior_values
             if reference_size < 500 else density_values)
    return score.rename('score')


def paired_score(beta_densities, alpha_densities, combined_priors):
    """Original equal-weight cohort rank sum; ties receive mean zero-based rank.

    Inputs are aligned observations, retaining repeated junctions and pair context.
    ``combined_priors`` is beta V/J/length plus alpha V/J, not a density product.
    """
    beta, alpha, prior = map(_score_series, (beta_densities, alpha_densities, combined_priors))
    if len(beta) != len(alpha) or len(beta) != len(prior):
        raise ValueError('paired components must describe the same aligned observations')
    if (beta < 0).any() or (alpha < 0).any():
        raise ValueError('density must be nonnegative')
    return (beta.rank('average') + alpha.rank('average') + prior.rank('average') - 3).rename('score')
