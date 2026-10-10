"""Reproduce the manuscript's single-chain PSSM density component.

This is a baseline, not the full manuscript scorer: sparse-reference germline priors
and cohort rank fusion are separate calculations. Reference representatives follow
source order within each pMHC/junction, as in the original per-epitope producer.
"""
from __future__ import annotations

from collections import defaultdict
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


def density(queries, reference, control, *, threads=1, pool_reference=False, control_counts=None):
    """Return keyed PSSM density and original radius<=1 significance component.

    Query IDs and contradictory labels remain distinct. Exact junctions are excluded
    symmetrically. Neither significance nor density is posterior specificity.
    ``pool_reference`` reproduces the original task-level reference union; the caller
    selects compatible MHC records first. Heterogeneous MHC output fields are null.
    """
    if threads < 1 or not len(control):
        raise ValueError("positive threads and a nonempty control are required")
    seqs = queries['cdr3'].unique(maintain_order=True).to_list()
    by_length = defaultdict(list)
    for s in seqs:
        by_length[len(s)].append(s)
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
    for key, obs in groups:
        # Reproduce the original per-epitope reference representative rule explicitly.
        ref = obs.unique('cdr3', keep='first', maintain_order=True)
        rs, rv = ref['cdr3'].to_list(), ref['v'].to_list()
        idx = Index.build(rs, 'aa')
        hits = {}
        for length, batch in by_length.items():
            p = SearchParams(max_subs=5, max_ins=0, max_dels=0, max_total_edits=5,
                             engine='seqtm')
            p.pos_matrix = significance_pssm(length)
            hits.update(zip(batch, idx.search_batch(batch, p, threads)))
        # Original significance used the nearest nonexact target radius, capped at one.
        unit = idx.edit_histogram_batch(seqs,first_hit.scope(1,1,1),threads,True)
        # Only target-supported radii enter the original formula. Narrower control
        # balls have the identical CDF at those radii; unsupported queries need none.
        radii={s:max([1 if hs[1] else 0]+[h.n_subs for h in hits[s] if rs[h.ref_id]!=s])
               for s,hs in zip(seqs,unit)}
        counts = control_histograms(control,seqs,threads=threads,radii=radii) if control_counts is None else control_counts
        significance = {s: evalue_result(hs[1],counts[s][1],len(rs),len(control))['p_enrichment']
                        if hs[1] else 1.0 for s,hs in zip(seqs,unit)}
        for query in queries.iter_rows(named=True):
            s, v = query['cdr3'], query['v']
            total = 0.0
            for hit in hits[s]:
                if rs[hit.ref_id] == s:
                    continue
                ref_v = rv[hit.ref_id]
                weight = 1.0 if gene_family(v) == gene_family(ref_v) else .25*vsim(v,ref_v)
                nc = counts[s][hit.n_subs]
                total += weight*math.exp(-hit.score/400.0)/max(len(rs)/len(control)*nc,.01)
            rows.append({'query_id':query['query_id'], **dict(zip(PMHC,key)),
                         'score':total, 'p_enrichment':significance[s],
                         'n_reference':len(rs),'control_size':len(control),
                         'estimator':'historical-pssm-density-v1'})
    schema={'query_id':queries.schema['query_id'], **{c:pl.String for c in PMHC},
            'score':pl.Float64,'p_enrichment':pl.Float64,'n_reference':pl.UInt64,
            'control_size':pl.UInt64,'estimator':pl.String}
    return pl.DataFrame(rows,schema=schema)


def germline_prior(queries, reference, background, *, include_length=True):
    """Original smoothed V/J[/junction-length] log likelihood ratio, keyed by query.

    Frames use ``cdr3,v,j``. Reference junctions use their first source-order V/J;
    background rows retain raw multiplicity. Alpha reproduction deliberately passes
    the original beta background with ``include_length=False``; no correction is
    inferred here. Unknown/missing-input disposition belongs to the caller.
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

    out = features(queries).select('query_id', *columns)
    ref, bg = features(ref), features(background)
    out = out.with_columns(pl.lit(0.0).alias('germline_lr'))
    for column in columns:
        rc = ref.group_by(column).len().rename({'len': '_ref_count'})
        bc = bg.group_by(column).len().rename({'len': '_bg_count'})
        out = (out.join(rc, on=column, how='left', maintain_order='left')
               .join(bc, on=column, how='left', maintain_order='left')
               .with_columns((pl.col('germline_lr') +
                   (((pl.col('_ref_count').fill_null(0) + .5) / (nr + 30)) /
                    ((pl.col('_bg_count').fill_null(0) + .5) / (nb + 30))).log())
                   .alias('germline_lr'))
               .drop('_ref_count', '_bg_count'))
    return out.select('query_id', 'germline_lr').with_columns(
        pl.lit(nr, dtype=pl.UInt64).alias('n_reference'),
        pl.lit(nb, dtype=pl.UInt64).alias('background_size'))


def _score_series(values):
    values = pl.Series(values, dtype=pl.Float64)
    if values.null_count() or not values.is_finite().all():
        raise ValueError('original fusion requires finite, nonmissing component scores')
    return values


def unpaired_score(densities, priors, reference_size):
    """Original sparse (<500 unique reference junctions) or dense unpaired score."""
    density_values, prior_values = _score_series(densities), _score_series(priors)
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
