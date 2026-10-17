"""Stream retained single-chain scores without allocating a query-by-label bank."""
from __future__ import annotations

import polars as pl

from ..aggregate.candidates import PMHC


def stream_components(queries, reference, control, *, species, locus, threads=1, consume):
    """Send ``(scorer, query_id/full-pMHC/score/supported frame)`` batches.

    Exact junctions are excluded from hits only; target counts and prior models
    remain fixed. ``current_original`` is historical ``unpaired_score``;
    ``current_original_log`` is its explicit log-scale variant. Under the same
    accepted control/model context, ``accepted_original`` and
    ``current_original`` are identical fresh evaluations for single chains; no
    stored predictions are loaded. Sparse targets already
    use log-density plus prior in both variants; dense targets use density in
    the raw variant. These are the retained formulas, not a new product score.

    Invalid or unavailable query/model cases may have no emitted row. The caller
    must retain its complete query axis and give such cases explicit dispositions.
    Original finite zero-density scores are retained with ``supported=False``.
    Paired cohort-rank fusion has no established all-label streaming identity;
    reject it explicitly rather than substitute a different marginal fusion.
    """
    from . import combined
    from ..evalue.control import _organism
    from ..match.geometry import precompute_v_loop_distances
    from ..match.historical import _prepare_prior_background, unpaired_score

    if threads < 1 or locus not in {'TRA', 'TRB', 'paired'}:
        raise ValueError('component scoring requires positive threads and TRA/TRB/paired locus')
    if locus == 'paired':
        raise NotImplementedError('streamed paired retained-component retrieval requires an explicit '
                                  'cohort-rank and linked/unlinked fusion contract')
    if not callable(consume):
        raise ValueError('streamed component scoring requires a batch consumer')
    if not isinstance(control, dict) or not {'density_index', 'prior_background'} <= control.keys():
        raise ValueError('accepted internal density_index/prior_background model context required')
    ids = queries['query_id']
    if (not ids.dtype.is_integer() or ids.null_count() or ids.n_unique() != queries.height or
            (queries.height and (ids.min() != 0 or ids.max() != queries.height - 1))):
        raise ValueError('component query identities must be contiguous integers starting at zero')
    if not {'cdr3', 'v', 'j', 'locus'} <= set(queries.columns):
        raise ValueError('single-chain component queries require cdr3/v/j/locus')
    # A caller must partition loci; silently rewriting another chain loses identity.
    if queries.filter(pl.col('locus').ne_missing(locus)).height:
        raise ValueError('component query locus differs from the requested population')
    organism = _organism(species)
    reference = combined._reference_population(reference, organism).filter(pl.col('gene') == locus)
    if not queries.height or not reference.height:
        return
    counts = reference.group_by(PMHC, maintain_order=True).agg(
        pl.col('cdr3').n_unique().alias('_fixed_reference_size'))
    prior_counts = control.get('_prior_counts')
    if prior_counts is None:
        prior_counts = _prepare_prior_background(control['prior_background'])
    contexts = control.get('v_loop_distances')
    if contexts is None:
        contexts = {matrix: precompute_v_loop_distances(species=organism, alphabet='full',
            matrix=matrix, v_regions='cdr12plus') for matrix in ('clipped', 'blosum62')}

    def original(frame):
        frame = frame.join(counts, on=PMHC, how='left', nulls_equal=True, validate='m:1')
        valid = frame.filter(pl.col('original_score').is_not_null())
        raw = valid.select('query_id', *PMHC).with_columns(unpaired_score(
            valid['original_density'], valid['germline_lr'], valid['_fixed_reference_size']).alias('score'))
        raw = frame.select('query_id', *PMHC, 'original_supported').join(raw,
            on=['query_id', *PMHC], how='left', nulls_equal=True, validate='1:1').with_columns(
            (pl.col('original_supported') & pl.col('score').is_not_null()).alias('supported')).drop('original_supported')
        consume('accepted_original', raw)
        consume('current_original', raw)
        consume('current_original_log', frame.select('query_id', *PMHC,
            pl.col('original_score').alias('score'),
            (pl.col('original_supported') & pl.col('original_score').is_not_null()).alias('supported')))

    combined._original(queries, reference, None, control['density_index'], control['prior_background'],
        locus=locus, threads=threads, exclude_exact=True, prior_counts=prior_counts, _consume=original)

    def kernel(name, frame):
        consume(name, frame.select('query_id', *PMHC,
            pl.when(pl.col('kernel_score') > 0).then(pl.col('kernel_score')).otherwise(None).alias('score'),
            (pl.col('kernel_score') > 0).fill_null(False).alias('supported')))

    combined._kernel(queries, reference, None, species=organism, threads=threads, paired=False,
        exclude_exact=True, contexts=contexts, _consume=kernel)

    from ..match.motifs import masked_motif_scores
    motif_queries=combined._valid_single(queries).filter(
        pl.col('v').str.strip_chars().ne('').fill_null(False))
    if motif_queries.height:
        masked_motif_scores(motif_queries,reference,exclude_exact=True,
            _consume=lambda frame: consume('motif',frame.select('query_id',*PMHC,
                pl.col('motif_score').alias('score'),(pl.col('supported_motifs')>0).alias('supported'))))
