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
from .vgene import precompute_v_similarities


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
            puncture_populations=False, v_similarities=None):
    """Return keyed PSSM density and original radius<=1 significance component.

    The original exp(-P/400) kernel is a tempered, self-normalised BLOSUM
    likelihood ratio: P=100*(A(q,q)+A(r,r)-2*A(q,r)), where A is the
    weighted raw similarity sum. Keep its aggregation and count P-value distinct.

    Query IDs and contradictory labels remain distinct. Exact junctions are excluded
    symmetrically. Neither significance nor density is posterior specificity.
    ``pool_reference`` reproduces the original task-level reference union; the caller
    selects compatible MHC records first. Heterogeneous MHC output fields are null.
    Candidate ``puncture_populations`` also removes exact matches from population
    sizes. An empty remaining population yields missing score/significance.
    """
    if threads < 1 or not len(control):
        raise ValueError("positive threads and a nonempty control are required")
    if v_similarities is None:
        v_similarities=precompute_v_similarities(queries['v'],reference['v'])
    seqs = queries['cdr3'].unique(maintain_order=True).to_list()
    control_sizes = dict.fromkeys(seqs, len(control))
    if puncture_populations:
        # Bin zero retains exact draw multiplicity. Native histograms require
        # radius >=1; only bin zero is used for population puncturing.
        exact = control.edit_histogram_batch(seqs, first_hit.scope(1,0,0), threads, False)
        control_sizes = {s:len(control)-int(row[0]) for s,row in zip(seqs,exact)}
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
        rv=[gene_family(v) for v in rv]
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
            s, v = query['cdr3'], gene_family(query['v'])
            nr, nc_population = reference_sizes[s], control_sizes[s]
            total = 0.0 if nr and nc_population else None
            for hit in hits[s] if total is not None else []:
                if rs[hit.ref_id] == s:
                    continue
                ref_v = rv[hit.ref_id]
                if v != ref_v and (v,ref_v) not in v_similarities:
                    raise ValueError('precomputed V similarities do not cover the query/reference families')
                weight = 1.0 if v == ref_v else .25*v_similarities[(v,ref_v)]
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


def legacy_components(queries, reference, control, background, *, threads=1,
                      v_similarities=None):
    """Owning held-out NED, separate V/J/length odds and global first-hit gate.

    Reference junctions take their first union-level source representative. The
    caller supplies the owning raw query projection and sampled prior background.
    Exact junction hits are excluded, while population/model counts remain fixed.
    Empirical negative-tail/Sidak fusion and accuracy summaries are external.
    """
    import numpy as np
    from collections import Counter
    if threads < 1 or not len(control) or not background.height:
        raise ValueError('legacy components require positive threads and nonempty control/background')
    if queries['query_id'].null_count() or queries['query_id'].n_unique()!=queries.height:
        raise ValueError('query identities must be unique and non-null')
    reference=reference.unique('cdr3',keep='first',maintain_order=True)
    if not reference.height:
        raise ValueError('legacy components require a nonempty reference union')
    seqs=queries['cdr3'].unique(maintain_order=True).to_list()
    rs=reference['cdr3'].to_list()
    hits,unit,radii=_density_neighbours(seqs,rs,threads)
    counts=control_histograms(control,seqs,threads=threads,radii=radii)
    rv=[gene_family(v) for v in reference['v']]
    if v_similarities is None:
        v_similarities=precompute_v_similarities(queries['v'],reference['v'])
    epitopes=reference['epitope'].unique(maintain_order=True).to_list()
    ei={e:i for i,e in enumerate(epitopes)}
    labels=reference['epitope'].to_list()
    # Native sparse edges are reused for every raw observation of a junction.
    edges={s:[h for h in hits[s] if rs[h.ref_id]!=s] for s in seqs}
    penalty={s:np.fromiter((math.exp(-h.score/400.)
                          for h in edges[s]),dtype=np.float64) for s in seqs}
    ned=np.zeros((queries.height,len(epitopes)),dtype=np.float64)
    for i,(seq,v) in enumerate(zip(queries['cdr3'],queries['v'])):
        family=gene_family(v)
        weights=[]
        for h in edges[seq]:
            other=rv[h.ref_id]
            if family!=other and (family,other) not in v_similarities:
                raise ValueError('precomputed V similarities do not cover the query/reference families')
            weights.append(1. if family==other else .25*v_similarities[(family,other)])
        # add.at retains the native edge addition order of the owning scalar sum.
        np.add.at(ned[i],[ei[labels[h.ref_id]] for h in edges[seq]],
                  np.asarray(weights)*penalty[seq]/np.asarray([1+counts[seq][h.n_subs] for h in edges[seq]]))
    gates={}
    for seq,curve in zip(seqs,unit):
        n=int(curve[1]);c=int(counts[seq][1])
        result=(evalue_result(n,c,reference.height,len(control)) if n else
                first_hit.pvalue([],[],reference.height,len(control)))
        gates[seq]={'first_hit_radius':1 if n else None,'first_hit_n_target':n,
                    'first_hit_n_control':c if n else 0,'first_hit_E':result['E'],
                    'first_hit_p_enrichment':result['p_enrichment'],
                    'first_hit_gate':result['p_enrichment']<.001}
    qfeatures={'v_logodds':[gene_family(v) for v in queries['v']],
               'j_logodds':[gene_family(j) for j in queries['j']],
               'length_logodds':[len(s) for s in queries['cdr3']]}
    bfeatures={'v_logodds':[gene_family(v) for v in background['v']],
               'j_logodds':[gene_family(j) for j in background['j']],
               'length_logodds':[len(s) for s in background['cdr3']]}
    bc={name:Counter(values) for name,values in bfeatures.items()}
    pieces=[]
    for epitope in epitopes:
        ref=reference.filter(pl.col('epitope')==epitope)
        available=ref.height>=20  # exact owning build_models eligibility
        rfeatures={'v_logodds':[gene_family(v) for v in ref['v']],
                   'j_logodds':[gene_family(j) for j in ref['j']],
                   'length_logodds':[len(s) for s in ref['cdr3']]}
        columns=[]
        for name,values in rfeatures.items():
            rc=Counter(values);bb=bc[name]
            nr,nb=ref.height+len(rc),background.height+len(bb)
            odds={value:(math.log if name=='length_logodds' else np.log)(((rc[value]+1)/nr)/((bb[value]+1)/nb))
                  for value in rc.keys()|bb.keys()}
            columns.append(pl.Series(name,[odds.get(value,0.) if available else None
                for value in qfeatures[name]],dtype=pl.Float64))
        frame=queries.select('query_id').with_columns(pl.lit(epitope,dtype=pl.String).alias('epitope'),
            pl.Series('ned_score',ned[:,ei[epitope]]),*columns,
            pl.lit(ref.height,dtype=pl.UInt64).alias('n_reference'),
            pl.lit(reference.height,dtype=pl.UInt64).alias('reference_union_size'),
            pl.lit(len(control),dtype=pl.UInt64).alias('control_size'),
            pl.lit(background.height,dtype=pl.UInt64).alias('prior_background_size'),
            pl.lit(available).alias('feature_model_available'))
        gate_fields={'first_hit_radius':pl.UInt32,'first_hit_n_target':pl.UInt64,
                     'first_hit_n_control':pl.UInt64,'first_hit_E':pl.Float64,
                     'first_hit_p_enrichment':pl.Float64,'first_hit_gate':pl.Boolean}
        frame=frame.with_columns(*[pl.Series(name,[gates[s][name] for s in queries['cdr3']],dtype=dtype)
                                  for name,dtype in gate_fields.items()])
        pieces.append(frame)
    return pl.concat(pieces).with_columns(pl.lit('legacy-heldout-components-v1').alias('estimator'))


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

    if not exclude_exact:
        # Preserve the accepted scalar division and inline null-fill operation order.
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


def unpaired_score(densities, priors, reference_size, *, log_scale=False):
    """Original score; log_scale puts dense targets on the same log-density units.

    Germline terms retain the original sparse-reference activation. Historical
    defaults and sparse scores are unchanged; dense ordering is unchanged.
    """
    density_values, prior_values = _score_series(densities), _score_series(priors)
    if isinstance(reference_size, pl.Series):
        sizes = _score_series(reference_size)
        if len(sizes) != len(density_values) or (sizes < 1).any() or (sizes != sizes.floor()).any():
            raise ValueError('aligned components and a positive integer reference size are required')
        if len(density_values) != len(prior_values) or (density_values < 0).any():
            raise ValueError('aligned components and nonnegative density are required')
        return pl.DataFrame({'density': density_values, 'prior': prior_values, 'size': sizes}).select(
            pl.when(pl.col('size') < 500).then((pl.col('density') + 1e-6).log() + pl.col('prior'))
            .otherwise((pl.col('density')+1e-6).log() if log_scale else pl.col('density'))
            .alias('score'))['score']
    if len(density_values) != len(prior_values) or reference_size < 1:
        raise ValueError('aligned components and a positive reference size are required')
    if (density_values < 0).any():
        raise ValueError('density must be nonnegative')
    score = (density_values+1e-6).log() if log_scale or reference_size<500 else density_values
    if reference_size<500:
        score=score+prior_values
    return score.rename('score')


def paired_score(beta_densities, alpha_densities, combined_priors, *, log_scale=False):
    """Original equal-weight cohort rank sum; ties receive mean zero-based rank.

    Inputs are aligned observations, retaining repeated junctions and pair context.
    ``combined_priors`` is beta V/J/length plus alpha V/J. Optional log_scale
    adds original chain log-densities and this same prior, without cohort ranks.
    It is factorised original evidence, not calibrated joint binding probability.
    """
    beta, alpha, prior = map(_score_series, (beta_densities, alpha_densities, combined_priors))
    if len(beta) != len(alpha) or len(beta) != len(prior):
        raise ValueError('paired components must describe the same aligned observations')
    if (beta < 0).any() or (alpha < 0).any():
        raise ValueError('density must be nonnegative')
    if log_scale:
        return ((beta+1e-6).log()+(alpha+1e-6).log()+prior).rename('score')
    return (beta.rank('average') + alpha.rank('average') + prior.rank('average') - 3).rename('score')
