"""TCRdist3 default human/mouse alpha/beta distance using native gap-block scoring.

Reproduces 3*CDR3 + aligned CDR1/CDR2/CDR2.5, including TCRdist3's clipped
BLOSUM matrix and restricted gap positions. This distance is not calibrated evidence.
"""
from __future__ import annotations

import csv
from importlib import resources

import numpy as np
from seqtree import SubstitutionMatrix, amino_acids, gapblock


def _gap_prior(i, gap, width):
    return _gap_window(i, gap, width, 5)


def _gap_prior_trim3(i, gap, width):
    return _gap_window(i, gap, width, 6)


def _gap_window(i, gap, width, trimmed, ntrim=3):
    if not gap:
        return 0
    short = width-gap+trimmed  # restore full junction coordinates
    low, high = 5, short-5
    while low > high:
        low -= 1
        high += 1
    return 0 if low-ntrim <= i <= high-ntrim else gapblock.UNREACHABLE


def _full_gap_prior(i, gap, width):
    return _gap_window(i, gap, width, 0, ntrim=0)


def load_v_loops(species='human'):
    """Bundled, source-pinned aligned loop templates keyed by exact allele."""
    if species not in ('human','mouse'):
        raise ValueError('TCRdist species must be human or mouse')
    path=resources.files('vdjmatch.resources')/'tcrdist'/f'{species}_v_loops.tsv'
    with path.open() as f:
        return {r['allele']:r for r in csv.DictReader(f,delimiter='\t')}


def resolve_v_alleles(calls, model, *, locus=None):
    """Resolve exact/model-proven co-locus spellings, preserving callers' raw calls.

    Omitted allele means *01. Only a unique modeled co-locus gene can supply an
    omitted /DV suffix or slash; no gene-family prefix or first comma call is chosen.
    Unknown alleles, missing inputs, cross-locus calls and ambiguous aliases are null.
    Mouse -DV and /DV spellings merge only when their aligned templates agree.
    """
    import polars as pl
    aliases={};exact_targets={}
    for key,meta in model.items():
        gene,allele=key.split('*')
        canonical=gene.replace('-DV','/DV')+'*'+allele
        if canonical not in model or any(meta[c]!=model[canonical][c] for c in ('locus','cdr1','cdr2','cdr25')):
            canonical=key
        exact_targets[key]=canonical
        forms={gene}
        if '/DV' in canonical:
            full=canonical.split('*')[0]
            forms.update([full,full.replace('/',''),full.split('/DV')[0]])
        for form in forms:
            tokens=[form+'*'+allele]+([form] if allele=='01' else [])
            for token in tokens:
                aliases.setdefault(token,set()).add(canonical)
    values=pl.Series(calls,dtype=pl.String)
    unique=values.drop_nulls().unique()
    resolved=[]
    for raw in unique:
        token=raw.strip()
        exact=token if '*' in token else token+'*01'
        candidates={exact_targets[exact]} if exact in model else aliases.get(token,set())
        target=next(iter(candidates)) if len(candidates)==1 else None
        resolved.append(target if target is not None and (locus is None or model[target]['locus']==locus) else None)
    return values.replace_strict(unique,pl.Series(resolved,dtype=pl.String),default=None,return_dtype=pl.String)


def _substitution_costs():
    alphabet=amino_acids(); b=SubstitutionMatrix.blosum62()
    return np.asarray([[0 if a==c else max(0,min(4,4-b.similarity(a,c)))
                        for c in alphabet] for a in alphabet],dtype=np.int32)


def _cdr3_options(threads=1, costs=None, ctrim=2):
    """Native CDR3 geometry; native scores equal2/3 of weighted TCRdist CDR3."""
    costs=_substitution_costs() if costs is None else costs
    return dict(matrix=SubstitutionMatrix.from_similarity((-costs).tolist()),
                gap_open=8,gap_extend=8,
                gap_prior=_gap_prior if ctrim==2 else _gap_prior_trim3,threads=threads)


def _distance_inputs(queries, references, query_v, reference_v, species, threads, ctrim, *, dense=False):
    """Shared junction/gene availability contract for dense distances and counts."""
    if isinstance(ctrim,bool) or not isinstance(ctrim,int) or ctrim not in (2,3):
        raise ValueError('ctrim must be 2 (TCRdist3) or 3 (symmetric trim)')
    if isinstance(threads,bool) or not isinstance(threads,int) or threads<1:
        raise ValueError('threads must be positive')
    queries, references, query_v, reference_v = map(list,(queries,references,query_v,reference_v))
    if len(queries)!=len(query_v) or len(references)!=len(reference_v):
        raise ValueError('each junction requires a V call')
    if dense and 4*len(queries)*len(references)>64*1024**2:
        raise ValueError('distance output exceeds64MiB; use contiguous query batches')
    canonical=set('ACDEFGHIKLMNPQRSTVWY')
    if any(not isinstance(s,str) or len(s)<8 or not set(s)<=canonical for s in queries+references):
        raise ValueError('TCRdist junctions require canonical amino acids and length>=8')
    model=load_v_loops(species)
    genes=resolve_v_alleles(query_v+reference_v,model).to_list()
    if any(v not in model for v in genes):
        raise ValueError(f'unknown {species} V allele; handle missing gene evidence before scoring')
    if len({model[v]['locus'] for v in genes})>1:
        raise ValueError('distance batch must contain one locus')
    return queries, references, model, genes


def _v_loop_groups(genes, model, d):
    """Small aligned-loop table plus per-input group IDs, using modeled alleles."""
    if not genes:
        return np.zeros((1,1),dtype=np.int32), np.empty(0,dtype=np.int32)
    # TCRdist's unknown/dot/star row is zero. Score unique allele pairs once per call.
    canonical=set('ACDEFGHIKLMNPQRSTVWY')
    alphabet=amino_acids()
    lut=np.zeros((128,128),dtype=np.int32)
    for a in canonical:
        for c in canonical:lut[ord(a),ord(c)]=d[alphabet.index(a),alphabet.index(c)]
    alleles,indices=np.unique(genes,return_inverse=True)
    loops=[model[v]['cdr1']+model[v]['cdr2']+model[v]['cdr25'] for v in alleles]
    width=len(loops[0])
    codes=np.frombuffer(''.join(loops).encode('ascii'),dtype=np.uint8).reshape(len(alleles),width)
    loop_dist=np.zeros((len(alleles),len(alleles)),dtype=np.int32)
    for pos in range(width):loop_dist+=lut[codes[:,None,pos],codes[None,:,pos]]
    return loop_dist, indices


def distance_matrix(queries, references, query_v, reference_v, *, species='human', threads=1, return_cdr3=False, ctrim=2):
    """TCRdist3 default distances in a bounded native batch, maximum64MiB output.

    Exact alleles are used; a missing allele suffix explicitly means *01, as in the
    comparator. Model-proven co-locus aliases resolve without changing raw calls.
    Unknown alleles, ambiguous names, cross-locus inputs and junctions shorter than8 fail
    instead of silently losing gene evidence. Dots in aligned loops carry the original
    TCRdist3 unknown-symbol cost zero (also for stars). Human templates are28 aligned
    positions; mouse TRA templates are29. Callers handle missing calls explicitly.
    ``return_cdr3=True`` returns total and weighted CDR3 matrices separately.
    ``ctrim=3`` tests removal of three residues at both junction ends; the default
    remains the TCRdist3 three/two trim. Full junction identity is handled by callers.
    """
    queries, references, model, genes = _distance_inputs(
        queries, references, query_v, reference_v, species, threads, ctrim, dense=True)
    if not queries or not references:
        out=np.empty((len(queries),len(references)),dtype=np.int32)
        return (out,out.copy()) if return_cdr3 else out
    d=_substitution_costs()
    # Integer similarities -d produce exactly2*d under seqtree's Gram transform.
    cdr3=np.array(gapblock.score_matrix([s[3:-ctrim] for s in queries],[s[3:-ctrim] for s in references],
        **_cdr3_options(threads,d,ctrim)),copy=True)
    cdr3//=2
    cdr3*=3
    out=cdr3.copy() if return_cdr3 else cdr3
    loop_dist, indices = _v_loop_groups(genes, model, d)
    out+=loop_dist[indices[:len(queries),None],indices[None,len(queries):]]
    return (out,cdr3) if return_cdr3 else out


def total_distance_count_batch(queries, references, query_v, reference_v, thresholds, *,
                               species='human', threads=1, ctrim=2,
                               position_weighting='uniform', exclude_exact=False):
    """Count full-junction/V-loop balls without a dense query/reference matrix.

    ``thresholds`` contains one integer list per query, in native units200 times
    the ordinary TCRdist distance. Results preserve threshold order/duplicates.
    ``uniform`` reproduces distance_matrix with three N-terminal and ctrim (2/3)
    C-terminal positions zero-weighted. ``significance`` instead uses the shipped
    full end-anchored profile, quantized max(1,round(100*w)), including flanks;
    it rejects ctrim3 rather than combining trimming and decay. Both retain the
    full-coordinate TCRdist moving gap window and linear gap cost12 per residue.

    The clipped substitution matrix contributes threefold CDR3 distance; aligned
    CDR1/CDR2/CDR2.5 contribute unweighted V-loop distance. Profile weights apply
    only to substitutions, in the symmetric longer-junction frame. Native units
    stay common: 6*clipped_cost*position_weight, gap2400*length_difference,
    V-offset200*loop_distance. Exact exclusion punctures identical full encoded
    junctions even for different V calls; duplicate references count separately.
    Unknown/ambiguous/mixed-locus inputs fail as in distance_matrix. This supplies
    geometry, not a background population, calibrated P-value or paired count.
    """
    import inspect
    if 'group_distances' not in inspect.signature(gapblock.count_batch).parameters:
        raise RuntimeError('total-distance counts require seqtree gapblock group_distances support')
    if position_weighting not in ('uniform','significance'):
        raise ValueError('position_weighting must be uniform or significance')
    if position_weighting=='significance' and ctrim!=2:
        raise ValueError('significance uses full positional decay; ctrim must remain default2')
    queries, references, model, genes = _distance_inputs(
        queries, references, query_v, reference_v, species, threads, ctrim)
    costs=_substitution_costs()
    loop_dist, indices=_v_loop_groups(genes,model,costs)
    lengths={len(s) for s in queries+references}
    if position_weighting=='significance':
        from .regions import significance_weights
        weights={length:[max(1,round(100*w)) for w in significance_weights(length)] for length in lengths}
    else:
        weights={length:[0]*3+[100]*(length-3-ctrim)+[0]*ctrim for length in lengths}
    return gapblock.count_batch(queries,references,thresholds,
        matrix=SubstitutionMatrix.from_similarity((-3*costs).tolist()),
        gap_open=2400,gap_extend=2400,gap_prior=_full_gap_prior,
        position_weights_by_length=weights,
        query_group_ids=indices[:len(queries)].tolist(),
        reference_group_ids=indices[len(queries):].tolist(),
        group_distances=(200*loop_dist).tolist(),threads=threads,exclude_exact=exclude_exact)


def paired_distance_matrix(query_alpha, query_beta, reference_alpha, reference_beta,
                           query_v_alpha, query_v_beta, reference_v_alpha, reference_v_beta,
                           *, species='human', threads=1, ctrim=2):
    """Sum alpha/beta distances against the same ordered reference pairs.

    Each position identifies one linked receptor pair; independent marginal nearest
    neighbours must not be substituted for the reference pair distance.
    """
    parts=list(map(list,(query_alpha,query_beta,reference_alpha,reference_beta,
                        query_v_alpha,query_v_beta,reference_v_alpha,reference_v_beta)))
    qa,qb,ra,rb,qva,qvb,rva,rvb=parts
    if len({len(qa),len(qb),len(qva),len(qvb)})>1 or len({len(ra),len(rb),len(rva),len(rvb)})>1:
        raise ValueError('each linked pair requires alpha and beta junctions and V calls')
    model=load_v_loops(species)
    qva,rva=(resolve_v_alleles(v,model,locus='TRA').to_list() for v in [qva,rva])
    qvb,rvb=(resolve_v_alleles(v,model,locus='TRB').to_list() for v in [qvb,rvb])
    if any(v is None for v in qva+rva+qvb+rvb):
        raise ValueError('linked pair V calls must resolve uniquely to their declared alpha/beta locus')
    out=distance_matrix(qa,ra,qva,rva,species=species,threads=threads,ctrim=ctrim)
    out+=distance_matrix(qb,rb,qvb,rvb,species=species,threads=threads,ctrim=ctrim)
    return out
