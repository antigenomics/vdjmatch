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
    if not gap:
        return 0
    short = width-gap+5  # restore the three N- and two C-terminal trimmed residues
    low, high = 5, short-5
    while low > high:
        low -= 1
        high += 1
    return 0 if low-3 <= i <= high-3 else gapblock.UNREACHABLE


def load_v_loops(species='human'):
    """Bundled, source-pinned aligned loop templates keyed by exact allele."""
    if species not in ('human','mouse'):
        raise ValueError('TCRdist species must be human or mouse')
    path=resources.files('vdjmatch.resources')/'tcrdist'/f'{species}_v_loops.tsv'
    with path.open() as f:
        return {r['allele']:r for r in csv.DictReader(f,delimiter='\t')}


def _substitution_costs():
    alphabet=amino_acids(); b=SubstitutionMatrix.blosum62()
    return np.asarray([[0 if a==c else max(0,min(4,4-b.similarity(a,c)))
                        for c in alphabet] for a in alphabet],dtype=np.int32)


def _cdr3_options(threads=1, costs=None):
    """Native CDR3 geometry; native scores equal2/3 of weighted TCRdist CDR3."""
    costs=_substitution_costs() if costs is None else costs
    return dict(matrix=SubstitutionMatrix.from_similarity((-costs).tolist()),
                gap_open=8,gap_extend=8,gap_prior=_gap_prior,threads=threads)


def distance_matrix(queries, references, query_v, reference_v, *, species='human', threads=1, return_cdr3=False):
    """TCRdist3 default distances in a bounded native batch, maximum64MiB output.

    Exact alleles are used; a missing allele suffix explicitly means *01, as in the
    comparator. Unknown alleles, cross-locus inputs and junctions shorter than8 fail
    instead of silently losing gene evidence. Dots in aligned loops carry the original
    TCRdist3 unknown-symbol cost zero (also for stars). Human templates are28 aligned
    positions; mouse TRA templates are29. Callers handle missing calls explicitly.
    ``return_cdr3=True`` returns total and weighted CDR3 matrices separately.
    """
    if isinstance(threads,bool) or not isinstance(threads,int) or threads<1:
        raise ValueError('threads must be positive')
    queries, references, query_v, reference_v = map(list,(queries,references,query_v,reference_v))
    if len(queries)!=len(query_v) or len(references)!=len(reference_v):
        raise ValueError('each junction requires a V call')
    if 4*len(queries)*len(references)>64*1024**2:
        raise ValueError('distance output exceeds64MiB; use contiguous query batches')
    canonical=set('ACDEFGHIKLMNPQRSTVWY')
    if any(not isinstance(s,str) or len(s)<8 or not set(s)<=canonical for s in queries+references):
        raise ValueError('TCRdist junctions require canonical amino acids and length>=8')
    model=load_v_loops(species)
    genes=[v if isinstance(v,str) and '*' in v else str(v)+'*01' for v in query_v+reference_v]
    if any(v not in model for v in genes):
        raise ValueError(f'unknown {species} V allele; handle missing gene evidence before scoring')
    if len({model[v]['locus'] for v in genes})>1:
        raise ValueError('distance batch must contain one locus')
    if not queries or not references:
        out=np.empty((len(queries),len(references)),dtype=np.int32)
        return (out,out.copy()) if return_cdr3 else out
    alphabet=amino_acids(); d=_substitution_costs()
    # Integer similarities -d produce exactly2*d under seqtree's Gram transform.
    cdr3=np.array(gapblock.score_matrix([s[3:-2] for s in queries],[s[3:-2] for s in references],
        **_cdr3_options(threads,d)),copy=True)
    cdr3//=2
    cdr3*=3
    out=cdr3.copy() if return_cdr3 else cdr3
    # TCRdist's unknown/dot/star row is zero. Score unique allele pairs once per call.
    lut=np.zeros((128,128),dtype=np.int32)
    for a in canonical:
        for c in canonical:lut[ord(a),ord(c)]=d[alphabet.index(a),alphabet.index(c)]
    alleles,indices=np.unique(genes,return_inverse=True)
    loops=[model[v]['cdr1']+model[v]['cdr2']+model[v]['cdr25'] for v in alleles]
    width=len(loops[0])
    codes=np.frombuffer(''.join(loops).encode('ascii'),dtype=np.uint8).reshape(len(alleles),width)
    loop_dist=np.zeros((len(alleles),len(alleles)),dtype=np.int32)
    for pos in range(width):loop_dist+=lut[codes[:,None,pos],codes[None,:,pos]]
    out+=loop_dist[indices[:len(queries),None],indices[None,len(queries):]]
    return (out,cdr3) if return_cdr3 else out


def paired_distance_matrix(query_alpha, query_beta, reference_alpha, reference_beta,
                           query_v_alpha, query_v_beta, reference_v_alpha, reference_v_beta,
                           *, species='human', threads=1):
    """Sum alpha/beta distances against the same ordered reference pairs.

    Each position identifies one linked receptor pair; independent marginal nearest
    neighbours must not be substituted for the reference pair distance.
    """
    parts=list(map(list,(query_alpha,query_beta,reference_alpha,reference_beta,
                        query_v_alpha,query_v_beta,reference_v_alpha,reference_v_beta)))
    qa,qb,ra,rb,qva,qvb,rva,rvb=parts
    if len({len(qa),len(qb),len(qva),len(qvb)})>1 or len({len(ra),len(rb),len(rva),len(rvb)})>1:
        raise ValueError('each linked pair requires alpha and beta junctions and V calls')
    if any(not isinstance(v,str) or not v.startswith('TRAV') for v in qva+rva) or any(
            not isinstance(v,str) or not v.startswith('TRBV') for v in qvb+rvb):
        raise ValueError('linked pair V calls must have their declared alpha/beta locus')
    out=distance_matrix(qa,ra,qva,rva,species=species,threads=threads)
    out+=distance_matrix(qb,rb,qvb,rvb,species=species,threads=threads)
    return out
