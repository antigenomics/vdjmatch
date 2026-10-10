"""TCRdist3 default human alpha/beta distance using existing native gap-block scoring.

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


def load_v_loops():
    """Bundled, source-pinned human aligned loop templates keyed by exact allele."""
    path=resources.files('vdjmatch.resources')/'tcrdist'/'human_v_loops.tsv'
    with path.open() as f:
        return {r['allele']:r for r in csv.DictReader(f,delimiter='\t')}


def distance_matrix(queries, references, query_v, reference_v, *, threads=1):
    """Human TCRdist3 default distances in a bounded native batch, maximum64MiB output.

    Exact alleles are used; a missing allele suffix explicitly means *01, as in the
    comparator. Unknown alleles, cross-locus inputs and junctions shorter than8 fail
    instead of silently losing gene evidence. Dots in aligned loops carry the original
    TCRdist3 unknown-symbol cost zero. Callers must handle missing calls explicitly.
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
    model=load_v_loops()
    genes=[v if isinstance(v,str) and '*' in v else str(v)+'*01' for v in query_v+reference_v]
    if any(v not in model for v in genes):
        raise ValueError('unknown human V allele; handle missing gene evidence before scoring')
    if len({model[v]['locus'] for v in genes})>1:
        raise ValueError('distance batch must contain one locus')
    if not queries or not references:
        return np.empty((len(queries),len(references)),dtype=np.int32)
    alphabet=amino_acids(); b=SubstitutionMatrix.blosum62()
    d=np.asarray([[0 if a==c else max(0,min(4,4-b.similarity(a,c))) for c in alphabet] for a in alphabet],dtype=np.int32)
    # Integer similarities -d produce exactly2*d under seqtree's Gram transform.
    matrix=SubstitutionMatrix.from_similarity((-d).tolist())
    out=3*np.asarray(gapblock.score_matrix([s[3:-2] for s in queries],[s[3:-2] for s in references],
        matrix=matrix,gap_open=8,gap_extend=8,gap_prior=_gap_prior,threads=threads))//2
    # TCRdist's unknown/dot row is zero. Loop templates all have28 aligned positions.
    lut=np.zeros((128,128),dtype=np.int32)
    for a in canonical:
        for c in canonical:lut[ord(a),ord(c)]=d[alphabet.index(a),alphabet.index(c)]
    loops=[model[v]['cdr1']+model[v]['cdr2']+model[v]['cdr25'] for v in genes]
    codes=np.frombuffer(''.join(loops).encode('ascii'),dtype=np.uint8).reshape(len(genes),28)
    qa,ra=codes[:len(queries)],codes[len(queries):]
    for pos in range(28):out+=lut[qa[:,None,pos],ra[None,:,pos]]
    return out
