"""Opt-in residue/position/length and gap/length multipliers over a core matrix."""
from pathlib import Path
import json
import math
import numpy as np
import polars as pl
from ..db.cache import sha256


def load_matrix_bundle(path):
    """Read hash-pinned model tables; missing cells retain the core's neutral weight.

    Coordinates and lengths include conserved junction anchors. Species and chain
    tables are separate. No query-specific fitting or benchmark labels are read.
    """
    root=Path(path).expanduser().resolve()
    manifest=json.loads((root/'manifest.json').read_text())
    if not isinstance(manifest,dict) or manifest.get('format')!='vdjam-bundle-v1':
        raise ValueError('unsupported matrix bundle format')
    residue,gaps={},{}
    for name,columns in [('residue.tsv',['species','locus','length','position','aa','multiplier']),
                         ('gaps.tsv',['species','locus','length','gap','multiplier'])]:
        if manifest.get('outputs_sha256',{}).get(name)!=sha256(root/name):
            raise ValueError(f'matrix bundle hash mismatch: {name}')
        frame=pl.read_csv(root/name,separator='\t',infer_schema_length=None)
        if frame.columns!=columns or frame.null_count().sum_horizontal().item():
            raise ValueError(f'matrix bundle requires complete {columns} columns')
        seen=set()
        for row in frame.iter_rows(named=True):
            species,locus,length=row['species'],row['locus'],row['length']
            multiplier=row['multiplier']
            if species not in ('human','mouse') or locus not in ('TRA','TRB'):
                raise ValueError('matrix bundle species/locus is unsupported')
            if type(length) is not int or not 2<=length<=65535:
                raise ValueError('matrix bundle requires an integer junction length in [2,65535]')
            if type(multiplier) not in (int,float) or not math.isfinite(multiplier) or multiplier<0:
                raise ValueError('matrix multipliers must be finite and nonnegative')
            key=tuple(row[c] for c in columns[:-1])
            if key in seen:raise ValueError('duplicate matrix bundle cell')
            seen.add(key)
            if name=='residue.tsv':
                position,aa=row['position'],row['aa']
                if type(position) is not int or not 0<=position<length or not isinstance(aa,str) or len(aa)!=1 or aa not in 'ACDEFGHIKLMNPQRSTVWY':
                    raise ValueError('invalid residue matrix coordinate')
                if multiplier>(2**31-1)/100:
                    raise ValueError('residue multiplier exceeds native integer range')
                context=(species,locus,length)
                if context not in residue:residue[context]=np.ones((length,128),dtype=float)
                table=residue[context]
                table[position,ord(aa)]=multiplier
            else:
                gap=row['gap']
                if type(gap) is not int or gap<1 or length+gap>65535:
                    raise ValueError('gap matrix requires a positive integer gap within native length range')
                gaps[(species,locus,length,gap)]=float(multiplier)
    return dict(residue=residue,gaps=gaps,manifest=manifest,path=str(root),
                manifest_sha256=sha256(root/'manifest.json'))


def sequence_weights(bundle,species,locus,sequences,ctrim):
    """Vectorized residue lookup per length; rows use native integer scale 100."""
    lengths=np.array([len(s) for s in sequences]);rows=[None]*len(sequences)
    for length in np.unique(lengths):
        indices=np.flatnonzero(lengths==length)
        profile=bundle['residue'].get((species,locus,int(length)))
        if profile is None:weights=np.full((len(indices),length),100,dtype=np.int64)
        else:
            codes=np.frombuffer(''.join(sequences[i] for i in indices).encode('ascii'),dtype=np.uint8).reshape(-1,length)
            weights=np.rint(100*profile[np.arange(length)[None,:],codes]).astype(np.int64)
        weights[:,:3]=0
        if ctrim:weights[:,-ctrim:]=0
        for i,row in zip(indices,weights.tolist()):rows[i]=row
    return rows


def gap_prior(bundle,species,locus,baseline,gap_open,gap_extend):
    """Preserve placement constraints; apply the learned multiplier to block cost."""
    from seqtree import gapblock
    def prior(i,d,m):
        if not d:return 0
        placement=baseline(i,d,m) if baseline is not None else 0
        multiplier=bundle['gaps'].get((species,locus,m-d,d),1.)
        value=placement+round(gapblock.gap_cost(d,gap_open,gap_extend)*multiplier)
        if value>2**31-1:
            raise ValueError('matrix bundle gap cost exceeds native integer range')
        return value
    return prior


def apply_bundle(options,bundle,species,lanes,ctrim):
    """Apply single or paired lanes while retaining the native grouped scorer."""
    opening,extension=options['gap_open'],options['gap_extend']
    for locus,q,r,suffix in lanes:
        options.pop('position_weights_by_length'+suffix,None)
        options['query_position_weights'+suffix]=sequence_weights(bundle,species,locus,q,ctrim)
        options['reference_position_weights'+suffix]=sequence_weights(bundle,species,locus,r,ctrim)
        options['gap_prior'+suffix]=gap_prior(bundle,species,locus,
            options.get('gap_prior'+suffix,options.get('gap_prior')),opening,extension)
    options['gap_open']=options['gap_extend']=0
    if len(lanes)>1:options['gap_prior']=None
