"""Operation-local numeric components; expand source observations only in batches."""
from dataclasses import dataclass

import numpy as np
import polars as pl

from ..aggregate.candidates import PMHC

_BATCH_ROWS = 65536
_ARRAY_LIMIT_BYTES = 12 * 1024**3


@dataclass
class ComponentBank:
    mapping: pl.DataFrame
    targets: pl.DataFrame
    values: dict

    def frames(self, queries=None):
        """Reproduce requested source IDs, including repeated receptor observations."""
        mapping = self.mapping
        if queries is not None:
            mapping = queries.select('query_id').join(mapping,on='query_id',how='left',
                validate='1:1',maintain_order='left')
            if mapping['_query_row'].null_count():
                raise ValueError('unknown component query identities')
        groups = self.targets.height
        if not groups:
            return
        for part in mapping.iter_slices(max(1,_BATCH_ROWS//groups)):
            rows = part['_query_row'].to_numpy()
            yield part.drop('_query_row').join(self.targets,how='cross').with_columns(*[
                pl.Series(name,value[rows].reshape(-1),nan_to_null=True)
                for name,value in self.values.items()])


def score_component_bank(queries,reference,control,*,species,locus,threads=1,exclude_exact=False):
    """Score every label once on deduplicated raw receptor configurations.

    Raw arrays use float64 scores and Boolean support. No source label, Pgen
    stratum or assay outcome changes a model or scoring key. Source identities
    and repeated null draws are restored by ``frames``. The allocation guard
    includes paired temporary factors and bounded arithmetic scratch; it does
    not certify whole-process RSS or the manuscript resource gate.
    """
    from . import combined
    from ..evalue.control import _organism
    from ..evalue.paired import build_paired_ref
    from ..match.geometry import precompute_v_loop_distances
    from ..match.historical import _prepare_prior_background,paired_score

    if threads<1 or locus not in {'TRA','TRB','paired'}:
        raise ValueError('component scoring requires positive threads and TRA/TRB/paired locus')
    if not isinstance(control,dict) or not {'density_index','prior_background'}<=control.keys():
        raise ValueError('accepted internal density_index/prior_background model context required')
    if queries['query_id'].null_count() or queries['query_id'].n_unique()!=queries.height:
        raise ValueError('component query identities must be unique and non-null')
    paired = locus=='paired'
    if paired and (not isinstance(control['density_index'],dict) or
                   not {'TRA','TRB'}<=control['density_index'].keys()):
        raise ValueError('paired accepted density context requires TRA and TRB indexes')
    reference=combined._reference_population(reference,_organism(species))
    reference=reference.filter(pl.col('gene').is_in(['TRA','TRB'] if paired else [locus]))
    targets=reference.select(PMHC).unique(maintain_order=True)
    keys=['cdr3a','cdr3b','va','vb','ja','jb'] if paired else ['cdr3','v','j']
    unique=queries.select(keys).unique(maintain_order=True).with_row_index('query_id')
    mapping=queries.select('query_id',*keys).join(unique.rename({'query_id':'_query_row'}),
        on=keys,how='left',nulls_equal=True,validate='m:1',maintain_order='left').select('query_id','_query_row')
    shape=(unique.height,targets.height)
    required=shape[0]*shape[1]*(80 if paired else 32)+_BATCH_ROWS*8*32
    if required>_ARRAY_LIMIT_BYTES or targets.height>_BATCH_ROWS:
        raise ValueError('component array budget exceeded; reduce the receptor population')
    prior_counts=control.get('_prior_counts')
    if prior_counts is None:
        prior_counts=_prepare_prior_background(control['prior_background'])
    contexts=control.get('v_loop_distances')
    if contexts is None:
        contexts={matrix:precompute_v_loop_distances(species=species,alphabet='full',matrix=matrix,
            v_regions='cdr12plus') for matrix in ('clipped','blosum62')}
    values={}
    def value(name):
        if name not in values:
            boolean=name.startswith('original_supported')
            values[name]=np.full(shape,False if boolean else np.nan,dtype=bool if boolean else np.float64)
        return values[name]
    target_rows=targets.with_row_index('_bank_target')
    def receiver(columns,kernel=False):
        seen=np.zeros(shape,dtype=bool)
        def consume(frame):
            frame=frame.join(target_rows,on=PMHC,how='left',nulls_equal=True,validate='m:1')
            if frame['_bank_target'].null_count():
                raise RuntimeError('component batch contains unknown targets')
            q=frame['query_id'].to_numpy();t=frame['_bank_target'].to_numpy()
            if np.any(q>=shape[0]) or seen[q,t].any() or np.unique(q.astype(np.uint64)*shape[1]+t).size!=len(q):
                raise RuntimeError('component batch contains unknown or repeated cases')
            seen[q,t]=True
            for source,dest in columns:
                raw=frame[source].to_numpy()
                value(dest)[q,t]=np.where(raw>0,raw,np.nan) if kernel else raw
        return consume
    lanes=[('TRA','a','alpha'),('TRB','b','beta')] if paired else [(locus,'','')]
    for chain,suffix,name in lanes:
        query=(unique.select('query_id',*[pl.col(c+suffix).alias(c) for c in ('cdr3','v','j')])
            if paired else unique).with_columns(pl.lit(chain).alias('locus'))
        ref=reference.filter(pl.col('gene')==chain)
        columns=[('original_score','original_'+name+'_score' if paired else 'original_score'),
            ('original_supported','original_supported_'+name if paired else 'original_supported')]
        if paired:
            columns += [('original_density','_'+name+'_density'),('germline_lr','_'+name+'_prior')]
        consume=receiver(columns)
        combined._original(query,ref,None,control['density_index'][chain] if paired else control['density_index'],
            control['prior_background'],locus=chain,threads=threads,exclude_exact=exclude_exact,
            prior_counts=prior_counts,_consume=consume)
        del consume
    if paired:
        factors=[value('_'+name+'_'+kind).reshape(-1) for name,kind in
            [('beta','density'),('alpha','density'),('alpha','prior'),('beta','prior')]]
        result=value('original_score').reshape(-1)
        part=None
        for start in range(0,len(result),_BATCH_ROWS):
            part=[x[start:start+_BATCH_ROWS] for x in factors]
            valid=np.isfinite(np.column_stack(part)).all(axis=1)
            if valid.any():
                beta,alpha,pa,pb=[pl.Series(x[valid]) for x in part]
                indices=np.flatnonzero(valid)+start
                result[indices]=paired_score(beta,alpha,pa+pb,log_scale=True).to_numpy()
        values['original_supported']=value('original_supported_alpha') & value('original_supported_beta')
        del factors,result,part
        for name in ('_alpha_density','_alpha_prior','_beta_density','_beta_prior'):
            del values[name]
    for chain,suffix,name in lanes:
        query=(unique.select('query_id',*[pl.col(c+suffix).alias(c) for c in ('cdr3','v','j')])
            if paired else unique).with_columns(pl.lit(chain).alias('locus'))
        ref=reference.filter(pl.col('gene')==chain)
        consumers={}
        def consume(kind,frame):
            if kind not in consumers:
                consumers[kind]=receiver([('kernel_score',kind+'_'+name+'_score' if paired else kind+'_score')],kernel=True)
            consumers[kind](frame)
        combined._kernel(query,ref,None,species=species,threads=threads,paired=False,
            exclude_exact=exclude_exact,contexts=contexts,_consume=consume)
        del consumers,consume
    if paired:
        if 'complex_id' not in reference.columns:
            reference=reference.with_columns(pl.lit('0').alias('complex_id'))
        consumers={}
        def consume(kind,frame):
            if kind not in consumers:
                consumers[kind]=receiver([('kernel_score',kind+'_linked_score')],kernel=True)
            consumers[kind](frame)
        combined._kernel(unique,build_paired_ref(reference),None,species=species,threads=threads,paired=True,
            exclude_exact=exclude_exact,contexts=contexts,_consume=consume)
        del consumers,consume
    names=(['original_score','original_supported',*[c+'_'+axis+'_score' for axis in ('alpha','beta')
            for c in ('original','gapblock','stem4')],*[c+'_linked_score' for c in ('gapblock','stem4')],
            'original_supported_alpha','original_supported_beta'] if paired else
           ['original_score','original_supported','gapblock_score','stem4_score'])
    for name in names:
        value(name)
    return ComponentBank(mapping,targets,values)
