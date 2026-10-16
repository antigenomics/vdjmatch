"""Operation-local sparse components; reconstruct fixed priors only in batches."""
from dataclasses import dataclass

import numpy as np
import polars as pl

from ..aggregate.candidates import PMHC

_BATCH_ROWS = 65536
_ARRAY_LIMIT_BYTES = 12 * 1024**3


def _mark(seen,rows,targets,n_queries,n_targets):
    if np.any((rows<0)|(rows>=n_queries)) or np.any((targets<0)|(targets>=n_targets)):
        raise RuntimeError('randomized scoring contains unknown cases')
    codes=rows.astype(np.uint64)*n_targets+targets
    if len(np.unique(codes))!=len(codes):
        raise RuntimeError('randomized scoring repeats cases')
    offsets=(codes//8).astype(np.intp);bits=np.left_shift(np.uint8(1),(codes%8).astype(np.uint8))
    if np.any(seen[offsets]&bits):
        raise RuntimeError('randomized scoring repeats cases')
    np.bitwise_or.at(seen,offsets,bits)


@dataclass
class _SparseScores:
    codes: np.ndarray
    scores: np.ndarray
    groups: int
    default: float = np.nan

    def take(self,rows):
        codes=(rows.astype(np.uint64)[:,None]*self.groups+np.arange(self.groups,dtype=np.uint64)).reshape(-1)
        positions=np.searchsorted(self.codes,codes)
        found=np.flatnonzero(positions<len(self.codes))
        found=found[self.codes[positions[found]]==codes[found]]
        out=np.full(len(codes),self.default)
        out[found]=self.scores[positions[found]]
        return out


def _sparse(parts,groups,default=np.nan):
    codes=np.concatenate([p[0] for p in parts]) if parts else np.empty(0,dtype=np.uint64)
    scores=np.concatenate([p[1] for p in parts]) if parts else np.empty(0)
    parts.clear()
    order=np.argsort(codes)
    return _SparseScores(codes[order],scores[order],groups,default)


@dataclass
class _OriginalModel:
    density: _SparseScores
    sizes: np.ndarray
    available: np.ndarray
    features: list
    background_size: int

    def evaluate(self,rows):
        from ..match.historical import _prior_log_ratio,unpaired_score
        groups=len(self.sizes)
        density=self.density.take(rows)
        sizes=np.tile(self.sizes,len(rows))
        valid=np.repeat(self.available[rows],groups)&(sizes>0)
        prior=np.zeros(len(sizes))
        for codes,counts,background in self.features:
            selected=codes[rows]
            prior=pl.DataFrame({'prior':prior,'reference':counts[selected].reshape(-1),
                'background':np.repeat(background[selected],groups),'size':sizes}).select(
                    (pl.col('prior')+_prior_log_ratio(pl.col('reference'),pl.col('background'),
                        pl.col('size'),self.background_size)).alias('prior'))['prior'].to_numpy().copy()
        prior[~valid]=np.nan
        score=np.full(len(sizes),np.nan)
        if valid.any():
            score[valid]=unpaired_score(density[valid],prior[valid],pl.Series(sizes[valid]),log_scale=True).to_numpy()
        return density,prior,score,density>0


def _original_model(queries,reference,targets,background,include_length,reserve):
    """Fixed reference/control marginal counts; no query-specific model changes."""
    from ..match.historical import _prior_features
    features=_prior_features(queries)
    ref=_prior_features(reference.unique([*PMHC,'cdr3'],keep='first',maintain_order=True)).join(
        targets,on=PMHC,how='left',nulls_equal=True,validate='m:1')
    groups=targets.height
    reserve(groups*8+queries.height)
    sizes=np.zeros(groups,dtype=np.uint64)
    counts=ref.group_by('_bank_target').len()
    sizes[counts['_bank_target'].to_numpy()]=counts['len'].to_numpy()
    available=queries.select((pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$') &
        pl.all_horizontal(pl.col('v','j').cast(pl.String).str.strip_chars().ne(''))).fill_null(False)).to_series().to_numpy()
    prepared=[]
    nb,marginals=background
    for column in ['v','j']+(['length'] if include_length else []):
        axis=features.select(column).unique(maintain_order=True).with_row_index('_feature')
        reserve(queries.height*4+axis.height*(groups*4+8))
        codes=features.select(column).join(axis,on=column,how='left',nulls_equal=True,
            validate='m:1',maintain_order='left')['_feature'].to_numpy()
        matrix=np.zeros((axis.height,groups),dtype=np.uint32)
        # Ordinary joins intentionally do not match null feature values, as in germline_prior.
        rc=ref.group_by('_bank_target',column).len().join(axis,on=column,how='inner',validate='m:1')
        matrix[rc['_feature'].to_numpy(),rc['_bank_target'].to_numpy()]=rc['len'].to_numpy()
        bc=axis.join(marginals[column],on=column,how='left',validate='1:1',maintain_order='left')
        prepared.append((codes,matrix,bc['_bg_count'].fill_null(0).cast(pl.UInt64).to_numpy()))
    return _OriginalModel(_sparse([],groups,0.),sizes,available,prepared,nb)


@dataclass
class ComponentBank:
    mapping: pl.DataFrame
    targets: pl.DataFrame
    values: dict
    originals: dict | None = None

    def frames(self, queries=None):
        """Reproduce every source case, including zero-density continuous priors."""
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
            rows=part['_query_row'].to_numpy()
            data={name:(value.take(rows) if isinstance(value,_SparseScores) else value[rows].reshape(-1))
                for name,value in self.values.items() if value is not None}
            if self.originals is not None:
                original={name:model.evaluate(rows) for name,model in self.originals.items()}
                for name,(_,_,score,support) in original.items():
                    data['original_'+name+'_score' if name else 'original_score']=score
                    data['original_supported_'+name if name else 'original_supported']=support
                if 'alpha' in original:
                    from ..match.historical import paired_score
                    da,pa,_,sa=original['alpha'];db,pb,_,sb=original['beta']
                    valid=np.isfinite(pa)&np.isfinite(pb)
                    paired=np.full(len(pa),np.nan)
                    if valid.any():
                        paired[valid]=paired_score(pl.Series(db[valid]),pl.Series(da[valid]),
                            pl.Series(pa[valid])+pl.Series(pb[valid]),log_scale=True).to_numpy()
                    data['original_score']=paired;data['original_supported']=sa&sb
            yield part.drop('_query_row').join(self.targets,how='cross').with_columns(*[
                pl.Series(name,data[name],nan_to_null=True) for name in self.values])


def score_component_bank(queries,reference,control,*,species,locus,threads=1,exclude_exact=False):
    """Collect positive raw support and reconstruct unchanged fixed priors.

    Missing kernel evidence remains missing; original density defaults to zero.
    Original priors use fixed per-feature reference/control counts, with the
    original floating-point operation order and sparse-reference activation.
    Every source observation is restored in bounded full-label batches. The
    12 GiB allocation guard includes sparse-record copies and arithmetic scratch;
    it does not establish whole-process RSS or the manuscript resource gate.
    """
    from . import combined
    from ..evalue.control import _organism
    from ..evalue.paired import build_paired_ref
    from ..match.geometry import precompute_v_loop_distances
    from ..match.historical import _prepare_prior_background

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
    coverage=(shape[0]*shape[1]+7)//8
    used=0
    def reserve(size):
        nonlocal used
        used+=size
        if used>_ARRAY_LIMIT_BYTES or targets.height>_BATCH_ROWS:
            raise ValueError('component array budget exceeded; reduce the receptor population')
    reserve(2*coverage+max(_BATCH_ROWS*8*32,unique.height*8*32))
    prior_counts=control.get('_prior_counts')
    if prior_counts is None:
        prior_counts=_prepare_prior_background(control['prior_background'])
    contexts=control.get('v_loop_distances')
    if contexts is None:
        contexts={matrix:precompute_v_loop_distances(species=species,alphabet='full',matrix=matrix,
            v_regions='cdr12plus') for matrix in ('clipped','blosum62')}
    values={};originals={};parts={}
    target_rows=targets.with_row_index('_bank_target')
    def receiver(column,name):
        seen=np.zeros(coverage,dtype=np.uint8)
        collected=parts.setdefault(name,[])
        def consume(frame):
            frame=frame.join(target_rows,on=PMHC,how='left',nulls_equal=True,validate='m:1')
            if frame['_bank_target'].null_count():
                raise RuntimeError('component batch contains unknown targets')
            q=frame['query_id'].to_numpy();t=frame['_bank_target'].to_numpy()
            _mark(seen,q,t,*shape)
            raw=frame[column].to_numpy();keep=raw>0
            if keep.any():
                reserve(int(keep.sum())*48)  # retained records plus packing/sort copies
                collected.append((q[keep].astype(np.uint64)*shape[1]+t[keep],raw[keep]))
        return consume
    lanes=[('TRA','a','alpha'),('TRB','b','beta')] if paired else [(locus,'','')]
    for chain,suffix,name in lanes:
        query=(unique.select('query_id',*[pl.col(c+suffix).alias(c) for c in ('cdr3','v','j')])
            if paired else unique).with_columns(pl.lit(chain).alias('locus'))
        ref=reference.filter(pl.col('gene')==chain)
        model=_original_model(query,ref,target_rows,prior_counts,chain=='TRB',reserve)
        consume=receiver('original_density',name)
        combined._original(query,ref,None,control['density_index'][chain] if paired else control['density_index'],
            control['prior_background'],locus=chain,threads=threads,exclude_exact=exclude_exact,
            prior_counts=prior_counts,_consume=consume)
        del consume
        model.density=_sparse(parts.pop(name,[]),shape[1],0.)
        originals[name]=model
    for chain,suffix,name in lanes:
        query=(unique.select('query_id',*[pl.col(c+suffix).alias(c) for c in ('cdr3','v','j')])
            if paired else unique).with_columns(pl.lit(chain).alias('locus'))
        ref=reference.filter(pl.col('gene')==chain)
        consumers={}
        def consume(kind,frame):
            key=kind+'_'+name+'_score' if paired else kind+'_score'
            if key not in consumers:consumers[key]=receiver('kernel_score',key)
            consumers[key](frame)
        combined._kernel(query,ref,None,species=species,threads=threads,paired=False,
            exclude_exact=exclude_exact,contexts=contexts,_consume=consume)
        del consumers,consume
    if paired:
        if 'complex_id' not in reference.columns:
            reference=reference.with_columns(pl.lit('0').alias('complex_id'))
        consumers={}
        def consume(kind,frame):
            key=kind+'_linked_score'
            if key not in consumers:consumers[key]=receiver('kernel_score',key)
            consumers[key](frame)
        combined._kernel(unique,build_paired_ref(reference),None,species=species,threads=threads,paired=True,
            exclude_exact=exclude_exact,contexts=contexts,_consume=consume)
        del consumers,consume
    names=(['original_score','original_supported',*[c+'_'+axis+'_score' for axis in ('alpha','beta')
            for c in ('original','gapblock','stem4')],*[c+'_linked_score' for c in ('gapblock','stem4')],
            'original_supported_alpha','original_supported_beta'] if paired else
           ['original_score','original_supported','gapblock_score','stem4_score'])
    for name in names:
        values[name]=None if name.startswith('original') else _sparse(parts.pop(name,[]),shape[1])
    return ComponentBank(mapping,targets,values,originals)
