"""TCRnet-filtered flexible repertoire matching; descriptive counts, not repertoire P-values."""
import math
import polars as pl
from ..aggregate.candidates import PMHC
from ..io.columns import _integer_count


class TcrnetMatcher:
    """Own fixed background indexes and a once-filtered reference for a fresh cohort.

    TCRnet uses its existing unique-junction, locus-specific control test. Reference
    filtering families are individual pMHC/locus populations; sample families are
    individual sample/locus populations. Gene variants and pairing are restored
    after the junction filter. Available-chain matches count each receptor once.
    """
    def __init__(self,reference,controls,*,species='human',filter_fdr=.1,
                 scope='3,1,1,3',match_v=False,threads=1):
        from seqtree import Index
        from ..api import Annotator
        from ..db.schema import normalize
        from ..match.scope import parse_scope
        if not math.isfinite(filter_fdr) or not 0<filter_fdr<1 or type(threads) is not int or threads<1:
            raise ValueError('finite filter FDR in (0,1) and positive threads required')
        parse_scope(scope)
        self.species,self.filter_fdr,self.threads=species,filter_fdr,threads
        self.scope,self.match_v=scope,match_v
        reference=normalize(reference).filter(pl.col('reference_valid'))
        if not reference.height:raise ValueError('nonempty curated reference required')
        from ..evalue.control import _organism
        if reference['species'].null_count() or any(_organism(s)!=_organism(species) for s in reference['species'].unique()):
            raise ValueError('reference and control species must agree')
        if not set(reference['gene'])<={'TRA','TRB'}:raise ValueError('TRA/TRB reference loci required')
        self.inventory=reference.select(PMHC).unique().sort(PMHC,nulls_last=True)
        self.controls={}
        for locus,frame in controls.items():
            if locus not in {'TRA','TRB'} or not frame.height or 'junction_aa' not in frame.columns:
                raise ValueError('nonempty TRA/TRB junction controls required')
            sequences=frame['junction_aa'].unique().sort()
            if not sequences.str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False).all():
                raise ValueError('valid full-junction amino acid controls required')
            self.controls[locus]=Index.build(sequences.to_list(),alphabet='aa')
        retained=[];audits=[]
        for key,part in reference.partition_by([*PMHC,'gene'],as_dict=True).items():
            locus=key[-1]
            frame=part.select(pl.col('cdr3').alias('junction_aa'),'v','j').with_columns(pl.lit(1).alias('count'))
            result=self._filter(frame,locus).with_columns(*[pl.lit(value,dtype=reference.schema[name]).alias(name) for name,value in zip(PMHC,key[:-1])])
            audits.append(result)
            retained.append(part.join(result.filter(pl.col('retained')).select(pl.col('junction_aa').alias('cdr3')),on='cdr3',how='semi'))
        self.reference_filtering=pl.concat(audits,how='diagonal_relaxed')
        filtered=pl.concat(retained)
        self.reference_counts=self.inventory.join(filtered.group_by(PMHC).agg(
            pl.len().alias('reference_retained_records'),pl.struct('gene','cdr3','v','j').n_unique().alias('reference_retained_keys')),
            on=PMHC,how='left',nulls_equal=True).with_columns(pl.col('reference_retained_records','reference_retained_keys').fill_null(0))
        self.annotator=Annotator.from_frame(filtered)

    def _filter(self,frame,locus):
        from vdjtools.overlap import tcrnet
        if locus not in self.controls:raise ValueError(f'missing {locus} control index')
        canonical=frame.select('junction_aa',pl.col('v').alias('v_call'),pl.col('j').alias('j_call'),pl.col('count').alias('duplicate_count'))
        return tcrnet(canonical,control=self.controls[locus],locus=locus,species=self.species,
            scope='1,0,0,1',exclude_exact=True,threads=self.threads).with_columns(
                ((pl.col('q_value')<self.filter_fdr)&(pl.col('n_neighbors')>pl.col('E'))).alias('retained'))

    def match(self,observations):
        required=['receptor_id','locus','junction_aa','v','j','count']
        if not set(required)<=set(observations.columns) or not observations.height:
            raise ValueError('nonempty observations with receptor identities, chain keys and counts required')
        obs=observations.select(required).with_columns(_integer_count('count').alias('count'))
        if any(obs[c].null_count() for c in required) or not set(obs['locus'])<={'TRA','TRB'}:
            raise ValueError('non-null TRA/TRB receptor keys and integer counts required')
        if not obs['junction_aa'].str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').all():raise ValueError('invalid junction')
        totals=obs.group_by('receptor_id').agg(pl.col('count').n_unique().alias('_n'),pl.col('count').first())
        if totals['_n'].ne(1).any():raise ValueError('all chains of a receptor must have the same count')
        total_count=sum(totals['count'].to_list())
        if total_count>2**63-1:raise ValueError('total count exceeds int64')
        parts=[];audits=[]
        for (locus,),part in obs.partition_by('locus',as_dict=True).items():
            result=self._filter(part,locus);audits.append(result)
            parts.append(part.join(result.filter(pl.col('retained')).select('junction_aa'),on='junction_aa',how='semi'))
        retained=pl.concat(parts).with_row_index('query_id')
        hit_parts=[]
        for (locus,),part in retained.partition_by('locus',as_dict=True).items():
            if locus not in self.annotator.loci:continue
            hits=self.annotator.hits(part,locus=locus,scope=self.scope,match_v=self.match_v,
                match_j=False,exclude_exact=False,threads=self.threads)
            if hits.height:
                hit_parts.append(hits.join(part.select('query_id','receptor_id'),on='query_id',validate='m:1').select(
                    'receptor_id',*PMHC,pl.lit(locus).alias('locus'),pl.col('query_cdr3').alias('junction_aa'),
                    'db_cdr3','score').unique())
        schema={'receptor_id':obs.schema['receptor_id'],**{c:self.inventory.schema[c] for c in PMHC},
                'locus':pl.String,'junction_aa':pl.String,'db_cdr3':pl.String,'score':pl.Int64}
        hits=pl.concat(hit_parts) if hit_parts else pl.DataFrame(schema=schema)
        counts=hits.select('receptor_id',*PMHC).unique().join(totals.select('receptor_id','count'),on='receptor_id',validate='m:1').group_by(PMHC).agg(
            pl.len().alias('matched_receptors'),pl.col('count').sum().alias('matched_count'))
        summary=self.reference_counts.join(counts,on=PMHC,how='left',nulls_equal=True).with_columns(
            pl.col('matched_receptors','matched_count').fill_null(0),pl.lit(totals.height).alias('n_receptors'),
            pl.lit(total_count).alias('total_count'),pl.lit(retained['receptor_id'].n_unique()).alias('retained_receptors'))
        return summary,hits,pl.concat(audits)
