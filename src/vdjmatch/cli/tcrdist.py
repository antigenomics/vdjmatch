"""Native TCRdist3-compatible single-chain/linked-pair distance neighbourhoods."""
from pathlib import Path


def register(subparsers):
    p=subparsers.add_parser('tcrdist-neighbours',help='native TCRdist3-compatible distance balls')
    p.add_argument('sample')
    p.add_argument('--vdjdb',required=True)
    p.add_argument('--locus',choices=['TRA','TRB','paired'],required=True)
    p.add_argument('--species',default='HomoSapiens')
    p.add_argument('--epitope',action='append')
    p.add_argument('--mhc-a')
    p.add_argument('--mhc-b')
    p.add_argument('--mhc-match',choices=['exact','compatible'],default='exact')
    p.add_argument('--targets-from-sample',action='store_true',help='score only each row\'s declared assayed pMHC with compatible MHC resolution')
    p.add_argument('--threads',type=int,default=1)
    p.add_argument('--junction-ends',choices=['tcrdist','trim3'],default='tcrdist',
                   help='TCRdist3 trims 3 N/2 C residues; trim3 tests 3 residues at both ends')
    p.add_argument('--background',choices=['real','generative','vdjdb-other'],
                   help='experimental assayed-target ball enrichment; no P-value')
    p.add_argument('--control',help='raw junction/V control table for real or generative background')
    p.add_argument('--position-weighting',choices=['uniform','significance'],default='uniform',
                   help='enrichment geometry; significance uses the complete shipped flank/core decay profile')
    p.add_argument('--neighbour-weighting',choices=['ball','linear','local-rank'],default='ball',
                   help='ball/linear enrichment or experimental background-local rank evidence')
    p.add_argument('--rank-signal-prior',type=float,default=.5,
                   help='prior for the local-rank model only; not biological specificity prevalence')
    p.add_argument('--radius',type=int,default=90,help='maximum distance; paired uses the alpha+beta sum against one linked reference')
    p.add_argument('--exclude-exact',action='store_true')
    p.add_argument('--unbounded-nearest',action='store_true',help='exhaustive nearest reference for each sample-declared target; neighbour counts retain the radius')
    p.add_argument('--output-prefix',required=True)
    p.set_defaults(func=main)


def _paired_input(path, raw, model, organism):
    """Preserve raw pair/chain identities, including rows dropped during ingestion."""
    import polars as pl
    from .. import io
    from ..match.tcrdist import resolve_v_alleles
    from ..aggregate.candidates import PMHC
    if 'clone_id' not in raw.columns:
        raise ValueError('paired AIRR input requires clone_id linkage')
    locus=pl.col('locus').str.to_uppercase() if 'locus' in raw.columns else pl.col('v_call').str.slice(0,3).str.to_uppercase()
    raw=raw.with_columns(pl.col('clone_id').cast(pl.String).alias('pair_id'),locus.alias('_locus'))
    if raw.filter(pl.col('pair_id').is_null() | (pl.col('pair_id')=='')).height:
        raise ValueError('every paired row requires nonempty clone_id linkage')
    if raw.filter(~pl.col('_locus').is_in(['TRA','TRB']).fill_null(False)).height:
        raise ValueError('paired rows require explicit TRA/TRB roles when V calls do not establish locus')
    if raw.group_by('pair_id','_locus').len().filter(pl.col('len')>1).height:
        raise ValueError('ambiguous paired input: multiple rows for one clone/locus')
    metadata=[c for c in ['species',*PMHC] if c in raw.columns]
    if 'species' in metadata:
        raw=raw.with_columns(pl.col('species').str.strip_chars().str.to_lowercase().replace(
            {'human':'human','homosapiens':'human','homo sapiens':'human',
             'mouse':'mouse','musmusculus':'mouse','mus musculus':'mouse'}).alias('_species'))
    shared=['_species' if c=='species' else c for c in metadata]
    if shared and raw.group_by('pair_id').agg(*(pl.col(c).n_unique().alias(c) for c in shared)).filter(
            pl.any_horizontal(*(pl.col(c)>1 for c in shared))).height:
        raise ValueError('paired query has conflicting species or peptide/MHC annotations')
    v_column=io.columns._resolve(raw).get('v')
    pairs=raw.group_by('pair_id',maintain_order=True).agg(
        pl.col('clone_id').first(),*(pl.col(c).first() for c in metadata),
        *[pl.col(c).filter(pl.col('_locus')==locus).first().alias(c+suffix)
          for locus,suffix in [('TRA','a'),('TRB','b')]
          for c in ['query_id',*(['sequence_id'] if 'sequence_id' in raw.columns else [])]],
        *[(pl.col(v_column).filter(pl.col('_locus')==locus).first() if v_column else
           pl.lit(None,dtype=pl.String)).alias('v_call_'+name)
          for locus,name in [('TRA','alpha'),('TRB','beta')]])
    pairs=pairs.with_row_index('query_id')
    cells,ingestion=io.read_cell(path,link='clone_id',source='airr',valid_aa=False,return_report=True)
    # Chain IDs and V calls come from raw rows, including missing/invalid sequences.
    cells=cells.drop([c for c in cells.columns if c.startswith('query_id') or c.startswith('sequence_id')])
    pairs=pairs.join(cells,on='pair_id',how='left',validate='1:1',maintain_order='left')
    for suffix,locus,name in [('a','TRA','alpha'),('b','TRB','beta')]:
        allele='_allele'+suffix
        pairs=pairs.with_columns(resolve_v_alleles(pairs['v_call_'+name],model,locus=locus).alias(allele))
        allowed=[v for v,m in model.items() if m['locus']==locus]
        pairs=pairs.with_columns(
            pl.col('cdr3'+suffix).str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]{8,}$').fill_null(False).alias('_sequence_valid'+suffix),
            pl.col(allele).is_in(allowed).fill_null(False).alias('_gene_valid'+suffix))
        pairs=pairs.with_columns((pl.col('_sequence_valid'+suffix)&pl.col('_gene_valid'+suffix)).alias('_valid'+suffix))
    species_ok=(pl.col('species').str.strip_chars().str.to_lowercase().is_in(
        ['human','homosapiens','homo sapiens'] if organism=='human' else ['mouse','musmusculus','mus musculus'])
        if 'species' in pairs.columns else pl.lit(True))
    pairs=pairs.with_columns(
        (pl.col('_valida') & pl.col('_validb') & species_ok.fill_null(False)).alias('available'),
        pl.when(~species_ok.fill_null(False)).then(pl.lit('unselected_species'))
        .when(pl.col('query_ida').is_null()).then(pl.lit('missing_alpha'))
        .when(pl.col('query_idb').is_null()).then(pl.lit('missing_beta'))
        .when(~pl.col('_sequence_valida')).then(pl.lit('invalid_alpha_query'))
        .when(~pl.col('_sequence_validb')).then(pl.lit('invalid_beta_query'))
        .when(~pl.col('_gene_valida')).then(pl.lit('unknown_alpha_gene'))
        .when(~pl.col('_gene_validb')).then(pl.lit('unknown_beta_gene'))
        .otherwise(pl.lit('searched')).alias('status'))
    return pairs,ingestion


def main(a):
    from importlib.metadata import version
    import json
    import time
    import numpy as np
    import polars as pl
    from .. import db,io
    from ..api import _prepare
    from ..aggregate.candidates import PMHC
    from ..io.airr import _read_table
    from ..db.cache import sha256
    from seqtree import _core
    from ..match.tcrdist import distance_matrix,paired_distance_matrix,load_v_loops,resolve_v_alleles
    from ..evalue.control import _organism
    start=time.perf_counter()
    if a.threads<1 or a.radius<0:raise ValueError('positive threads and nonnegative radius required')
    if a.unbounded_nearest and not a.targets_from_sample:
        raise ValueError('--unbounded-nearest requires --targets-from-sample')
    if a.background and not a.targets_from_sample:
        raise ValueError('--background requires --targets-from-sample')
    if a.locus=='paired' and a.background and (a.background!='vdjdb-other' or a.position_weighting!='uniform'):
        raise ValueError('paired background requires vdjdb-other and uniform positional weights; generative/real joint controls are not defined')
    if a.neighbour_weighting=='linear' and (not a.background or a.radius<=0):
        raise ValueError('linear neighbour weighting requires a background and positive radius')
    if a.neighbour_weighting=='local-rank' and (not a.background or not 0<a.rank_signal_prior<1):
        raise ValueError('local-rank requires a background and prior strictly between zero and one')
    if bool(a.control) != (a.background in {'real','generative'}):
        raise ValueError('--control is required only for --background real or generative')
    if a.position_weighting!='uniform' and (not a.background or a.junction_ends!='tcrdist'):
        raise ValueError('significance weighting requires --background and excludes an additional trim3')
    raw=_read_table(a.sample)
    if 'query_id' not in raw.columns:raw=raw.with_row_index('query_id')
    if raw['query_id'].null_count() or raw['query_id'].n_unique()!=raw.height:
        raise ValueError('query_id must be non-null and unique')
    organism=_organism(a.species)
    ctrim=2 if a.junction_ends=='tcrdist' else 3
    paired=a.locus=='paired'
    model=load_v_loops(organism)
    if paired:
        pair_rows,ingestion=_paired_input(a.sample,raw,model,organism)
        q=pair_rows.filter(pl.col('available'))
        metadata=pair_rows.select('query_id',*[c for c in ['species',*PMHC] if c in pair_rows.columns])
    else:
        data,ingestion=io.read_rearrangement(a.sample,source='airr',valid_aa=False,return_report=True)
        _,q=_prepare(data)
        metadata=raw
    reference_species={'human':'HomoSapiens','mouse':'MusMusculus'}[organism]
    r=db.load(a.vdjdb,species=reference_species,gene=None if paired else a.locus,epitope=a.epitope,
              mhc_a=a.mhc_a,mhc_b=a.mhc_b,mhc_match=a.mhc_match)
    original=r.height
    allowed=[v for v,m in model.items() if m['locus']==a.locus]
    def prepare(frame):
        return frame.with_columns(resolve_v_alleles(frame['v'],model,locus=a.locus).alias('_allele')).filter(
                pl.col('_allele').is_in(allowed)&(pl.col('cdr3').str.len_chars()>=8))
    if paired:
        from ..evalue.paired import build_paired_ref
        r=build_paired_ref(r)
        paired_observations=r.height
        for suffix,locus in [('alpha','TRA'),('beta','TRB')]:
            allowed_chain=[v for v,m in model.items() if m['locus']==locus]
            r=r.with_columns(resolve_v_alleles(r['v_'+suffix],model,locus=locus).alias('_allele_'+suffix)).filter(
                    pl.col('_allele_'+suffix).is_in(allowed_chain) & (pl.col(suffix).str.len_chars()>=8))
        search_keys=['alpha','beta','_allele_alpha','_allele_beta']
        r=r.unique(subset=[*search_keys,*PMHC],maintain_order=True)
    else:
        q=prepare(q.filter(pl.col('locus')==a.locus))
        r=prepare(r.filter(pl.col('reference_valid'))).unique(subset=['cdr3','_allele',*PMHC],maintain_order=True)
        search_keys=['cdr3','_allele']
    if a.targets_from_sample:
        if a.epitope or a.mhc_a or a.mhc_b:
            raise ValueError('--targets-from-sample cannot be combined with explicit target selectors')
        from .historical import assayed_targets
        groups=assayed_targets(q.drop([c for c in ['species',*PMHC] if c in q.columns]),metadata,r,a.species)
    else:
        groups=[(None,q,r)]
    # Counts share one distance on references and controls; no dense control matrix.
    enrichment_rows=[]
    background_provenance=None
    background_costs={'control_load':0.0,'native_count_batches':0.0}
    if a.background:
        background_start=time.perf_counter()
        groups=list(groups)  # assayed_targets is streamed; the new arm uses it twice
        from ..match.tcrdist import total_distance_count_batch,paired_total_distance_mass_batch
        from ..evalue.control import raw_gene_background
        linear=a.neighbour_weighting=='linear'
        local=a.neighbour_weighting=='local-rank'
        query_keys=['cdr3a','cdr3b','_allelea','_alleleb'] if paired else ['cdr3','_allele']
        options=dict(species=organism,threads=a.threads,ctrim=ctrim,
                     position_weighting=a.position_weighting,exclude_exact=a.exclude_exact,linear_mass=linear)
        def counts(part,population):
            if not part.height:
                return [],[],np.zeros((0,a.radius+1),dtype=np.uint64) if local else None
            if not population.height:
                return [0]*part.height,[0.0]*part.height,np.zeros((part.height,a.radius+1),dtype=np.uint64) if local else None
            count_start=time.perf_counter()
            thresholds=[list(range(0,200*a.radius+1,200))]*part.height if local else [[200*a.radius]]*part.height
            if paired:
                native=paired_total_distance_mass_batch(
                    part['cdr3a'].to_list(),part['cdr3b'].to_list(),
                    population['alpha'].to_list(),population['beta'].to_list(),
                    part['_allelea'].to_list(),part['_alleleb'].to_list(),
                    population['_allele_alpha'].to_list(),population['_allele_beta'].to_list(),
                    thresholds,species=organism,threads=a.threads,
                    ctrim=ctrim,exclude_exact=a.exclude_exact)
                curve=np.asarray(native[0],dtype=np.uint64) if local else None
                result=[row[-1] for row in native[0]]
                mass=[row[0]/(200*a.radius) for row in native[1]] if linear else result
            else:
                native=total_distance_count_batch(part['cdr3'].to_list(),population['cdr3'].to_list(),
                    part['_allele'].to_list(),population['_allele'].to_list(),
                    thresholds,**options)
                curve=np.asarray(native,dtype=np.uint64) if local else None
                result=[row[-1] for row in (native[0] if linear else native)]
                mass=[row[0]/(200*a.radius) for row in native[1]] if linear else result
            background_costs['native_count_batches']+=time.perf_counter()-count_start
            return result,mass,curve
        def populations(part,population):
            if not a.exclude_exact:
                return [population.height]*part.height
            if paired:
                exact=population.select(pl.col('alpha').alias('cdr3a'),pl.col('beta').alias('cdr3b')).join(
                    part.select('cdr3a','cdr3b').unique(),on=['cdr3a','cdr3b'],how='semi')
                multiplicity=dict((tuple(row[:2]),row[2]) for row in exact.group_by('cdr3a','cdr3b').len().iter_rows())
                return [population.height-multiplicity.get(key,0) for key in part.select('cdr3a','cdr3b').iter_rows()]
            # Only query identities can be punctured; avoid a full-control Python map.
            exact=population.filter(pl.col('cdr3').is_in(part['cdr3'].unique().implode()))
            multiplicity=dict(exact.group_by('cdr3').len().iter_rows())
            return [population.height-multiplicity.get(s,0) for s in part['cdr3']]
        shared=None
        if a.background!='vdjdb-other':
            load_start=time.perf_counter()
            controls,background_provenance=raw_gene_background(a.control,a.locus,organism,deduplicate=a.background!='generative')
            background_costs['control_load']=time.perf_counter()-load_start
            # One native control batch, shared transiently across declared targets.
            parts=[part.select('cdr3','_allele') for _,part,_ in groups if part.height]
            usable=pl.concat(parts).unique(maintain_order=True) if parts else q.select('cdr3','_allele').head(0)
            m,w,mc=counts(usable,controls)
            shared=usable.with_columns(pl.Series('_m',m,dtype=pl.UInt64),pl.Series('control_weight',w,dtype=pl.Float64),
                                       pl.Series('_M',populations(usable,controls),dtype=pl.UInt64))
            if local:
                shared=shared.with_columns(pl.Series('_rank_control',mc))
            del mc
        else:
            background_provenance={'kind':'vdjdb_other_epitopes','reference':db.provenance(a.vdjdb),
                'target_removal':'all records of target peptide, before key deduplication',
                'shared_key_policy':'other epitope annotations retained; no arbitrary label selection',
                'observation_unit':('distinct linked full alpha/beta junctions + resolved V alleles' if paired else
                                    'distinct full junction + resolved V allele')+'; no abundance weights',
                'restriction':'same species/locus; other epitopes may have other MHC restrictions'}
        for task,part,reference_part in groups:
            unique=part.select(query_keys).unique(maintain_order=True)
            if not unique.height:
                continue
            reference_keys=reference_part.unique(search_keys,maintain_order=True)
            n,reference_weight,nc=counts(unique,reference_keys)
            N=populations(unique,reference_keys)
            if shared is not None:
                evidence=unique.join(shared,on=query_keys,how='left',validate='1:1',maintain_order='left')
            else:
                other=r.filter(pl.col('epitope')!=task['epitope']).unique(search_keys,maintain_order=True)
                m,w,mc=counts(unique,other)
                evidence=unique.with_columns(pl.Series('_m',m,dtype=pl.UInt64),pl.Series('control_weight',w,dtype=pl.Float64),
                                            pl.Series('_M',populations(unique,other),dtype=pl.UInt64))
                if local:
                    evidence=evidence.with_columns(pl.Series('_rank_control',mc))
                del mc
            evidence=evidence.with_columns(pl.Series('n_reference',n,dtype=pl.UInt64),
                                          pl.Series('reference_weight',reference_weight,dtype=pl.Float64),
                                          pl.Series('reference_population',N,dtype=pl.UInt64))
            evidence=evidence.rename({'_m':'n_control','_M':'control_population'})
            available=(pl.col('reference_population')>0)&(pl.col('control_population')>0)
            evidence=evidence.with_columns(
                pl.when(available).then(pl.col('reference_population')*(pl.col('control_weight')+.5)/(pl.col('control_population')+1)).alias('expected_count'),
                pl.when(available).then(pl.lit('available')).otherwise(pl.lit('empty_population')).alias('enrichment_status'))
            evidence=evidence.with_columns((pl.col('reference_weight')/pl.col('expected_count')).alias('enrichment'),
                (pl.col('n_control')==0).alias('zero_control_hits'),(pl.col('control_weight')==0).alias('zero_control_weight'))
            if local:
                from ..evalue.local_rank import local_rank_evidence
                rank=local_rank_evidence(nc,evidence['_rank_control'].to_numpy(),N,
                                         evidence['control_population'].to_numpy(),prior=a.rank_signal_prior)
                evidence=evidence.drop('_rank_control').with_columns(
                    *(pl.Series(k,v) for k,v in rank.items()),
                    pl.lit(a.rank_signal_prior).alias('rank_signal_prior'),
                    pl.when(available).then(pl.lit('conditional_rank_model_only')).otherwise(pl.lit('empty_population')).alias('rank_model_status'))
                evidence=evidence.with_columns(
                    pl.when(available).then(pl.col('rank_bayes_factor')).alias('rank_bayes_factor'),
                    *(pl.when(available).then(pl.col(c)).alias(c) for c in ['p_rank_bound','p_nearest_rank']))
            del nc
            enrichment_rows.append(part.select('query_id',*query_keys).join(evidence,
                on=query_keys,how='left',validate='m:1',maintain_order='left').drop(query_keys)
                .with_columns(*(pl.lit(task[c],dtype=pl.String).alias(c) for c in PMHC)))
        background_costs['total']=time.perf_counter()-background_start
    rows=[]
    dispositions={}
    for task,query_part,reference_part in groups:
        if task is not None:
            reference_part=reference_part.unique(search_keys,maintain_order=True)
            for query_id in query_part['query_id']:
                dispositions[query_id]=('invalid_target' if any(not isinstance(v,str) or not v.strip() for v in task.values())
                                        else 'searched' if reference_part.height else 'no_reference')
        if not reference_part.height or not query_part.height:
            continue
        # Paired calls hold two matrices plus native/gather temporaries; bound them together.
        matrix_budget=32*1024**2 if paired else 64*1024**2
        bytes_per_pair=16 if paired else 4
        if bytes_per_pair*reference_part.height>matrix_budget:
            raise ValueError('reference exceeds the matrix budget for one query; use a narrower target')
        batch=max(1,matrix_budget//(bytes_per_pair*reference_part.height))
        if not paired:
            rs,rv=reference_part['cdr3'].to_list(),reference_part['_allele'].to_list()
        for part in query_part.iter_slices(batch):
            if paired:
                qa,qb=part['cdr3a'].to_list(),part['cdr3b'].to_list()
                ra,rb=reference_part['alpha'].to_list(),reference_part['beta'].to_list()
                distances=paired_distance_matrix(qa,qb,ra,rb,part['_allelea'].to_list(),part['_alleleb'].to_list(),
                    reference_part['_allele_alpha'].to_list(),reference_part['_allele_beta'].to_list(),species=organism,threads=a.threads,ctrim=ctrim)
            else:
                seqs=part['cdr3'].to_list()
                distances=distance_matrix(seqs,rs,part['_allele'].to_list(),rv,species=organism,threads=a.threads,ctrim=ctrim)
            keep=distances<=a.radius
            nonexact=True
            if a.exclude_exact:
                if paired:
                    nonexact=~((np.asarray(qa)[:,None]==np.asarray(ra)[None,:]) &
                               (np.asarray(qb)[:,None]==np.asarray(rb)[None,:]))
                else:
                    nonexact=np.asarray(seqs)[:,None]!=np.asarray(rs)[None,:]
                keep&=nonexact
            if a.unbounded_nearest:
                sentinel=np.iinfo(np.int32).max
                nearest=np.min(distances,axis=1,where=nonexact,initial=sentinel)
                available=nearest!=sentinel
                h=part.filter(pl.Series(available)).select('query_id').with_columns(
                    pl.Series('distance',nearest[available]),
                    pl.Series('n_neighbours',keep.sum(axis=1,dtype=np.uint32)[available]))
                rows.append(h.with_columns(*(pl.lit(task[c],dtype=pl.String).alias(c) for c in PMHC)))
                continue
            qi,ri=np.nonzero(keep)
            h=pl.DataFrame({'_qi':qi,'_ri':ri,'distance':distances[qi,ri]})
            h=h.join(part.with_row_index('_qi').select('_qi','query_id'),on='_qi')
            if task is None:
                h=h.join(reference_part.with_row_index('_ri').select('_ri',*PMHC),on='_ri')
            else:
                h=h.with_columns(*(pl.lit(task[c],dtype=pl.String).alias(c) for c in PMHC))
            rows.append(h.group_by('query_id',*PMHC).agg(pl.col('distance').min(),pl.len().alias('n_neighbours')))
    schema={'query_id':q.schema['query_id'],**{c:pl.String for c in PMHC},'distance':pl.Int32,'n_neighbours':pl.UInt32}
    out=pl.concat(rows) if rows else pl.DataFrame(schema=schema)
    out=out.sort('query_id','distance',*PMHC,nulls_last=True)
    if paired:
        status=pair_rows.select(*[c for c in ['query_id','pair_id','clone_id','query_ida','query_idb',
            'sequence_ida','sequence_idb','species',*PMHC] if c in pair_rows.columns],'available','status',
            'v_call_alpha','v_call_beta',
            pl.col('_allelea').alias('resolved_v_allele_alpha'),pl.col('_alleleb').alias('resolved_v_allele_beta'))
        hit_ids=set(out['query_id'])
        status=status.with_columns(pl.Series('status',[
            old if old!='searched' else dispositions.get(query_id) if dispositions.get(query_id) in {'invalid_target','no_reference'} else
            'matched' if query_id in hit_ids else 'no_neighbours'
            for query_id,old in zip(status['query_id'],status['status'])]))
    else:
        status=raw.select([c for c in ('query_id','sequence_id') if c in raw.columns]).join(q.select('query_id').with_columns(pl.lit(True).alias('available')),on='query_id',how='left')
    if not paired:
        v_column=io.columns._resolve(raw).get('v')
        raw_calls=raw[v_column] if v_column else pl.Series([None]*raw.height,dtype=pl.String)
        status=status.with_columns(raw_calls.alias('v_call'),
            resolve_v_alleles(raw_calls,model,locus=a.locus).alias('resolved_v_allele'))
    status=status.with_columns(pl.col('available').fill_null(False))
    if a.targets_from_sample and not paired:
        valid_target_ids=set(dispositions)
        hit_ids=set(out['query_id'])
        species_ids=set(raw.filter(pl.col('species').str.strip_chars().str.to_lowercase().is_in(
            ['human','homosapiens','homo sapiens'] if organism=='human' else ['mouse','musmusculus','mus musculus']))['query_id'])
        status=status.with_columns(pl.Series('status',[
            'unselected_species' if query_id not in species_ids else
            'invalid_query_or_gene' if query_id not in valid_target_ids else
            dispositions[query_id] if dispositions[query_id] in {'no_reference','invalid_target'} else
            'matched' if query_id in hit_ids else 'no_neighbours'
            for query_id in status['query_id']]))
        status=status.with_columns((pl.col('available') & (pl.col('status')!='unselected_species')).alias('available'))
    prefix=Path(a.output_prefix);prefix.parent.mkdir(parents=True,exist_ok=True)
    out.write_csv(str(prefix)+'.candidates.tsv',separator='\t')
    status.write_csv(str(prefix)+'.queries.tsv',separator='\t')
    if a.background:
        schema={'query_id':q.schema['query_id'],**{c:pl.String for c in PMHC},
                **{c:pl.UInt64 for c in ['n_reference','reference_population','n_control','control_population']},
                'reference_weight':pl.Float64,'control_weight':pl.Float64,
                'expected_count':pl.Float64,'enrichment_status':pl.String,'enrichment':pl.Float64,
                'zero_control_hits':pl.Boolean,'zero_control_weight':pl.Boolean}
        enrichment=pl.concat(enrichment_rows) if enrichment_rows else pl.DataFrame(schema=schema)
        if a.neighbour_weighting=='local-rank' and not enrichment_rows:
            enrichment=enrichment.with_columns(
                *(pl.lit(None,dtype=pl.Float64).alias(k) for k in ['rank_bayes_factor','p_rank_bound','p_nearest_rank','posterior_rank_signal','rank_signal_prior']),
                pl.lit(None,dtype=pl.String).alias('rank_model_status'))
        enrichment.sort('query_id',*PMHC).write_csv(str(prefix)+'.enrichment.tsv',separator='\t')
    Path(str(prefix)+'.manifest.json').write_text(json.dumps({'sample_sha256':sha256(Path(a.sample)),
        'software':{'vdjmatch':version('vdjmatch'),'seqtree':version('seqtree'),
                    'seqtree_native_sha256':sha256(Path(_core.__file__)),
                    'scorer_source_sha256':sha256(Path(distance_matrix.__code__.co_filename)),
                    'cli_source_sha256':sha256(Path(__file__))},
        'reference':db.provenance(a.vdjdb),'ingestion':ingestion,'reference_observations':original,
        'reference_usable_keys':r.height,'query_rows':raw.height,'query_usable_rows':2*q.height if paired else q.height,
        **({'query_pairs':pair_rows.height,'query_usable_pairs':q.height,'reference_complete_pairs':paired_observations,
            'query_status_counts':dict(status.group_by('status').len().iter_rows()),
            'reference_key_unit':'distinct linked alpha/beta junction+V-allele+pMHC keys',
            'distance_unit':'alpha+beta sum to the same linked reference pair',
            'exact_exclusion':'both full junctions equal; V calls do not change the puncture',
            'matrix_work_budget_bytes':32*1024**2,
            'independent_marginal_baseline':'not emitted; requires separate full compatible single-chain reference axes'} if paired else {}),
        'parameters':{k:v for k,v in vars(a).items() if k!='func'},
        'selection':('exhaustive nearest nonexcluded reference; n_neighbours counts only within radius'
                     if a.unbounded_nearest else 'radius ball'),
        'distance':('TCRdist3 default' if ctrim==2 else 'experimental symmetric-trim')+
                   ' 3*CDR3+CDR1+CDR2+CDR2.5; '+organism+' combo_xcr_2024-03-05',
        'junction_end_trim':{'n_terminal':3,'c_terminal':ctrim},
        **({'background_enrichment':{'background':background_provenance,'stage_wall_seconds':background_costs,
            'geometry':'full-profile decay' if a.position_weighting=='significance' else f'trim3/{ctrim}',
            'formula':'reference_weight/[reference_population*(control_weight+0.5)/(control_population+1)]',
            'kernel':'max(0,1-distance/radius)' if a.neighbour_weighting=='linear' else 'distance<=radius',
            'denominators':'reference distinct receptor/V-allele keys; controls use declared population units; both after optional full-junction exclusion',
            'calibration':('none; posterior mean kernel mass under unit Dirichlet prior split equally at distance0/outside-radius; ranking only'
                           if a.neighbour_weighting=='linear' else 'none; Jeffreys background predictive smoothing, not a Bayes factor or P-value'),
            'nearest_output':'distance and native neighbour counts retain their unweighted comparator geometry'}} if a.background else {}),
        **({'local_rank_model':{
            'score_column':'rank_bayes_factor','distance_bins':'ceil(native distance / 200), fixed integer cutoffs 0..radius',
            'rank_weight':'(M+1)/(H_(M+1)*rank); ties averaged; outside-radius references have conditional null mean',
            'null':'exchangeable reference/control distance labels conditional on pooled distances',
            'alternative':'null label law tilted by the mean normalized reciprocal reference rank',
            'p_rank_bound':'min(1,1/rank_bayes_factor); conservative Markov bound, not exact tail',
            'p_nearest_rank':'hypergeometric first-target rank tail; upper-bin ties conservative; separate statistic',
            'posterior_rank_signal':'prior*BF/(1-prior+prior*BF); model probability, not biological Prob(TP)',
            'prior':a.rank_signal_prior,'sampling_status':'exchangeability is not established by source declaration; generated draws versus distinct references may violate it',
            'selection':'fixed radius, geometry and target; no selected-route or multiple-target calibration'}}
           if a.neighbour_weighting=='local-rank' else {}),
        'calibration':'none','wall_seconds':time.perf_counter()-start},indent=2)+'\n')
    return 0
