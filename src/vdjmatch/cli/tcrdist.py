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
    raw=_read_table(a.sample)
    if 'query_id' not in raw.columns:raw=raw.with_row_index('query_id')
    if raw['query_id'].null_count() or raw['query_id'].n_unique()!=raw.height:
        raise ValueError('query_id must be non-null and unique')
    organism=_organism(a.species)
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
                    reference_part['_allele_alpha'].to_list(),reference_part['_allele_beta'].to_list(),species=organism,threads=a.threads)
            else:
                seqs=part['cdr3'].to_list()
                distances=distance_matrix(seqs,rs,part['_allele'].to_list(),rv,species=organism,threads=a.threads)
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
        'distance':'TCRdist3 default3*CDR3+CDR1+CDR2+CDR2.5; '+organism+' combo_xcr_2024-03-05',
        'calibration':'none','wall_seconds':time.perf_counter()-start},indent=2)+'\n')
    return 0
