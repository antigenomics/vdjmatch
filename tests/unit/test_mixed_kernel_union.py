"""Mixed compact union batches against independent fresh single/pair services."""
from itertools import product

import polars as pl
import pytest
from polars.testing import assert_frame_equal

from vdjmatch.evalue.kernel import KERNEL_SCORERS
from vdjmatch.runner.annotation import annotate_repertoires


def inputs(tmp_path):
    meta=dict(species='HomoSapiens',epitope='A',mhc_a='HLA-A*02:01',mhc_b='B2M',mhc_class='MHCI')
    suffixes=[''.join(s) for s in list(product('ACDEFGHIKLMNPQRSTVWY',repeat=2))[:70]]
    chain={'TRA':['CAVRDS'+s+'QLIW' for s in suffixes], 'TRB':['CASSL'+s+'AYEQYF' for s in suffixes]}
    singles={};paired=[]
    for locus in ['TRA','TRB']:
        rows=[]
        for i,seq in enumerate(chain[locus]):
            v='TRAV12-1' if locus=='TRA' else 'TRBV19'
            rows.append(dict(query_id='q'+str(i),sequence_id='seq'+str(i),clone_id='single'+str(i),
                cell_id='cell'+str(i),junction_aa=seq,v_call=v,j_call=locus+'J1',locus=locus,**meta))
            # Raw V source quirks remain independent; same AA/J has a different V.
            pv='TRAV12-2' if locus=='TRA' and i%3==0 else 'TRBV7-9' if locus=='TRB' and i%3==0 else v
            paired.append(dict(query_id='p'+str(i)+locus,sequence_id='pair-seq'+str(i)+locus,clone_id='pair'+str(i),
                cell_id='pair-cell'+str(i),junction_aa=seq,v_call=pv,j_call=locus+'J2',locus=locus,**meta))
        rows += [rows[0]|{'query_id':'different-J','sequence_id':'different-J','clone_id':'different-J','cell_id':'different-J',
                         'j_call':locus+'J2'},
                 rows[1]|{'query_id':'invalid','sequence_id':'invalid','clone_id':'invalid','cell_id':'invalid','v_call':'UNKNOWN'},
                 rows[2]|{'query_id':'nohit','sequence_id':'nohit','clone_id':'nohit','cell_id':'nohit',
                          'junction_aa':'C'+'W'*90+'F'}]
        path=tmp_path/(locus+'.tsv');pl.DataFrame(rows).write_csv(path,separator='\t');singles[locus]=path
    paired += [paired[0]|{'query_id':'partial-a','sequence_id':'partial-a','clone_id':'partial','cell_id':'partial','v_call':'UNKNOWN'},
               paired[70]|{'query_id':'partial-b','sequence_id':'partial-b','clone_id':'partial','cell_id':'partial','v_call':'TRBV19'}]
    pairpath=tmp_path/'paired.tsv';pl.DataFrame(paired).write_csv(pairpath,separator='\t')
    refs=[]
    for i in [0,1,3,5]:
        for e,mhc in [('A','HLA-A*02:01'),('B','HLA-A*02:01'),('A','HLA-A*03:01')]:
            for locus in ['TRA','TRB']:
                refs.append(dict(gene=locus,cdr3=chain[locus][i],v='TRAV12-1' if locus=='TRA' else 'TRBV19',
                    complex_id=f'{i}-{e}-{mhc}',**(meta|{'epitope':e,'mhc_a':mhc})))
    ref=tmp_path/'reference.tsv';pl.DataFrame(refs).write_csv(ref,separator='\t')
    return [singles['TRA'],singles['TRB'],pairpath],ref


def comparable(frame):
    # Service/global numeric IDs differ between independent invocations; restored source IDs do not.
    frame=frame.drop('query_id','query_ida','query_idb',strict=False)
    keys=[c for c in ['sample_id','source_query_id','source_query_ida','source_query_idb','source_pair_id'] if c in frame.columns]
    return frame.sort(keys).select(sorted(frame.columns))


@pytest.mark.parametrize('scorer',KERNEL_SCORERS)
@pytest.mark.parametrize('exclude',[False,True])
def test_union_equals_separate_fresh_services_and_searches_three_populations(tmp_path,monkeypatch,scorer,exclude):
    paths,reference=inputs(tmp_path)
    options=dict(reference=reference,kernel_scorer=scorer,exclude_exact=exclude,threads=2)
    pieces={}
    for path in paths:
        out=annotate_repertoires(path,output_prefix=tmp_path/('separate-'+path.stem),**options)
        for kind in out['output_paths']:pieces.setdefault(kind,[]).append(out[kind])
    from vdjmatch.match import tcrdist
    original=tcrdist.grouped_total_distance_count_batch;calls=[]
    def count(q,r,cutoffs,**kwargs):
        calls.append((kwargs['paired'],q.height))
        return original(q,r,cutoffs,**kwargs)
    monkeypatch.setattr(tcrdist,'grouped_total_distance_count_batch',count)
    out=annotate_repertoires(paths,output_prefix=tmp_path/'union',**options)
    assert [paired for paired,_ in calls]==[False,False,True]
    assert calls[0][1]>64 and calls[1][1]>64
    assert set(out['output_paths'])==set(pieces)
    for kind,frames in pieces.items():
        assert_frame_equal(comparable(out[kind]),comparable(pl.concat(frames,how='diagonal_relaxed')),
            check_exact=True,check_row_order=True)
    assert out['queries']['query_id'].n_unique()==out['queries'].height
    assert out['manifest']['mixed_chain_union']
    assert out['manifest']['population_service_calls']==3
    assert [x['query_chain_rows'] for x in out['manifest']['populations']]==[73,73,142]
    assert [x['service_chain_rows'] for x in out['manifest']['service_batches']]==[215,215,142]
    nohit=out['selection'].filter(pl.col('source_query_id')=='nohit')
    assert nohit['kernel_score'].null_count()==2 and nohit['status'].unique().to_list()==['no_neighbours']
    partial=out['selection'].filter((pl.col('locus')=='paired') & (pl.col('clone_id')=='partial'))
    assert partial.height==1 and not partial['available'][0]
    assert partial['has_support'][0]


@pytest.mark.parametrize('edge',['assayed-missing-target','empty-reference','unlinked-reference','all-invalid','missing-single-TRA'])
def test_union_edge_dispositions_equal_separate_services(tmp_path,edge):
    paths,reference=inputs(tmp_path)
    options=dict(reference=reference,exclude_exact=True,threads=2)
    if edge=='assayed-missing-target':
        options['targets_from_sample']=True
        for path in paths:
            frame=pl.read_csv(path,separator='\t')
            frame.with_columns(pl.when(pl.col('junction_aa').is_in(['CAVRDSAAQLIW','CASSLAAAYEQYF']))
                .then(None).otherwise(pl.col('epitope')).alias('epitope')).write_csv(path,separator='\t')
    elif edge=='empty-reference':
        pl.read_csv(reference,separator='\t').head(0).write_csv(reference,separator='\t')
    elif edge=='unlinked-reference':
        pl.read_csv(reference,separator='\t').with_columns(pl.lit('0').alias('complex_id')).write_csv(reference,separator='\t')
    elif edge=='all-invalid':
        for path in paths:
            pl.read_csv(path,separator='\t').with_columns(pl.lit('UNKNOWN').alias('v_call')).write_csv(path,separator='\t')
    else:
        paths=paths[1:]
    pieces={}
    for path in paths:
        out=annotate_repertoires(path,output_prefix=tmp_path/('separate-'+path.stem),**options)
        for kind in out['output_paths']:pieces.setdefault(kind,[]).append(out[kind])
    out=annotate_repertoires(paths,output_prefix=tmp_path/'union',**options)
    for kind,frames in pieces.items():
        assert_frame_equal(comparable(out[kind]),comparable(pl.concat(frames,how='diagonal_relaxed')),check_exact=True)
    if edge=='unlinked-reference':
        paired=out['selection'].filter(pl.col('locus')=='paired')
        assert paired['has_support'].any()
        assert paired['linked_selection_status'].unique().to_list()==['no_support']
    assert out['manifest']['mixed_chain_union']


@pytest.mark.parametrize('weighting',['kernel','radius'])
def test_mixed_detailed_dispatch_remains_separate(tmp_path,monkeypatch,weighting):
    import json
    from pathlib import Path
    from vdjmatch.cli import tcrdist
    paths,reference=inputs(tmp_path);calls=[]
    def main(a, *, v_loop_distances):
        assert v_loop_distances['matrix']=='blosum62'
        raw=pl.read_csv(a.sample,separator='\t');calls.append((a.locus,raw.height))
        q=(raw.select('clone_id').unique(maintain_order=True).with_row_index('query_id')
           if a.locus=='paired' else raw.select('query_id'))
        q.with_columns(pl.lit('no_neighbours').alias('status')).write_csv(a.output_prefix+'.queries.tsv',separator='\t')
        q.with_columns(pl.lit('no_support').alias('selection_status')).write_csv(a.output_prefix+'.selection.tsv',separator='\t')
        Path(a.output_prefix+'.manifest.json').write_text(json.dumps({'locus':a.locus}))
    monkeypatch.setattr(tcrdist,'main',main)
    out=annotate_repertoires(paths,reference=reference,output_prefix=tmp_path/'out',threads=2,
        selection_only=False,neighbour_weighting=weighting)
    assert calls==[('TRA',73),('TRB',73),('paired',142)]
    assert 'mixed_chain_union' not in out['manifest']
