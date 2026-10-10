"""Independent tcrdist3==0.3 default-model oracle, generated2026-10-10."""
import numpy as np
import pytest
from vdjmatch.match.tcrdist import distance_matrix, load_v_loops, paired_distance_matrix


def test_native_matches_tcrdist3_germline_and_cdr3_oracle():
    junctions=['CTSKGMKF','CPNKNINGEEPQLCF','CGNKELKNMHITGNARHAGFF','CSTNRTCCVAFAATHEDQYMWMNQTSYRWF']
    genes=['TRBV19*01','TRBV7-9*01','TRBV5-1*01','TRBV20-1*01']
    expected=[[0,153,239,360],[153,0,214,357],[239,214,0,331],[360,357,331,0]]
    a=distance_matrix(junctions,junctions,genes,genes,threads=1)
    b=distance_matrix(junctions,junctions,genes,genes,threads=2)
    assert np.array_equal(a,expected) and np.array_equal(a,b)
    with pytest.raises(ValueError,match='unknown human V'):
        distance_matrix(junctions[:1],junctions[:1],['unknown'],genes[:1])


@pytest.mark.parametrize('species,genes,expected',[
    ('mouse',['TRAV1*01','TRAV2*01','TRAV4D-2*01','TRAV5-2*01'],
     [[0,166,221,358],[166,0,228,339],[221,228,0,319],[358,339,319,0]]),
    ('mouse',['TRBV1*01','TRBV2*01','TRBV3*01','TRBV4*01'],
     [[0,166,238,363],[166,0,223,322],[238,223,0,317],[363,322,317,0]]),
    ('human',['TRAV1-1*01','TRAV1-2*01','TRAV12-1*01','TRAV14/DV4*01'],
     [[0,118,234,353],[118,0,228,349],[234,228,0,322],[353,349,322,0]]),
])
def test_species_templates_match_installed_oracle(species,genes,expected):
    # TCRrep 0.3, combo_xcr_2024-03-05.tsv, cpus=1, deduplicate=False.
    # Mouse alpha includes dot/star templates and its29-position alignment.
    seq=['CTSKGMKF','CPNKNINGEEPQLCF','CGNKELKNMHITGNARHAGFF','CSTNRTCCVAFAATHEDQYMWMNQTSYRWF']
    assert np.array_equal(distance_matrix(seq,seq,genes,genes,species=species),expected)
    # Repeated/out-of-order alleles exercise the small lookup's rectangular gather.
    qi=[3,0,3];ri=[2,0]
    assert np.array_equal(distance_matrix([seq[i] for i in qi],[seq[i] for i in ri],
        [genes[i].split('*')[0] for i in qi],[genes[i] for i in ri],species=species,threads=2),
        np.asarray(expected)[np.ix_(qi,ri)])


def test_paired_distance_uses_linked_reference_pairs():
    a='CASSLGQAYEQYF';b='CASSLGRAYEQYF'
    # Each marginal chain has an exact reference; neither linked pair is exact.
    result=paired_distance_matrix([a],[a],[a,b],[b,a],['TRAV1*01'],['TRBV1*01'],
                                  ['TRAV1*01']*2,['TRBV1*01']*2,species='mouse')
    assert result.tolist()==[[9,9]]
    with pytest.raises(ValueError,match='linked pair'):
        paired_distance_matrix([a],[],[a],[a],['TRAV1*01'],[],['TRAV1*01'],['TRBV1*01'],species='mouse')
    with pytest.raises(ValueError,match='locus'):
        paired_distance_matrix([a],[a],[a],[a],['TRBV1*01'],['TRBV1*01'],
                               ['TRBV1*01'],['TRBV1*01'],species='mouse')
    with pytest.raises(ValueError,match='species'):
        load_v_loops('rat')


def test_tcrdist_cli_preserves_unavailable_queries(tmp_path):
    import polars as pl
    from vdjmatch.cli.__main__ import main
    q=tmp_path/'query.tsv';r=tmp_path/'reference.tsv';prefix=tmp_path/'out'
    pl.DataFrame({'sequence_id':['a','b','c','d','e'],
                  'junction_aa':['CASSLGQAYEQYF',None,'','INVALID!', 'CASSLGQAYEQYF'],
                  'v_call':['TRBV19*01']*4+[None],'locus':['TRB']*5}).write_csv(q,separator='\t')
    pl.DataFrame({'gene':['TRB'],'cdr3':['CASSLGRAYEQYF'],'v':['TRBV19*01'],
                  'species':['HomoSapiens'],'epitope':['E']}).write_csv(r,separator='\t')
    assert main(['tcrdist-neighbours',str(q),'--vdjdb',str(r),'--locus','TRB',
                 '--output-prefix',str(prefix),'--exclude-exact'])==0
    scores=pl.read_csv(str(prefix)+'.candidates.tsv',separator='\t')
    status=pl.read_csv(str(prefix)+'.queries.tsv',separator='\t')
    assert scores['distance'].to_list()==[9]
    assert status['query_id'].to_list()==list(range(5))
    assert status['sequence_id'].to_list()==list('abcde')
    assert status['available'].to_list()==[True,False,False,False,False]


def test_assayed_tcrdist_union_counts_keys_and_preserves_search_disposition(tmp_path):
    import polars as pl
    from vdjmatch.cli.__main__ import main
    q=tmp_path/'q.tsv';r=tmp_path/'r.tsv';prefix=tmp_path/'out'
    pl.DataFrame({'junction_aa':['CASSLGQAYEQYF']*3,'v_call':['TRBV19*01']*3,
        'locus':['TRB']*3,'species':['human','human','mouse'],
        'epitope':['E','absent','E'],'mhc_a':['HLA-A*02']*3,
        'mhc_b':['B2M']*3,'mhc_class':['MHCI']*3}).write_csv(q,separator='\t')
    pl.DataFrame({'cdr3':['CASSLGRAYEQYF']*2,'gene':['TRB']*2,'v':['TRBV19*01']*2,
        'species':['HomoSapiens']*2,'epitope':['E']*2,
        'mhc_a':['HLA-A*02','HLA-A*02:01'],'mhc_b':['B2M']*2,
        'mhc_class':['MHCI']*2}).write_csv(r,separator='\t')
    assert main(['tcrdist-neighbours',str(q),'--vdjdb',str(r),'--locus','TRB',
        '--species','human','--targets-from-sample','--exclude-exact','--output-prefix',str(prefix)])==0
    scores=pl.read_csv(str(prefix)+'.candidates.tsv',separator='\t')
    status=pl.read_csv(str(prefix)+'.queries.tsv',separator='\t')
    assert scores['n_neighbours'].to_list()==[1]
    assert scores['mhc_a'].to_list()==['HLA-A*02']
    assert status['status'].to_list()==['matched','no_reference','unselected_species']


def test_paired_tcrdist_cli_same_reference_and_full_junction_puncture(tmp_path):
    import json
    import polars as pl
    from vdjmatch.cli.__main__ import main
    a='CASSLGQAYEQYF';b='CASSLGRAYEQYF'
    q=tmp_path/'q.tsv';r=tmp_path/'r.tsv';prefix=tmp_path/'out'
    pl.DataFrame({'query_id':['raw-a','raw-b'],'sequence_id':['sa','sb'],'clone_id':['007']*2,
        'locus':['TRA','TRB'],'junction_aa':[a,a],'v_call':['TRAV1-1*01','TRBV19*01']}).write_csv(q,separator='\t')
    # Independent marginal minima are zero, but each linked pair distance is9.
    pl.DataFrame({'complex_id':['1','1','2','2'],'gene':['TRA','TRB']*2,'cdr3':[a,b,b,a],
        'v':['TRAV1-1*01','TRBV19*01']*2,'epitope':['E']*4,'species':['HomoSapiens']*4}).write_csv(r,separator='\t')
    args=['tcrdist-neighbours',str(q),'--vdjdb',str(r),'--locus','paired','--output-prefix',str(prefix)]
    assert main(args+['--exclude-exact'])==0
    scores=pl.read_csv(str(prefix)+'.candidates.tsv',separator='\t')
    status=pl.read_csv(str(prefix)+'.queries.tsv',separator='\t',infer_schema_length=0)
    assert scores['distance'].to_list()==[9]
    assert scores['n_neighbours'].to_list()==[2]
    assert status['clone_id'].to_list()==['007']
    assert status['query_ida'].to_list()==['raw-a'] and status['query_idb'].to_list()==['raw-b']
    assert main(args+['--radius','0'])==0
    assert pl.read_csv(str(prefix)+'.candidates.tsv',separator='\t').height==0
    # Both full junctions matching are excluded even with different V templates.
    pl.DataFrame({'complex_id':['3','3'],'gene':['TRA','TRB'],'cdr3':[a,a],
        'v':['TRAV1-2*01','TRBV7-9*01'],'epitope':['E']*2,'species':['HomoSapiens']*2}).write_csv(r,separator='\t')
    assert main(args+['--exclude-exact'])==0
    assert pl.read_csv(str(prefix)+'.candidates.tsv',separator='\t').height==0
    assert pl.read_csv(str(prefix)+'.queries.tsv',separator='\t')['status'].to_list()==['no_neighbours']
    manifest=json.loads((tmp_path/'out.manifest.json').read_text())
    assert manifest['query_pairs']==1 and manifest['query_rows']==2
    assert manifest['query_status_counts']=={'no_neighbours':1}


def test_paired_tcrdist_cli_unavailable_raw_rows_and_assayed_union(tmp_path,capsys):
    import polars as pl
    from vdjmatch.cli.__main__ import main
    q=tmp_path/'q.tsv';r=tmp_path/'r.tsv';prefix=tmp_path/'out';a='CASSLGQAYEQYF';b='CASSLGRAYEQYF'
    pl.DataFrame({'query_id':['a','b','c','d','e','f','g','h','i'],
        'clone_id':['hit','hit','missing','invalid','invalid','unknown','unknown','empty','empty'],
        'locus':['TRA','TRB','TRA','TRA','TRB','TRA','TRB','TRA','TRB'],
        'junction_aa':[a,a,a,a,'INVALID!',a,a,a,None],
        'v_call':['TRAV1-1*01','TRBV19*01','TRAV1-1*01','TRAV1-1*01','TRBV19*01','unknown','TRBV19*01','TRAV1-1*01','TRBV19*01'],
        'species':['human']*9,'epitope':['E']*9,'mhc_a':['HLA-A*02']*9,
        'mhc_b':['B2M']*9,'mhc_class':['MHCI']*9}).write_csv(q,separator='\t')
    pl.DataFrame({'complex_id':['1','1','2','2'],'gene':['TRA','TRB']*2,'cdr3':[a,b]*2,
        'v':['TRAV1-1*01','TRBV19*01']*2,'epitope':['E']*4,'species':['HomoSapiens']*4,
        'mhc_a':['HLA-A*02']*2+['HLA-A*02:01']*2,'mhc_b':['B2M']*4,'mhc_class':['MHCI']*4}).write_csv(r,separator='\t')
    args=['tcrdist-neighbours',str(q),'--vdjdb',str(r),'--locus','paired',
          '--targets-from-sample','--exclude-exact','--output-prefix',str(prefix)]
    assert main(args)==0
    scores=pl.read_csv(str(prefix)+'.candidates.tsv',separator='\t')
    status=pl.read_csv(str(prefix)+'.queries.tsv',separator='\t')
    assert scores['distance'].to_list()==[9] and scores['n_neighbours'].to_list()==[1]
    assert scores['mhc_a'].to_list()==['HLA-A*02']
    assert status['status'].to_list()==['matched','missing_beta','invalid_beta_query','unknown_alpha_gene','invalid_beta_query']
    assert status['query_ida'].to_list()==['a','c','d','f','h']
    assert status['query_idb'].to_list()==['b',None,'e','g','i']
    raw=pl.read_csv(q,separator='\t').with_columns(pl.when(pl.col('query_id')=='b')
        .then(pl.lit('mouse')).otherwise(pl.col('species')).alias('species'))
    raw.write_csv(q,separator='\t')
    with pytest.raises(SystemExit) as error:
        main(args)
    assert error.value.code==2 and 'conflicting species' in capsys.readouterr().err


@pytest.mark.parametrize('locus', ['TRB', 'paired'])
def test_tcrdist_unbounded_nearest_outside_ball_and_exact_only(tmp_path, capsys, monkeypatch, locus):
    import polars as pl
    from vdjmatch.cli.__main__ import main
    from vdjmatch.match import tcrdist as native
    function='paired_distance_matrix' if locus=='paired' else 'distance_matrix'
    original=getattr(native,function)
    def read_only(*args, **kwargs):
        matrix=original(*args, **kwargs)
        matrix.setflags(write=False)
        return matrix
    monkeypatch.setattr(native,function,read_only)
    a='CASSLGQAYEQYF';b='CASSLGRAYEQYF'
    sample=tmp_path/'q.tsv';reference=tmp_path/'r.tsv';prefix=tmp_path/'out'
    if locus=='paired':
        raw=pl.DataFrame({'query_id':['raw-a','raw-b'],'clone_id':['007']*2,
            'locus':['TRA','TRB'],'junction_aa':[a,a],'v_call':['TRAV1-1*01','TRBV19*01']})
        ref=pl.DataFrame({'complex_id':['1','1'],'gene':['TRA','TRB'],'cdr3':[a,b],
                         'v':['TRAV1-1*01','TRBV19*01']})
    else:
        raw=pl.DataFrame({'query_id':['hit','unknown','invalid','other'],
            'locus':['TRB']*4,'junction_aa':[a,a,'INVALID!',a],
            'v_call':['TRBV19*01','unknown','TRBV19*01','TRBV19*01']})
        ref=pl.DataFrame({'gene':['TRB'],'cdr3':[b],'v':['TRBV19*01']})
    raw=raw.with_columns(pl.lit('human').alias('species'),pl.lit('E').alias('epitope'),
        pl.lit('HLA-A*02').alias('mhc_a'),pl.lit('B2M').alias('mhc_b'),pl.lit('MHCI').alias('mhc_class'))
    if locus=='TRB':
        raw=raw.with_columns(pl.when(pl.col('query_id')=='other').then(pl.lit('mouse'))
                            .otherwise(pl.col('species')).alias('species'))
    raw.write_csv(sample,separator='\t')
    ref=ref.with_columns(pl.lit('HomoSapiens').alias('species'),pl.lit('E').alias('epitope'),
        pl.lit('HLA-A*02:01').alias('mhc_a'),pl.lit('B2M').alias('mhc_b'),pl.lit('MHCI').alias('mhc_class'))
    ref.write_csv(reference,separator='\t')
    args=['tcrdist-neighbours',str(sample),'--vdjdb',str(reference),'--locus',locus,
          '--radius','0','--exclude-exact','--output-prefix',str(prefix)]
    with pytest.raises(SystemExit) as error:
        main(args+['--unbounded-nearest'])
    assert error.value.code==2 and 'requires --targets-from-sample' in capsys.readouterr().err
    args+=['--targets-from-sample']
    assert main(args)==0
    assert pl.read_csv(str(prefix)+'.candidates.tsv',separator='\t').height==0
    assert main(args+['--unbounded-nearest'])==0
    scores=pl.read_csv(str(prefix)+'.candidates.tsv',separator='\t')
    assert scores['distance'].to_list()==[9] and scores['n_neighbours'].to_list()==[0]
    assert scores['mhc_a'].to_list()==['HLA-A*02']
    status=pl.read_csv(str(prefix)+'.queries.tsv',separator='\t',infer_schema_length=0)
    if locus=='paired':
        assert status['query_ida'].to_list()==['raw-a'] and status['query_idb'].to_list()==['raw-b']
    else:
        assert status['query_id'].to_list()==['hit','unknown','invalid','other']
        assert status['status'].to_list()==['matched','invalid_query_or_gene','invalid_query_or_gene','unselected_species']
    # The only reference is now a full junction identity (a full pair in paired mode).
    ref.with_columns(pl.lit(a).alias('cdr3')).write_csv(reference,separator='\t')
    assert main(args+['--unbounded-nearest'])==0
    assert pl.read_csv(str(prefix)+'.candidates.tsv',separator='\t').height==0
    status=pl.read_csv(str(prefix)+'.queries.tsv',separator='\t')
    assert status['status'][0]=='no_neighbours'
