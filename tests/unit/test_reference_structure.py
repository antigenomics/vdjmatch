"""Every small graph puncture agrees with an independent component traversal."""
from itertools import combinations
import math
import pytest
from vdjmatch.aggregate.structure import graph_profiles


def brute(n,edges,removed=None):
    nodes=set(range(n))-{removed}; sizes=[]
    while nodes:
        todo=[nodes.pop()]; component=set(todo)
        while todo:
            u=todo.pop()
            for a,b in edges:
                v=b if a==u else a if b==u else None
                if v in nodes:
                    nodes.remove(v);component.add(v);todo.append(v)
        sizes.append(len(component))
    return len(sizes),sum(s==1 for s in sizes),max(sizes,default=0),sum(s*s for s in sizes)


def test_all_four_vertex_graphs_and_articulation_punctures():
    candidates=list(combinations(range(4),2))
    for mask in range(1<<len(candidates)):
        edges=[e for i,e in enumerate(candidates) if mask&(1<<i)]
        full,punctures=graph_profiles(4,edges)
        for removed,row in [(None,full),*enumerate(punctures)]:
            c,i,largest,square=brute(4,edges,removed)
            assert (row['components'],row['isolated_junctions'],row['largest_component']) == (c,i,largest)
            assert row['size_weighted_component_mean'] == square/row['junctions']
            assert row['graph_edges'] == sum(removed not in e for e in edges)
    _,rows=graph_profiles(5,[(0,i) for i in range(1,5)])
    assert rows[0]['components']==4 and rows[0]['isolated_junctions']==4
    assert math.isclose(rows[0]['component_entropy'],math.log(4))


def test_empty_singleton_and_bad_edges():
    full,rows=graph_profiles(0,[])
    assert not rows and full['component_entropy'] is None
    _,rows=graph_profiles(1,[])
    assert rows[0]['junctions']==0 and rows[0]['largest_component_fraction'] is None
    with pytest.raises(ValueError):graph_profiles(2,[(0,2)])


def test_composition_punctures_all_v_annotations():
    import polars as pl
    from vdjmatch.cli.structure import composition
    records=pl.DataFrame({'cdr3':['CAAF','CAAF','CAAF','CAAAF','CAAAAAF'],
        '_allele':['TRBV19*01','TRBV20-1*01','TRBV19*01','TRBV19*02','TRBV20-1*01']})
    junctions=records.select('cdr3').unique(maintain_order=True)
    full,rows=composition(records,junctions)
    assert full['clonotype_keys']==4 and full['reference_records']==5
    for sequence,row in zip(junctions['cdr3'],rows):
        kept=records.filter(pl.col('cdr3')!=sequence)
        keys=kept.unique(['cdr3','_allele']).with_columns(pl.col('_allele').str.split('*').list.first().alias('gene'))
        counts=keys.group_by('gene').len()['len'].to_list()
        expected=-sum(c/keys.height*math.log(c/keys.height) for c in counts)
        assert math.isclose(row['v_entropy'],expected,abs_tol=1e-14)
        assert row['clonotype_keys']==keys.height
        assert row['reference_records']==kept.height
        lengths=[len(s) for s in kept['cdr3'].unique()]
        assert row['length_mean']==sum(lengths)/len(lengths)


def test_cli_uses_one_batch_for_shared_target_junctions(tmp_path,monkeypatch):
    import json
    import polars as pl
    from vdjmatch.cli.__main__ import main
    import vdjmatch.cluster as cluster
    reference=tmp_path/'reference.tsv'
    pl.DataFrame({'cdr3':['CASSLAAAYEQYF','CASSLACAYEQYF','CASSLAAAYEQYF','CASSLADAYEQYF'],
        'v':['TRBV19']*4,'j':['TRBJ2-7']*4,'epitope':['A','A','B','B'],
        'species':['HomoSapiens']*4,'gene':['TRB']*4,
        'mhc_a':['HLA-A*02:01']*4,'mhc_b':['B2M']*4,'mhc_class':['MHCI']*4}).write_csv(reference,separator='\t')
    calls=[];native=cluster.overlap
    def batched(*args,**kwargs):
        calls.append(len(args[0]));return native(*args,**kwargs)
    monkeypatch.setattr(cluster,'overlap',batched)
    prefix=tmp_path/'profile'
    assert main(['reference-structure','--vdjdb',str(reference),'--locus','TRB',
                 '--threads','1','--output-prefix',str(prefix)])==0
    full=pl.read_csv(str(prefix)+'.structures.tsv',separator='\t')
    puncture=pl.read_csv(str(prefix)+'.punctures.tsv',separator='\t')
    assert calls==[3] and full['junctions'].to_list()==[2,2]
    assert puncture.height==4 and puncture['junctions'].to_list()==[1]*4
    assert puncture['graph_edges'].to_list()==[0]*4
    degrees=pl.read_csv(str(prefix)+'.degrees.tsv',separator='\t')
    assert degrees['degree'].to_list()==[1,1]
    assert degrees['junction_count'].to_list()==[2,2]
    assert full['hamming_mean_degree'].to_list()==[1.,1.]
    assert full['hamming_median_degree'].to_list()==[1.,1.]
    assert full['hamming_isolated_fraction'].to_list()==[0.,0.]
    manifest=json.loads((tmp_path/'profile.manifest.json').read_text())
    assert manifest['native_batches']==1
    assert 'profile.degrees.tsv' in manifest['tables']
    assert all(len(manifest[key])==64 for key in ('native_sha256','v_model_sha256','loader_sha256','search_source_sha256'))
    targets=tmp_path/'targets.tsv'
    pl.DataFrame({'query_id':[0,1],'epitope':['A','B'],'species':['human']*2,
        'mhc_a':['HLA-A*02']*2,'mhc_b':['B2M']*2,'mhc_class':['MHCI']*2}).write_csv(targets,separator='\t')
    calls.clear()
    assert main(['reference-structure','--vdjdb',str(reference),'--locus','TRB',
        '--targets-from-sample',str(targets),'--epitope','A','--output-prefix',str(prefix)])==0
    selected=pl.read_csv(str(prefix)+'.structures.tsv',separator='\t')
    assert calls==[2] and selected.height==1
    assert selected['epitope'].to_list()==['A'] and selected['reference_records'].to_list()==[2]
    assert selected['mhc_a'].to_list()==['HLA-A*02']


def test_empty_reference_exports_typed_weights_and_zero_node_histogram(tmp_path):
    import polars as pl
    from vdjmatch.cli.__main__ import main
    reference=tmp_path/'invalid.tsv';prefix=tmp_path/'empty'
    pl.DataFrame({'cdr3':['C?F'], 'v':['TRBV19'], 'j':['TRBJ2-7'],
        'epitope':['E'], 'species':['HomoSapiens'], 'gene':['TRB']}).write_csv(reference,separator='\t')
    assert main(['reference-structure','--vdjdb',str(reference),'--locus','TRB',
        '--output-prefix',str(prefix)])==0
    punctures=pl.read_csv(str(prefix)+'.punctures.tsv',separator='\t')
    assert punctures.is_empty()
    assert {'original_weight','gapblock_weight','stem4_weight','mixture_only_puncture'} <= set(punctures.columns)
    degrees=pl.read_csv(str(prefix)+'.degrees.tsv',separator='\t')
    assert degrees['junctions'].item()==0 and degrees['junction_count'].item()==0
    assert degrees['degree'].item() is None


def test_puncture_identity_survives_unordered_node_join(tmp_path, monkeypatch):
    import polars as pl
    from vdjmatch.cli.__main__ import main
    center='CASSLAAAYEQYF';leaves=['CASSLCAAYEQYF','CASSLACAYEQYF']
    reference=tmp_path/'path.tsv';prefix=tmp_path/'profile'
    pl.DataFrame({'cdr3':[center,*leaves], 'v':['TRBV19']*3, 'j':['TRBJ2-7']*3,
        'epitope':['E']*3, 'species':['HomoSapiens']*3, 'gene':['TRB']*3}).write_csv(reference,separator='\t')
    join=pl.DataFrame.join
    def unordered(self, other, *args, **kwargs):
        result=join(self,other,*args,**kwargs)
        if kwargs.get('on')=='cdr3' and '_local' in self.columns and '_node' in other.columns:
            return result.reverse()
        return result
    monkeypatch.setattr(pl.DataFrame,'join',unordered)
    assert main(['reference-structure','--vdjdb',str(reference),'--locus','TRB',
        '--output-prefix',str(prefix)])==0
    punctures=pl.read_csv(str(prefix)+'.punctures.tsv',separator='\t')
    hub=punctures.filter(pl.col('excluded_junction')==center)
    assert hub['graph_edges'].item()==0 and hub['components'].item()==2
    ends=punctures.filter(pl.col('excluded_junction').is_in(leaves))
    assert ends['graph_edges'].to_list()==[1,1] and ends['components'].to_list()==[1,1]
