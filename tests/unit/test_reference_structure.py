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
    from vdjmatch.aggregate.structure import composition
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
    assert not (tmp_path/'profile.punctures.tsv').exists()
    assert calls==[3] and full['junctions'].to_list()==[2,2]
    assert full['graph_edges'].to_list()==[1,1]
    degrees=pl.read_csv(str(prefix)+'.degrees.tsv',separator='\t')
    assert degrees['degree'].to_list()==[1,1]
    assert degrees['junction_count'].to_list()==[2,2]
    assert full['hamming_mean_degree'].to_list()==[1.,1.]
    assert full['hamming_median_degree'].to_list()==[1.,1.]
    assert full['hamming_isolated_fraction'].to_list()==[0.,0.]
    manifest=json.loads((tmp_path/'profile.manifest.json').read_text())
    assert manifest['native_batches']==1
    genes=pl.read_csv(str(prefix)+'.genes.tsv',separator='\t')
    assert set(genes['segment'])=={'V','J'}
    assert genes.filter((pl.col('segment')=='J') & (pl.col('length_scope')=='all_lengths'))['reference_records'].sum()==4
    assert set(manifest['tables'])=={'profile.structures.tsv','profile.degrees.tsv','profile.lengths.tsv','profile.genes.tsv','profile.graph-acquisition.tsv','profile.clonotype-degrees.tsv'}
    assert all(len(manifest[key])==64 for key in ('native_sha256','v_model_sha256','loader_sha256','search_source_sha256'))
    control = tmp_path/'control.tsv'
    pl.DataFrame({'junction_aa':['CASSLAAAYEQYF','CASSLACAYEQYF'],
        'v_call':['TRBV19','TRBV20-1'], 'j_call':['TRBJ2-7']*2,
        'occurrences':[1,3]}).write_csv(control,separator='\t')
    controlled = tmp_path/'controlled'
    calls.clear()
    assert main(['reference-structure','--vdjdb',str(reference),'--locus','TRB',
        '--control',str(control),'--output-prefix',str(controlled)])==0
    assert calls==[3]
    assert pl.read_csv(str(controlled)+'.structures.tsv',separator='\t').equals(full)
    compared=pl.read_csv(str(controlled)+'.genes.tsv',separator='\t')
    assert compared.select(genes.columns).equals(genes)
    v=compared.filter((pl.col('segment')=='V')&(pl.col('length_scope')=='all_lengths'))
    assert v['key_to_control_event_ratio'].to_list()==[4.,4.]
    m=json.loads((tmp_path/'controlled.manifest.json').read_text())
    assert m['control']['control_source_rows']==2 and m['control']['control_valid_rows']==2
    assert m['control']['event_keys']==2 and m['control']['total_events']==4
    assert 'controlled.control-genes.tsv' in m['tables']
    assert 'controlled.information.tsv' in m['tables']
    information = pl.read_csv(str(controlled)+'.information.tsv', separator='\t')
    assert information.height == 8
    assert information.filter(pl.col('segment')=='V')['gene_information_bits'].to_list() == [1.] * 4
    assert information.filter(pl.col('segment')=='J')['gene_information_bits'].to_list() == [0.] * 4
    assert compared.filter(pl.col('segment')=='V')['gene_bits'].to_list() == [1.] * 4
    assert information['length_information_bits'].to_list() == [0.] * 8
    assert m['selection_information']['units'] == 'log2 likelihood-ratio bits per junction/gene key'
    motifs = pl.read_csv(str(controlled)+'.motifs.tsv', separator='\t')
    assert 'controlled.motifs.tsv' in m['tables'] and motifs.height == m['motif_information']['rows']
    assert set(motifs['epitope']) == {'A','B'}
    assert motifs['reference_keys'].unique().to_list() == [2]
    assert motifs['control_keys'].unique().to_list() == [1]
    assert 'motif_information_bits' in motifs.columns
    assert motifs['total_reference_keys'].unique().to_list() == [2]
    assert motifs['total_control_keys'].unique().to_list() == [2]
    assert motifs['stratum_bits'].to_list() == [1.] * motifs.height
    finite = motifs.filter(pl.col('joint_presence_bits').is_not_null())
    assert finite['joint_presence_bits'].to_list() == pytest.approx(
        finite.select(pl.col('stratum_bits') + pl.col('presence_bits')).to_series().to_list())
    mapping = pl.read_csv(str(controlled)+'.germline-map.tsv', separator='\t')
    allele_column = mapping.filter((pl.col('v_gene')=='TRBV19') & (pl.col('loop')=='cdr2') &
        (pl.col('alignment_column')==2) & (pl.col('alphabet')=='full'))
    assert allele_column['projection_status'].item()=='allele_dependent'
    assert allele_column['residue'].item() is None
    factors = pl.read_csv(str(controlled)+'.germline.tsv', separator='\t')
    germline = pl.read_csv(str(controlled)+'.germline-information.tsv', separator='\t')
    assert factors.height==m['germline_information']['factor_rows']
    assert germline.height==m['germline_information']['information_rows']
    assert set(germline['reference_v_keys'])=={2}
    assert set(germline['control_v_keys'])=={2}
    assert germline.filter((pl.col('loop')=='cdr2')&(pl.col('alignment_column')==2))['consensus_information_bits'].to_list()==[0.]*8
    assert {str(controlled.name)+s for s in ('.germline-map.tsv','.germline.tsv','.germline-information.tsv')}<=set(m['tables'])
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
    assert not (tmp_path/'empty.punctures.tsv').exists()
    full=pl.read_csv(str(prefix)+'.structures.tsv',separator='\t')
    assert full['junctions'].item()==0
    assert {'original_weight','gapblock_weight','stem4_weight'} <= set(full.columns)
    degrees=pl.read_csv(str(prefix)+'.degrees.tsv',separator='\t')
    assert degrees['junctions'].item()==0 and degrees['junction_count'].item()==0
    assert degrees['degree'].item() is None
    control=tmp_path/'control.tsv'
    pl.DataFrame({'junction_aa':['CASSLAAAYEQYF'],'v_call':['TRBV19'],'j_call':['TRBJ2-7']}).write_csv(control,separator='\t')
    controlled=tmp_path/'empty-control'
    assert main(['reference-structure','--vdjdb',str(reference),'--locus','TRB','--control',str(control),
                 '--output-prefix',str(controlled)])==0
    for suffix in ('.germline.tsv','.germline-information.tsv'):
        assert pl.read_csv(str(controlled)+suffix,separator='\t').is_empty()


def test_full_graph_identity_survives_unordered_node_join(tmp_path, monkeypatch):
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
    assert not (tmp_path/'profile.punctures.tsv').exists()
    full=pl.read_csv(str(prefix)+'.structures.tsv',separator='\t')
    assert full['graph_edges'].item()==2 and full['components'].item()==1
    degrees=pl.read_csv(str(prefix)+'.degrees.tsv',separator='\t')
    assert degrees.select('degree','junction_count').rows()==[(1,2),(2,1)]


def test_length_coverage_keeps_short_unknown_v_and_duplicate_records(tmp_path, monkeypatch):
    import polars as pl
    import vdjmatch.cluster as cluster
    from vdjmatch.cli.__main__ import main
    reference = tmp_path / 'reference.tsv'
    sequences = ['CAAF', 'CAAF', 'CACF', 'CASSLAAAYEQYF', 'CASSLCAAYEQYF',
                 'CASSLACAYEQYF', 'CASSLAAAAYEQYF', 'C?F']
    pl.DataFrame({'cdr3': sequences, 'v': ['TRBV19'] * 4 + ['unknown'] + ['TRBV19'] * 3,
        'j': ['TRBJ2-7'] * 8, 'epitope': ['E'] * 8,
        'species': ['HomoSapiens'] * 8, 'gene': ['TRB'] * 8}).write_csv(reference, separator='\t')
    calls = []; native = cluster.overlap
    def batched(*args, **kwargs):
        calls.append(len(args[0])); return native(*args, **kwargs)
    monkeypatch.setattr(cluster, 'overlap', batched)
    # Reverse joins to ensure degrees are aligned by local node identity.
    join = pl.DataFrame.join
    def unordered(self, other, *args, **kwargs):
        result = join(self, other, *args, **kwargs)
        return result.reverse() if kwargs.get('on') == 'cdr3' and '_local' in self.columns else result
    monkeypatch.setattr(pl.DataFrame, 'join', unordered)
    prefix = tmp_path / 'profile'
    assert main(['reference-structure', '--vdjdb', str(reference), '--locus', 'TRB',
                 '--output-prefix', str(prefix)]) == 0
    lengths = pl.read_csv(str(prefix) + '.lengths.tsv', separator='\t').sort('junction_length')
    assert calls == [6]
    assert lengths.select('junction_length', 'valid_reference_records', 'modeled_reference_records',
                          'junctions', 'graph_edges', 'isolated_junctions').rows() == [
        (4, 3, 0, 2, 1, 0), (13, 3, 3, 3, 2, 0), (14, 1, 1, 1, 0, 1)]
    assert lengths['mean_degree'].to_list() == [1., 4/3, 0.]
    assert lengths['isolated_fraction'].to_list() == [0., 0., 1.]
    full = pl.read_csv(str(prefix) + '.structures.tsv', separator='\t')
    assert full['v_geometry_imputed_records'].item() == 1
    assert lengths['junctions'].sum() == full['hamming_junctions'].item()
    assert lengths['graph_edges'].sum() == full['hamming_graph_edges'].item()


def test_reference_gene_usage_preserves_records_keys_lengths_and_unresolved_calls():
    import polars as pl
    from vdjmatch.aggregate.structure import reference_gene_usage
    records = pl.DataFrame({'cdr3': ['CAAF', 'CAAF', 'CCAF', 'CAAFF', 'CCAF'],
        'v': ['TRBV1*01', 'TRBV1*02', 'TRBV2', None, 'TRBV1,TRBV2'],
        'j': ['TRBJ1-1', 'TRBJ1-1', 'TRBJ2-1', 'TRBJ2-1', None]})
    result = reference_gene_usage(records, locus='TRB')
    full = result.filter(pl.col('length_scope') == 'all_lengths')
    for segment in ('V', 'J'):
        part = full.filter(pl.col('segment') == segment)
        assert part['reference_records'].sum() == 5
        assert part['record_fraction'].sum() == pytest.approx(1.)
        assert part['key_fraction'].sum() == pytest.approx(1.)
    v = full.filter((pl.col('segment') == 'V') & (pl.col('gene') == 'TRBV1'))
    assert v['reference_records'].item() == 2 and v['junction_gene_keys'].item() == 1
    unknown = full.filter((pl.col('segment') == 'V') & pl.col('gene').is_null())
    assert unknown['reference_records'].item() == unknown['junction_gene_keys'].item() == 2
    lengths = result.filter(pl.col('length_scope') == 'exact_length')
    assert lengths.filter(pl.col('segment') == 'V')['reference_records'].sum() == 5
    assert set(lengths['junction_length']) == {4, 5}
    assert reference_gene_usage(records.head(0), locus='TRB').is_empty()


def test_gene_usage_control_preserves_event_weights_unknowns_and_missing_lengths():
    import polars as pl
    from vdjmatch.aggregate.structure import reference_gene_usage, control_gene_usage, compare_gene_usage
    from vdjmatch.evalue.mock_reference import control_events
    raw = pl.DataFrame({'junction_aa': ['CAAF', 'CAAF', 'CCAF', 'CAAFF'],
        'v_call': ['TRBV1*01', 'TRBV1*02', 'TRBV2', 'TRBV2'],
        'j_call': ['TRBJ1-1'] * 4, 'occurrences': [1, 2, 3, 10]})
    events = control_events(raw, species='human', locus='TRB')
    control = control_gene_usage(events, locus='TRB')
    records = pl.DataFrame({'cdr3': ['CAAF', 'CAAF', 'CCAF', 'CAAFF', 'CAAAFF'],
        'v': ['TRBV1*01', 'TRBV1*02', 'TRBV3', None, 'TRBV1'], 'j': ['TRBJ1-1'] * 5})
    usage = reference_gene_usage(records, locus='TRB')
    result = compare_gene_usage(usage, control)
    assert result.select(usage.columns).equals(usage)
    full = result.filter((pl.col('segment') == 'V') & (pl.col('length_scope') == 'all_lengths'))
    v1 = full.filter(pl.col('gene') == 'TRBV1').row(0, named=True)
    assert v1['control_events'] == 3 and v1['total_control_events'] == 16
    assert v1['control_junction_gene_keys'] == 1 and v1['total_control_junction_gene_keys'] == 3
    assert v1['key_to_control_event_ratio'] == pytest.approx(8 / 3)
    assert v1['record_to_control_event_ratio'] == pytest.approx(3.2)
    absent = full.filter(pl.col('gene') == 'TRBV3').row(0, named=True)
    assert absent['control_events'] == 0 and absent['control_event_fraction'] == 0
    assert absent['key_to_control_event_ratio'] is None and absent['control_status'] == 'gene_not_observed'
    unknown = full.filter(pl.col('gene').is_null()).row(0, named=True)
    assert unknown['reference_records'] == 1 and unknown['total_control_events'] == 16
    assert unknown['control_events'] is None and unknown['control_status'] == 'unresolved_reference_gene'
    exact = result.filter((pl.col('segment') == 'V') & (pl.col('junction_length') == 4) & (pl.col('gene') == 'TRBV1'))
    assert exact['control_events'].item() == 3 and exact['total_control_events'].item() == 6
    assert exact['key_to_control_event_ratio'].item() == 1
    missing = result.filter((pl.col('segment') == 'V') & (pl.col('junction_length') == 6)).row(0, named=True)
    assert missing['control_events'] is None and missing['control_status'] == 'no_control_length'
    with pytest.raises(ValueError, match='positive integer'):
        control_gene_usage(events.with_columns(pl.lit(1.5).alias('event_weight')), locus='TRB')
    from vdjmatch.aggregate.structure import gene_selection_information
    contextual = result.with_columns(pl.lit('human').alias('species'), pl.lit('TRB').alias('locus'),
        pl.lit('TEST').alias('epitope'), pl.lit(None, dtype=pl.String).alias('mhc_a'),
        pl.lit('B2M').alias('mhc_b'), pl.lit('MHCI').alias('mhc_class'))
    factors, information = gene_selection_information(contextual)
    assert factors.height == result.height
    v = information.filter(pl.col('segment') == 'V')
    assert v['reference_keys'].to_list() == [4, 4]
    assert v['unresolved_comparison_keys'].to_list() == [2, 3]
    assert v['unavailable_length_keys'].to_list() == [0, 1]
    assert v['gene_information_bits'].null_count()==0 and v['gene_information_bits'].is_finite().all()
    assert v['joint_information_bits'].null_count()==0 and v['joint_information_bits'].is_finite().all()
    assert v['information_status'].to_list() == ['available_neutral_imputation'] * 2
    assert factors['gene_bits_imputed'].any()
    assert factors.filter(pl.col('gene_bits_imputed'))['gene_bits'].to_list()==[0.]*int(factors['gene_bits_imputed'].sum())


def test_selection_information_matches_count_likelihood_and_conditional_chain_rule():
    import polars as pl
    from vdjmatch.aggregate.structure import (reference_gene_usage, control_gene_usage,
        compare_gene_usage, gene_selection_information)
    probabilities = ([.1, .3, .2, .4], [.3, .2, .4, .1])
    populations = []
    for weights in probabilities:
        rows = []
        for cell, weight in enumerate(weights):
            length = 10 + cell // 2
            for i in range(round(10 * weight)):
                rows.append({'cdr3': 'C' + 'A' * (length - 3) + 'ACDE'[i] + 'F',
                    'v': 'TRBV' + str(1 + cell % 2), 'j': 'TRBJ1-1', 'event_weight': 1 + i})
        populations.append(pl.DataFrame(rows))
    usage = reference_gene_usage(populations[0], locus='TRB')
    control = control_gene_usage(populations[1].rename({'cdr3': 'junction_aa'}), locus='TRB')
    comparison = compare_gene_usage(usage, control).with_columns(pl.lit('human').alias('species'),
        pl.lit('TRB').alias('locus'), pl.lit('TEST').alias('epitope'),
        pl.lit(None, dtype=pl.String).alias('mhc_a'), pl.lit('B2M').alias('mhc_b'), pl.lit('MHCI').alias('mhc_class'))
    factors, information = gene_selection_information(comparison)
    assert factors.select(comparison.columns).sort('segment', 'length_scope', 'junction_length', 'gene').equals(
        comparison.sort('segment', 'length_scope', 'junction_length', 'gene'))
    p, q = probabilities
    joint = sum(a * math.log2(a/b) for a,b in zip(p,q))
    length = sum(sum(p[i:i+2]) * math.log2(sum(p[i:i+2])/sum(q[i:i+2])) for i in (0,2))
    marginal = sum(sum(p[i::2]) * math.log2(sum(p[i::2])/sum(q[i::2])) for i in (0,1))
    row = information.filter((pl.col('segment') == 'V') & (pl.col('length_scope') == 'exact_length')).row(0, named=True)
    assert row['reference_keys'] == 10
    assert row['length_information_bits'] == pytest.approx(length)
    assert row['joint_information_bits'] == pytest.approx(joint)
    assert row['gene_information_bits'] == pytest.approx(joint - length)
    full = information.filter((pl.col('segment') == 'V') & (pl.col('length_scope') == 'all_lengths')).row(0, named=True)
    assert full['gene_information_bits'] == pytest.approx(marginal)
    assert full['joint_information_bits'] == pytest.approx(marginal)
    assert factors.filter(pl.col('gene_bits') < 0).height > 0
    assert information.height == 4 and information['mhc_a'].null_count() == 4
    assert information['unresolved_comparison_keys'].to_list() == [0] * 4
    # Identical key distributions carry zero information despite event-weight imbalance.
    identical = compare_gene_usage(usage, control_gene_usage(
        populations[0].rename({'cdr3': 'junction_aa'}), locus='TRB'))
    for c in ('species', 'locus', 'epitope', 'mhc_a', 'mhc_b', 'mhc_class'):
        identical = identical.with_columns(pl.lit(comparison[c][0], dtype=comparison.schema[c]).alias(c))
    _, zero = gene_selection_information(identical)
    assert zero['joint_information_bits'].to_list() == pytest.approx([0.] * 4)
    with pytest.raises(ValueError, match='unique'):
        gene_selection_information(pl.concat([comparison, comparison.head(1)]))
