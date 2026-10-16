"""Germline projections preserve V-key mass and all modeled allele ambiguity."""
import math
import polars as pl
import pytest


def test_germline_projection_matches_hand_bits_without_merging_v_keys():
    from vdjmatch.aggregate.germline_information import germline_projection_map, iter_germline_information
    from vdjmatch.aggregate.structure import reference_gene_usage, control_gene_usage, compare_gene_usage, gene_selection_information
    model = {
        'TRBV1*01': dict(locus='TRB', cdr1='AR.', cdr2='Q', cdr25='D'),
        'TRBV1*02': dict(locus='TRB', cdr1='AR.', cdr2='H', cdr25='D'),
        'TRBV2*01': dict(locus='TRB', cdr1='AR.', cdr2='Q', cdr25='D'),
        'TRBV3*01': dict(locus='TRB', cdr1='ST.', cdr2='H', cdr25='D'),
        'TRBV4*01': dict(locus='TRB', cdr1='...', cdr2='H', cdr25='D'),
        'TRBV5*01': dict(locus='TRB', cdr1='A*.', cdr2='H', cdr25='D')}
    mapping = germline_projection_map(model, locus='TRB')
    one = mapping.filter((pl.col('v_gene')=='TRBV1') & (pl.col('loop')=='cdr2') & (pl.col('alphabet')=='full'))
    assert one['residue'].item() is None and one['projection_status'].item() == 'allele_dependent'
    assert one['modeled_alleles'].item() == 2
    assert mapping.filter((pl.col('v_gene')=='TRBV4') & (pl.col('loop')=='cdr1'))['projection_status'].unique().to_list() == ['missing_loop_template']
    assert mapping.filter((pl.col('v_gene')=='TRBV5') & (pl.col('loop')=='cdr1') &
        (pl.col('alignment_column')==2))['projection_status'].unique().to_list()==['non_amino_acid_template']
    reference = pl.DataFrame({'cdr3':['CAAF','CAAF'], 'v':['TRBV1*01','TRBV2'], 'j':['TRBJ1-1']*2})
    controls = pl.DataFrame({'junction_aa':['C'+aa+'AF' for aa in 'ACDEFGHI'],
        'v':['TRBV1']+['TRBV2']*3+['TRBV3']*4,'j':['TRBJ1-1']*8,'event_weight':[2]*8})
    usage = reference_gene_usage(reference,locus='TRB').with_columns(pl.lit('human').alias('species'),pl.lit('TRB').alias('locus'),
        pl.lit('TEST').alias('epitope'),pl.lit(None,dtype=pl.String).alias('mhc_a'),pl.lit('B2M').alias('mhc_b'),pl.lit('MHCI').alias('mhc_class'))
    background = control_gene_usage(controls,locus='TRB')
    gene_rows,gene_info = gene_selection_information(compare_gene_usage(usage,background))
    projections = list(iter_germline_information(usage,background,mapping))
    points, info = next((p,i) for p,i in projections if p['loop'][0]=='cdr1' and p['alignment_column'][0]==2 and p['alphabet'][0]=='full')
    row=points.filter(pl.col('length_scope')=='all_lengths').row(0,named=True)
    assert row['residue']=='R' and row['reference_v_keys']==row['reference_stratum_v_keys']==2
    assert row['control_v_keys']==4 and row['control_stratum_v_keys']==8
    assert row['consensus_bits']==row['residue_bits_lower']==pytest.approx(1.)
    assert info['consensus_information_bits'].to_list()==pytest.approx([1.,1.])
    marginal=gene_info.filter((pl.col('segment')=='V')&(pl.col('length_scope')=='all_lengths'))['gene_information_bits'].item()
    assert marginal==pytest.approx(1.+.5*math.log2(4/3)) and marginal>1
    assert info.select((pl.col('joint_information_bits')==pl.col('length_information_bits')+pl.col('consensus_information_bits')).all()).item()
    ambiguous = next(i for p,i in projections if p['loop'][0]=='cdr2' and p['alphabet'][0]=='full')
    assert ambiguous['unavailable_comparison_v_keys'].to_list()==[1,1]
    assert ambiguous['consensus_information_bits'].null_count()==2
    uncertain = next(p for p,i in projections if p['loop'][0]=='cdr2' and p['alphabet'][0]=='full')
    q = uncertain.filter((pl.col('length_scope')=='all_lengths')&(pl.col('residue')=='Q')).row(0,named=True)
    assert q['control_v_keys']==3 and q['control_possible_v_keys']==4
    assert q['consensus_bits']==pytest.approx(math.log2(4/3)) and q['residue_bits_lower']==0.
    for p,i in projections:
        assert i['reference_v_keys'].to_list()==[2,2]
        assert i['control_v_keys'].to_list()==[8,8]
        assert p['consensus_bits'].drop_nulls().is_finite().all()
    unknown = controls.vstack(pl.DataFrame({'junction_aa':['CVAF'],'v':['TRBV999'],
        'j':['TRBJ1-1'],'event_weight':[2]}))
    weakened = next(p for p,i in iter_germline_information(usage, control_gene_usage(unknown,locus='TRB'), mapping)
        if p['loop'][0]=='cdr1' and p['alignment_column'][0]==2 and p['alphabet'][0]=='full')
    row = weakened.filter(pl.col('length_scope')=='all_lengths').row(0,named=True)
    assert row['control_possible_v_keys']==5 and row['control_stratum_v_keys']==9
    assert row['residue_bits_lower']==pytest.approx(math.log2(9/5))
    assert row['residue_bits_lower']<1.
    empty=list(iter_germline_information(usage.head(0),background,mapping))
    assert len(empty)==len(projections)
    for p,i in empty:
        assert p.is_empty() and i.is_empty()
        p.write_csv();i.write_csv()
