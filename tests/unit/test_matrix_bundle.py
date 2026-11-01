"""M4 bundle invariants: neutral geometry, residue/gap effects, and asset validation."""
import hashlib
import json
import numpy as np
import pytest
from vdjmatch.match.matrix_bundle import load_matrix_bundle, sequence_weights, gap_prior
from test_position_profile import frames, counts


def bundle(path, multiplier=1.):
    path.mkdir()
    (path/'residue.tsv').write_text('species\tlocus\tlength\tposition\taa\tmultiplier\n'
        f'human\tTRB\t13\t4\tL\t{multiplier}\n')
    (path/'gaps.tsv').write_text('species\tlocus\tlength\tgap\tmultiplier\n'
        f'human\tTRB\t13\t1\t{multiplier}\n')
    (path/'manifest.json').write_text(json.dumps({'format':'vdjam-bundle-v1','outputs_sha256':{
        name:hashlib.sha256((path/name).read_bytes()).hexdigest() for name in ['residue.tsv','gaps.tsv']}}))
    return load_matrix_bundle(path)


@pytest.mark.parametrize('paired',[False,True])
def test_neutral_bundle_is_exact_for_single_and_paired_native_counts(tmp_path,paired):
    model=bundle(tmp_path/'model')
    q,r=frames(paired)
    cuts=list(range(0,12001,37))
    baseline=counts(q,r,cuts,paired,threads=2)
    result=counts(q,r,cuts,paired,threads=2,matrix_bundle=model)
    np.testing.assert_array_equal(result,baseline)
    assert baseline.sum()>0


def test_residue_and_gap_terms_have_distinct_effects_and_neutral_fallback(tmp_path):
    model=bundle(tmp_path/'model',2.)
    q='CASSLGQAYEQYF'
    weights=sequence_weights(model,'human','TRB',[q,'CASSLGQF'],3)
    assert weights[0]==[0,0,0,100,200,100,100,100,100,100,0,0,0]
    assert weights[1]==[0,0,0,100,100,0,0,0]
    prior=gap_prior(model,'human','TRB',lambda i,d,m:0,2,1)
    assert prior(4,1,14)==4
    assert prior(4,2,15)==3  # undeclared gap retains its baseline cost
    assert prior(4,0,13)==0


def test_bundle_rejects_negative_or_unverified_assets(tmp_path):
    with pytest.raises(ValueError,match='nonnegative'):
        bundle(tmp_path/'negative',-1.)
    modeldir=tmp_path/'hash'
    bundle(modeldir)
    (modeldir/'residue.tsv').write_text('changed')
    with pytest.raises(ValueError,match='hash'):
        load_matrix_bundle(modeldir)


@pytest.mark.parametrize('paired',[False,True])
def test_cli_neutral_bundle_and_explicit_core_preserve_evidence(tmp_path,paired):
    import polars as pl
    from polars.testing import assert_frame_equal
    from test_kernel_evidence_only import fixture
    from vdjmatch.cli.__main__ import main
    from vdjmatch.match.geometry import _substitution_costs,amino_acids
    sample,reference=fixture(tmp_path,paired)
    bundle(tmp_path/'model')
    aa=amino_acids();cost=_substitution_costs('full','blosum62')
    table=tmp_path/'core.tsv'
    table.write_text('aa.1\taa.2\tcost\n'+''.join(f'{x}\t{y}\t{cost[i,j]}\n' for i,x in enumerate(aa) for j,y in enumerate(aa) if i<j))
    args=['evidence',str(sample),'--vdjdb',str(reference),'--locus','paired' if paired else 'TRB','--all-targets','--evidence-only','--exclude-exact']
    for name,extra in [('baseline',[]),('bundle',['--matrix-bundle',str(tmp_path/'model'),'--distance-matrix',str(table)])]:
        assert main(args+extra+['--output-prefix',str(tmp_path/name)])==0
    for lane in ['', '.TRA','.TRB'] if paired else ['']:
        for kind in ['queries','selection','enrichment']:
            left=pl.read_csv(str(tmp_path/'baseline')+lane+'.'+kind+'.tsv',separator='\t',infer_schema_length=None)
            right=pl.read_csv(str(tmp_path/'bundle')+lane+'.'+kind+'.tsv',separator='\t',infer_schema_length=None)
            assert_frame_equal(left,right)
    manifest=json.loads((tmp_path/'bundle.manifest.json').read_text())
    assert manifest['matrix_bundle']['default_changed'] is False
    assert manifest['core_matrix_sha256']==hashlib.sha256(table.read_bytes()).hexdigest()
