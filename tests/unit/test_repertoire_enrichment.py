import numpy as np
import polars as pl
import pytest
from vdjmatch.aggregate.repertoire_enrichment import matched_enrichment


def test_exact_bin_count_matching_and_both_endpoints():
    observations=pl.DataFrame({"sample_id":["s"]*3,"query_id":[0,1,2],"count":[7,3,5],"pgen_bin":[-8,-8,-12]})
    control=pl.DataFrame({"query_id":[3,4,5],"event_weight":[1,1,1],"pgen_bin":[-8,-8,-12]})
    targets=pl.DataFrame({"epitope":["X","Y"]})
    # Duplicate reference observations must not multiply a matched receptor.
    hits=pl.DataFrame({"query_id":[0,0,1,2,3],"epitope":["X","X","X","X","X"]})
    summary,mocks,audit=matched_enrichment(observations,control,hits,targets,replicates=200,seed=7)
    x=summary.filter(pl.col("epitope")=="X").row(0,named=True)
    y=summary.filter(pl.col("epitope")=="Y").row(0,named=True)
    assert (x["matched_records"],x["matched_reads"])==(3,15)
    assert x["p_records"]==x["p_reads"]==1/201
    assert x["excess_records"]==2 and x["enrichment_records"]==3
    assert y["matched_records"]==y["matched_reads"]==0 and y["p_records"]==y["p_reads"]==1
    assert mocks.height==400 and set(mocks["mock_records"])=={0,1}
    assert set(mocks.filter(pl.col("epitope")=="X")["mock_reads"])=={3,7}
    assert audit["requested_records"].sum()==3 and audit["exact_histogram_and_counts"].all()
    assert mocks.equals(matched_enrichment(observations,control,hits,targets,replicates=200,seed=7)[1])
    with pytest.raises(ValueError,match="Pgen bin -12"):
        matched_enrichment(observations,control.filter(pl.col("pgen_bin")==-8),hits,targets,replicates=200)
    with pytest.raises(ValueError,match="Pgen bin -8"):
        matched_enrichment(observations,control.filter(pl.col("query_id")!=4),hits,targets,replicates=200)


@pytest.mark.parametrize("mismatches",[0,1])
def test_repertoire_enrichment_cli_end_to_end(tmp_path,mismatches):
    import os,subprocess,sys,json
    from pathlib import Path
    sample=tmp_path/"sample.tsv"
    sample.write_text("junction_aa\tv_call\tj_call\tlocus\tduplicate_count\nCASSF\tTRBV2\tTRBJ1-1\tTRB\t3\nCASSF\tTRBV2\tTRBJ1-1\tTRB\t4\n")
    metadata=tmp_path/"samples.tsv";metadata.write_text("sample_id\tpath\nS0\tsample.tsv\n")
    control=tmp_path/"control.tsv";control.write_text("junction_aa\tv\tj\toccurrences\nCASSF\tTRBV2\tTRBJ1-1\t100\n")
    reference=tmp_path/"reference.tsv";reference.write_text("cdr3\tgene\tv\tj\tspecies\tepitope\tmhc_a\tmhc_b\tmhc_class\nCASSF\tTRB\tTRBV2\tTRBJ1-1\tHomoSapiens\tLLWNGPMAV\tHLA-A*02:01\tB2M\tMHCI\n")
    prefix=tmp_path/"fresh"
    command=[sys.executable,"-m","vdjmatch.cli","repertoire-enrichment",str(metadata),"--vdjdb",str(reference),"--control",str(control),"--threads","1","--pgen-mismatches",str(mismatches),"--output-prefix",str(prefix)]
    r=subprocess.run(command,text=True,capture_output=True,timeout=60)
    assert r.returncode==0,r.stdout+r.stderr
    summary=pl.read_csv(str(prefix)+".summary.tsv",separator="\t")
    assert summary["matched_records"].item()==2 and summary["matched_reads"].item()==7
    assert summary["p_records"].item()==summary["p_reads"].item()==1
    manifest=json.loads(Path(str(prefix)+".manifest.json").read_text())
    assert manifest["source_frozen"] and manifest["inputs_frozen"] and manifest["native_search_batches"]==1
    assert manifest["automatic_detection_threshold"] is None
    assert manifest["pgen_mismatches"]==mismatches


def test_tcr_source_counts_are_exact(tmp_path):
    from vdjmatch import io
    p=tmp_path/"source.tsv"
    p.write_text("Read.count\tCDR3.nucleotide.sequence\tCDR3.amino.acid.sequence\tbestVGene\tbestJGene\n9007199254740993\tTGTGCCAGCAGCTTT\tCASSF\tTRBV2\tTRBJ1-1\n1\tTGTGCCAGCAGCTTT\tCASSF\tTRBV2,TRBV3\t\n")
    frame=io.read_rearrangement(p,valid_aa=False)
    assert frame["count"].to_list()==[9007199254740993,1]
    assert frame["v"].to_list()==["TRBV2","TRBV2,TRBV3"]


def test_pgen_batch_respects_requested_mass(monkeypatch):
    from vdjmatch.evalue.mock_reference import compute_pgen
    from vdjtools.model import native
    observed=[]
    def batch(model,sequences,*,mismatches,threads):
        observed.append((sequences,mismatches,threads))
        return [0.001*(mismatches+1)]*len(sequences)
    monkeypatch.setattr(native,"pgen_aa_batch",batch)
    keys=pl.DataFrame({"junction_aa":["CASSF","CASSF"]})
    for mm in (0,1):
        result=compute_pgen(keys,keys.head(0),mismatches=mm,threads=2)
        assert result.height==1 and result["pgen"].item()==0.001*(mm+1)
    assert observed==[(["CASSF"],0,2),(["CASSF"],1,2)]


def test_v_ambiguity_policy():
    from vdjmatch.cli.repertoire_enrichment import _single_v
    assert _single_v(pl.DataFrame({"v":["TRBV2*01,TRBV2*02"]}),"v")["v"].item()=="TRBV2"
    for value in ["TRBV2,TRBV3","TRBV2*01,TRBV3*01"]:
        with pytest.raises(ValueError,match="unresolved multi-gene V"):
            _single_v(pl.DataFrame({"v":[value]}),"v")
