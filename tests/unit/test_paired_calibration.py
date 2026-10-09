"""Small paired calibration contracts: counting units, finite controls and linkage."""

import polars as pl
import pytest
from seqtree import Index
from seqtree.evalue import evalue_result

from vdjmatch.evalue import paired
from vdjmatch.match import PairedVdjdbIndex, search_params


PARAMS = search_params("0,0,0,0")


def observations(ids=("R1", "R2"), epitopes=("Z", "A")):
    return pl.DataFrame({"complex_id": [cid for cid in ids for _ in range(2)],
                         "gene": ["TRA", "TRB"] * len(ids),
                         "cdr3": ["CAVF", "CASSF"] * len(ids),
                         "epitope": [e for e in epitopes for _ in range(2)]})


def controls():
    return Index.build(["CWWF", "CFFF", "CYYF", "CMMF"], "aa")


def query():
    return pl.DataFrame({"cdr3a": ["CAVF"], "cdr3b": ["CASSF"]})


@pytest.mark.parametrize("ids", [(1, 2), ("R1", "R2")])
def test_unique_pairs_preserve_observation_support_and_ties(ids):
    ref = observations(ids)
    idx = PairedVdjdbIndex.build(ref)
    assert idx.n_pairs == 1 and idx.n_observations == 2
    assert idx.reference.equals(ref)
    out = idx.annotate_pairs(query(), controls(), controls(), PARAMS)
    assert out["n_joint"].to_list() == [1]
    assert out["n_alpha"].to_list() == [1]
    assert out["n_beta"].to_list() == [1]
    assert out["epitope"].to_list() == ["A"]
    assert out["E"][0] == pytest.approx(9 / 16)
    assert 0 < out["p_joint"][0] < 1
    reverse = PairedVdjdbIndex.build(ref.reverse()).annotate_pairs(query(), controls(), controls(), PARAMS)
    assert out.equals(reverse)


def test_exact_duplicated_rows_do_not_multiply_observations():
    ref = observations(ids=("R1",), epitopes=("E",))
    idx = PairedVdjdbIndex.build(pl.concat([ref, ref]))
    assert idx.n_pairs == idx.n_observations == 1


def test_cardinality_and_conflicting_annotations():
    ref = observations(ids=("R1",), epitopes=("E",))
    extra = ref.head(1).with_columns(pl.lit("CAYF").alias("cdr3"))
    with pytest.raises(ValueError, match="multiple distinct chains"):
        PairedVdjdbIndex.build(pl.concat([ref, extra]))
    ref = ref.with_columns(pl.Series("epitope", ["E1", "E2"]))
    with pytest.raises(ValueError, match="conflicting"):
        paired.build_paired_ref(ref)


def test_orphan_and_empty_reference():
    ref = observations().head(1)
    idx = PairedVdjdbIndex.build(ref)
    assert idx.n_pairs == 0
    out = idx.annotate_pairs(query(), controls(), controls(), PARAMS)
    assert out["n_joint"][0] == 0 and out["p_joint"][0] == 1
    assert out["epitope"][0] is None
    idx = PairedVdjdbIndex.build(ref.clear())
    assert idx.n_pairs == 0
    assert idx.annotate_pairs(query().clear(), controls(), controls(), PARAMS).height == 0


def test_empty_controls_error_in_all_paths():
    empty = Index.build([], "aa")
    idx = PairedVdjdbIndex.build(observations())
    with pytest.raises(ValueError, match="nonempty"):
        idx.annotate_pairs(query(), empty, controls(), PARAMS)
    with pytest.raises(ValueError, match="nonempty"):
        paired.paired_scan(paired.build_paired_ref(observations()), empty, controls(), [("CAVF", "CASSF")])
    with pytest.raises(ValueError, match="nonempty"):
        paired.pvalue([], [], [], 1, 0, 4)


def test_first_hit_finite_controls_radius_and_existing_count():
    out = paired.pvalue([(2, "E"), (0, "E"), (0, "E")], [1], [], N=10, Ma=100, Mb=200)
    assert out["radius"] == 0 and out["n_pair"] == 2
    assert out["E"] == pytest.approx(10 * 3 / 100 * 3 / 200)
    assert out["p_enrichment"] == evalue_result(2, 9, 10, 100 * 200)["p_enrichment"]
    assert 0 < out["p_enrichment"] < 1
    assert out["rule_of_three_alpha"] and out["rule_of_three_beta"]
    nonzero = paired.pvalue([(0, "E")], [0], [], N=10, Ma=100, Mb=200)
    assert nonzero["E"] == pytest.approx(10 * 1 / 100 * 3 / 200)
    assert not nonzero["rule_of_three_alpha"] and nonzero["rule_of_three_beta"]


def test_first_hit_identity_is_opt_in_and_deduplicates():
    ref = paired.build_paired_ref(observations())
    old, ca, cb = paired.paired_scan(ref, controls(), controls(), [("CAVF", "CASSF")], params=PARAMS)
    assert old == [[(0, "A"), (0, "Z")]]
    new, _, _ = paired.paired_scan(ref, controls(), controls(), [("CAVF", "CASSF")], params=PARAMS, include_identity=True)
    assert all(len(h) == 3 for h in new[0])
    out = paired.pvalue(new[0], ca[0], cb[0], N=1, Ma=4, Mb=4)
    assert out["n_pair"] == 1
    specific = paired.pvalue(new[0], ca[0], cb[0], N=1, Ma=4, Mb=4, epitope="Z")
    assert specific["n_pair"] == 1
    assert paired.pvalue(old[0], ca[0], cb[0], N=2, Ma=4, Mb=4)["n_pair"] == 2


def test_no_hit_and_count_validation():
    out = paired.pvalue([], [], [], N=1, Ma=4, Mb=4)
    assert out["radius"] is None and out["n_pair"] == 0 and out["p_enrichment"] == 1
    with pytest.raises(ValueError, match="counts"):
        paired.pvalue([(0, "E"), (0, "E")], [], [], N=1, Ma=4, Mb=4)
