"""Unit tests for vdjmatch core: scope parsing, CIGAR, column aliasing, scoring, engine, cluster."""
import polars as pl
import pytest

from vdjmatch.match import scope, cigar, scoring, VdjdbIndex, search_params
from vdjmatch.io import columns
from vdjmatch import cluster, aggregate


# --- scope ---
def test_parse_scope():
    assert scope.parse_scope("1") == (1, 0, 0, 1)
    assert scope.parse_scope("2,1,2") == (2, 1, 1, 2)        # s, indel, total
    assert scope.parse_scope("2,1,1,3") == (2, 1, 1, 3)
    assert scope.parse_scope("3,0,0,3") == (3, 0, 0, 3)


def test_search_params_sets_pos_matrix():
    p = search_params("2,0,0,2", engine="seqtm")
    assert p.max_subs == 2 and p.max_total_edits == 2 and p.engine == "seqtm"


# --- cigar ---
def test_cigar_roundtrip():
    assert cigar.to_cigar("MMMMM") == "5="
    assert cigar.to_cigar("MMMSMM") == "3=1X2="
    assert cigar.to_cigar("MMIID") == "2=2I1D"
    assert cigar.to_cigar("") == ""
    assert cigar.match_line("MMSMI") == "|| | "          # | for match, space otherwise
    assert cigar.counts("MMSID") == {"matches": 2, "subs": 1, "ins": 1, "dels": 1}


# --- column aliasing ---
def test_normalize_query_aliases_and_locus():
    df = pl.DataFrame({"junction_aa": ["CASSF", "cassl", "CXZ*"], "v_gene": ["TRBV2", "TRBV5", "TRBV1"],
                       "j_gene": ["TRBJ2-1"] * 3})
    out = columns.normalize_query(df)
    assert out.columns == ["cdr3", "v", "j", "locus", "count", "pair_id"]
    assert out["locus"].to_list()[0] == "TRB"            # derived from v gene
    assert out["cdr3"].to_list() == ["CASSF", "CASSL"]   # upper-cased; non-AA 'CXZ*' dropped
    assert out["count"].to_list() == [1, 1]              # default count


def test_normalize_query_requires_cdr3():
    with pytest.raises(ValueError):
        columns.normalize_query(pl.DataFrame({"v_gene": ["TRBV2"]}))


# --- scoring ---
def test_load_vdjam_builds():
    m = scoring.load_vdjam()
    assert m.size() == 24


# --- engine end-to-end (tiny) ---
def _tiny_vdjdb():
    return pl.DataFrame({
        "gene": ["TRB", "TRB", "TRB"],
        "cdr3": ["CASSIRSSYEQYF", "CASSIRSSYEQYY", "CASSPGTGYEQFF"],
        "v": ["TRBV19", "TRBV19", "TRBV5"], "j": ["TRBJ2-7"] * 3,
        "epitope": ["GILGFVFTL", "GILGFVFTL", "NLVPMVATV"],
        "mhc_a": ["HLA-A*02:01"] * 3, "mhc_b": ["B2M"] * 3, "mhc_class": ["MHCI"] * 3,
        "antigen_gene": ["M1", "M1", "pp65"], "antigen_species": ["Flu", "Flu", "CMV"],
        "vdjdb_score": [2, 1, 2], "complex_id": [0, 0, 0], "species": ["HomoSapiens"] * 3,
    })


def test_engine_annotate_and_aggregate():
    idx = VdjdbIndex.build(_tiny_vdjdb(), species="HomoSapiens")
    assert idx.genes == ["TRB"]
    q = pl.DataFrame({"cdr3": ["CASSIRSSYEQYF"], "v": ["TRBV19"], "j": ["TRBJ2-7"],
                      "locus": ["TRB"], "count": [5]})
    hits = idx.annotate(q, search_params("1,0,0,1", matrix=scoring.load_vdjam()),
                        gene="TRB", align=True)
    # exact self + the 1-substitution neighbour, both GILGFVFTL
    assert set(hits["db_cdr3"]) == {"CASSIRSSYEQYF", "CASSIRSSYEQYY"}
    assert hits.filter(pl.col("n_subs") == 0)["epitope"].to_list()[0] == "GILGFVFTL"
    assert "cigar" in hits.columns
    call = aggregate.best_call(hits)
    assert call["epitope"].to_list()[0] == "GILGFVFTL"


# --- cluster ---
def test_overlap_within_drops_self():
    seqs = ["CASSF", "CASSL", "CWWWW"]
    pairs = cluster.overlap(seqs, scope="1,0,0,1")
    assert pairs.filter(pl.col("a_idx") == pl.col("b_idx")).height == 0   # no self-pairs
    # CASSF ~ CASSL (1 sub) should pair
    got = {(r["a_cdr3"], r["b_cdr3"]) for r in pairs.iter_rows(named=True)}
    assert ("CASSF", "CASSL") in got
