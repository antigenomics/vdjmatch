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


# --- paired alpha/beta E-value ---
def test_paired_evalue():
    from seqtree import Index
    from vdjmatch.match import PairedVdjdbIndex, search_params, load_vdjam
    vdjdb = pl.DataFrame({
        "species": ["HomoSapiens"] * 4,
        "complex_id": [1, 1, 2, 2],
        "gene": ["TRA", "TRB", "TRA", "TRB"],
        "cdr3": ["CAVRDSNYQLIW", "CASSIRSSYEQYF", "CAGHTGNQFYF", "CASSPGTGGYGYTF"],
        "epitope": ["GILGFVFTL", "GILGFVFTL", "NLVPMVATV", "NLVPMVATV"],
    })
    pidx = PairedVdjdbIndex.build(vdjdb, species="HomoSapiens")
    assert pidx.n_pairs == 2
    # controls of unrelated random clonotypes (low background)
    ctrl = Index.build(["CAAAAAAAAF", "CBBBBBBBBF".replace("B", "G"), "CWWWWWWWWF"], "aa")
    pairs = pl.DataFrame({"cdr3a": ["CAVRDSNYQLIW"], "cdr3b": ["CASSIRSSYEQYF"]})
    res = pidx.annotate_pairs(pairs, ctrl, ctrl, search_params("1,0,0,1", matrix=load_vdjam()))
    assert res["n_joint"][0] >= 1                 # the true complex matches both chains
    assert res["epitope"][0] == "GILGFVFTL"
    assert res["p_joint"][0] < 1.0                # enriched over background


# --- region-aware scoring (germline-retention weighting) ---
def test_region_weights_downweight_germline_flanks():
    from vdjmatch.match import regions
    ret = regions.load_retention()
    # CASSIRSSYEQYF (len 13), TRBV19 / TRBJ2-7: CASS prefix + (Y)EQYF suffix are germline-retained
    w = regions.position_weights(13, "TRBV19", "TRBJ2-7", "TRB", ret)
    assert len(w) == 13
    assert w[0] < 0.2 and w[1] < 0.2          # C, A of CASS -> germline -> low weight
    assert w[-1] < 0.2                         # terminal F (J anchor) -> germline -> low weight
    assert max(w[4:8]) > 0.7                   # NDN core -> high weight
    # unknown genes -> full weight everywhere (no germline credit)
    assert all(x == 1.0 for x in regions.position_weights(10, "NOPE", "NOPE", "TRB", ret))


def test_significance_weights_centre_heavy_and_pssm_builds():
    from seqtree import Index, SearchParams
    from vdjmatch.match import regions
    w = regions.significance_weights(12)
    assert len(w) == 12
    c = len(w) // 2
    assert w[c] > w[0] and w[c] > w[-1]        # central substitutions weighted above V/J borders
    # native seqtree PSSM path: build + search must run with the positional matrix attached
    pm = regions.significance_pssm(12)
    idx = Index.build(["CASSIRSSYEQYF", "CASSIRSAYEQYF"], "aa")
    p = SearchParams(max_subs=2, max_total_edits=2, max_penalty=10 ** 9, engine="seqtrie")
    p.pos_matrix = pm
    hits = idx.search_batch(["CASSIRSSYEQYF"], p, 0)[0]
    assert any(h.n_subs == 0 for h in hits)    # exact self is found


def test_region_aware_engine_path():
    idx = VdjdbIndex.build(_tiny_vdjdb(), species="HomoSapiens")
    q = pl.DataFrame({"cdr3": ["CASSIRSSYEQYF"], "v": ["TRBV19"], "j": ["TRBJ2-7"],
                      "locus": ["TRB"], "count": [1]})
    hits = idx.annotate(q, search_params("1,0,0,1", matrix=scoring.load_vdjam()),
                        gene="TRB", align=True, region_aware=True)
    assert "region_score" in hits.columns
    assert all(s >= 0 for s in hits["region_score"])
