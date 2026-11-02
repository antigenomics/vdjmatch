"""Exact count identities and origin/record accounting, without scientific reruns."""

import numpy as np
import polars as pl
import pytest
from vdjmatch.match.scoring_audit import column_pair_counts, origin_block_matrices


def test_origin_pair_counter_matches_explicit_distinct_pairs():
    encoded = np.array([[0, 0], [1, 0], [0, 0]], dtype=np.uint8)
    weights = np.array([[0.2, 0], [0.5, 0], [0.7, 0]])
    counts, _, stats = column_pair_counts(encoded, np.arange(3), weights)
    assert np.allclose(counts[:2, :2], [[0.28, 0.45], [0.45, 0]])
    assert np.isclose(stats["column_pairs"], 0.59)
    assert column_pair_counts(encoded, np.array([0, 0, 1]), weights)[0][
        0, 1
    ] == pytest.approx(0.175)
    with pytest.raises(ValueError):
        column_pair_counts(encoded, np.arange(3), weights * np.nan)
    assert (
        column_pair_counts(
            np.zeros((0, 2), dtype=int), np.array([], dtype=int), np.zeros((0, 2))
        )[0].sum()
        == 0
    )


def _inputs():
    seqs = ["AAA", "AAC", "CAA", "CAC", "DDD"]
    clones = pl.DataFrame(
        dict(
            query_id=list(range(5)),
            cdr3=seqs,
            v=["V"] * 5,
            j=["J"] * 5,
            epitope=["E"] * 5,
        )
    )
    rows = [
        dict(
            query_id=i,
            position=p,
            n_probability=(None if i == 4 else 0.0 if p == 1 else 1.0),
        )
        for i in range(5)
        for p in range(3)
    ]
    return clones, pl.DataFrame(rows)


def test_origin_blocks_remove_fixed_identities_and_retain_missing_cells():
    clones, origins = _inputs()
    matrix, blocks, support, pairs, _ = origin_block_matrices(clones, origins)
    assert support.filter(~pl.col("origin_supported"))["query_id"].to_list() == [4]
    assert blocks.height == 1 and blocks["clonotypes"][0] == 4
    assert pairs["n_weighted_ordered_site_pairs"].sum() == 24
    assert origin_block_matrices(clones.reverse(), origins.reverse())[0].equals(matrix)
    pooled = matrix.filter(pl.col("weighting") == "pooled")
    # Two N columns, each with 2 A and 2 C: ordered AA=CC=4, AC=CA=8.
    assert (
        pooled.filter((pl.col("a") == "A") & (pl.col("b") == "A"))["pair_weight"][0]
        == 4
    )
    ac = pooled.filter((pl.col("a") == "A") & (pl.col("b") == "C")).row(0, named=True)
    assert ac["pair_weight"] == 8 and ac["supported"] and ac["distance_bits"] == 0
    assert (
        pooled.filter((pl.col("a") == "D") & (pl.col("b") == "D"))["distance_bits"][0]
        is None
    )
    assert (
        np.isclose(pooled["q"].sum(), 1.0)
        and blocks["n_weighted_unordered_site_pairs"][0] == 12
    )
    # A smaller epitope has AA=2, AC=CA=1: equal-epitope and pooled differ.
    second = clones.filter(pl.col("query_id") < 2).with_columns(
        pl.lit("F").alias("epitope")
    )
    repeated, _, _, _, _ = origin_block_matrices(pl.concat([clones, second]), origins)
    aa = repeated.filter((pl.col("a") == "A") & (pl.col("b") == "A"))
    assert aa.filter(pl.col("weighting") == "equal_epitope")["q"][0] == pytest.approx(
        1 / 3
    )
    assert aa.filter(pl.col("weighting") == "pooled")["q"][0] == pytest.approx(3 / 14)
    with pytest.raises(ValueError, match="every codon"):
        origin_block_matrices(clones, origins.slice(1))
    with pytest.raises(ValueError, match="unique"):
        origin_block_matrices(pl.concat([clones, clones.slice(0, 1)]), origins)


def test_origin_cli_preserves_source_units_and_zero_support(tmp_path, monkeypatch):
    import argparse
    import json
    import vdjtools.model
    from vdjmatch import db
    from vdjmatch.cli.scoring_audit import register

    seqs = [
        "CAVRDSNYQLIW",
        "CAVRDSNYQLIW",
        "CAVRESNYQLIW",
        "CAVRDSNYQLIF",
        "CAVRDSNYQLIW",
    ]
    frame = pl.DataFrame(
        dict(
            cdr3=seqs,
            v=["TRAV1-1*01"] * 4 + ["TRAV1-1*01,TRAV2*01"],
            j=["TRAJ1*01"] * 5,
            species=["HomoSapiens"] * 5,
            gene=["TRA"] * 5,
            epitope=["E"] * 5,
            mhc_a=["A"] * 5,
            mhc_b=["B"] * 5,
            mhc_class=["MHCI"] * 5,
            reference_valid=[True] * 5,
        )
    )
    short = frame.slice(0, 1).with_columns(pl.lit("C").alias("cdr3"))
    frame = pl.concat([frame, short])
    tail = frame[[0, 2]].with_columns(pl.lit("TAIL").alias("epitope"))
    complete = pl.concat([frame, tail])
    monkeypatch.setattr(db, "load", lambda *args, **kwargs: complete)
    import vdjmatch.match.scoring_audit as owner

    original = owner.origin_block_matrices

    def count_with_complete_labels(clones, origins, known_labels):
        assert known_labels.filter(pl.col("epitope") == "TAIL").height == 2
        return original(clones, origins, known_labels)

    monkeypatch.setattr(owner, "origin_block_matrices", count_with_complete_labels)
    monkeypatch.setattr(db, "provenance", lambda p: dict(test_raw=str(p)))

    def posterior(model, cdr3_aas, **kwargs):
        assert (
            len(cdr3_aas) == 4
        )  # Duplicate source records do not become extra origin queries.
        assert all("," not in v for v in kwargs["v"])
        return pl.DataFrame(
            [
                dict(
                    row=i,
                    position=p,
                    n_probability=(
                        None if s.endswith("F") else 0.5 if 0 < p < len(s) - 1 else 0.0
                    ),
                    supported=not s.endswith("F"),
                )
                for i, s in enumerate(cdr3_aas)
                for p in range(len(s))
            ]
        )

    monkeypatch.setattr(
        vdjtools.model, "n_region_posterior_batch", posterior, raising=False
    )
    raw = tmp_path / "raw.tsv"
    raw.write_text("fixture source")
    parser = argparse.ArgumentParser()
    register(parser.add_subparsers())
    prefix = tmp_path / "origin"
    args = parser.parse_args(
        [
            "scoring-audit",
            "--task",
            "origin-matrix",
            "--pmhc",
            str(raw),
            "--locus",
            "TRA",
            "--min-records",
            "3",
            "--output-prefix",
            str(prefix),
        ]
    )
    assert args.func(args) == 0
    support = pl.read_csv(str(prefix) + ".support.tsv", separator="\t")
    assert support["source_records"].sum() == 6
    observed = {r["status"]: r["source_records"] for r in support.to_dicts()}
    assert observed == dict(
        origin_supported=4, zero_sequence_likelihood=1, unresolved_call=1
    )
    calls = pl.read_csv(str(prefix) + ".calls.tsv", separator="\t")
    assert calls.filter(pl.col("status") == "ambiguous_call").height == 1
    manifest = json.loads((tmp_path / "origin.provenance.json").read_text())
    assert manifest["origin_queries"] == 4 and manifest["source_records"] == 6
    assert pl.read_csv(str(prefix) + ".matrix.tsv", separator="\t").height == 800
    band = pl.read_csv(str(prefix) + ".length-band-support.tsv", separator="\t")
    assert band.filter(pl.col("legacy_length_lt8"))["source_records"].sum() == 1

    def neighbour_route(model, queries, membership, origins, known_labels, *, source_membership, threads, pair_data):
        assert queries.height == 4 and source_membership.height == 4
        assert source_membership['reference_id'].null_count() == 4
        assert known_labels.filter(pl.col('epitope') == 'TAIL').height == 2
        return {'neighbour-route': pl.DataFrame({'queries': [queries.height]})}

    monkeypatch.setattr(owner, 'origin_neighbour_tables', neighbour_route)
    rate_calls = []
    def rate_route(model, queries, membership, origins, *, threads, pair_data):
        assert pair_data is not None
        rate_calls.append(queries.height)
        return {'origin-rate-counts': pl.DataFrame({'queries': [queries.height]})}
    monkeypatch.setattr(owner, 'origin_rate_tables', rate_route)
    args.origin_rates = True
    args.origin_neighbours = True
    args.output_prefix = str(tmp_path / 'neighbours')
    assert args.func(args) == 0
    assert pl.read_csv(str(args.output_prefix) + '.neighbour-route.tsv', separator='\t')['queries'][0] == 4
    assert rate_calls == [4]


def test_matched_background_uses_disjoint_complete_labels_and_same_N_opportunity():
    clones = pl.DataFrame(
        dict(
            query_id=[0, 1, 2, 3],
            cdr3=["AAA", "CAA", "DAA", "EAA"],
            v=["V"] * 4,
            j=["J"] * 4,
            epitope=["E", "E", "F", "G"],
        )
    )
    weights = [0.2, 0.5, 0.7, 0.9]
    origins = pl.DataFrame(
        [
            dict(query_id=q, position=p, n_probability=weights[q] if p == 0 else 0.0)
            for q in range(4)
            for p in range(3)
        ]
    )
    labels = pl.concat(
        [
            clones.select("query_id", "epitope"),
            pl.DataFrame(dict(query_id=[0, 2], epitope=["TAIL", "TAIL"])),
        ]
    )
    *_, contrast = origin_block_matrices(clones, origins, labels)
    pairs = contrast["recognition-pairs"].filter(pl.col("epitope") == "E")
    assert pairs["positive_pair_weight"].sum() == pytest.approx(0.2)
    assert pairs["background_pair_weight"].sum() == pytest.approx(0.2)
    assert (
        pairs.filter((pl.col("a") == "A") & (pl.col("b") == "D"))[
            "background_pair_weight"
        ][0]
        == 0
    )
    assert pairs.filter((pl.col("a") == "C") & (pl.col("b") == "D"))[
        "background_pair_weight"
    ][0] == pytest.approx(0.35 * 0.2 / 1.96)
    support = contrast["recognition-support"]
    assert support.height == 1 and support["position"][0] == 0
    assert support["background_ordered_pair_weight"][0] == pytest.approx(1.96)
    # A common known label forbids all background pairs, retaining unmatched mass.
    common = pl.concat(
        [
            labels,
            clones.select("query_id").with_columns(pl.lit("COMMON").alias("epitope")),
        ]
    )
    *_, unmatched = origin_block_matrices(clones, origins, common)
    assert not unmatched["recognition-support"]["matched"].any()
    assert unmatched["recognition-support"][
        "positive_ordered_pair_weight"
    ].sum() == pytest.approx(0.2)
    assert unmatched["recognition-pairs"]["positive_pair_weight"].sum() == 0
    assert unmatched["recognition-pairs"]["background_pair_weight"].sum() == 0


def test_log_role_symmetry_removes_independent_composition_effect():
    from vdjmatch.match.scoring_audit import (
        similarity_from_counts,
        costs_from_similarity,
    )

    p = np.array([0.8, 0.2])
    g = np.array([0.2, 0.8])
    positive = np.outer(p, p)
    directed = np.outer(p, g)
    unordered, _ = similarity_from_counts(
        positive, (directed + directed.T) / 2, lambda_bits=1.0, round_scores=False
    )
    oriented, _ = similarity_from_counts(
        positive, directed, lambda_bits=1.0, round_scores=False
    )
    assert costs_from_similarity(unordered)[0, 1] > 2.0
    assert np.allclose(
        costs_from_similarity((oriented + oriented.T) / 2), 0.0, atol=1e-14
    )


def test_exact_flank_opportunity_uses_native_batch_and_preserves_zero_mass():
    from vdjtools.model.io import from_germline
    from vdjmatch.match.scoring_audit import n_region_residue_opportunities

    germline = pl.DataFrame([
        dict(allele=a, gene=a.split('*')[0], segment=k, sequence=s,
             cdr3_anchor=-1, functional=True, full_germline='')
        for a, k, s in [('TOYV1*01', 'V', 'TGT'), ('TOYJ1*01', 'J', 'TTT')]
    ])
    model = from_germline(germline, locus='TOY', ins_max=3, palindrome_max=0)
    rho = n_region_residue_opportunities(
        model,
        ["CKF", "CKF", "AAA"],
        [1, 0, 1],
        ["TOYV1*01"] * 3,
        ["TOYJ1*01"] * 3,
        threads=1,
    )
    degeneracy = np.array([4, 2, 2, 2, 2, 4, 2, 3, 2, 6, 1, 2, 4, 2, 6, 6, 4, 4, 1, 2])
    assert np.allclose(rho[0], degeneracy / 61)
    assert np.isnan(
        rho[1:]
    ).all()  # Fixed germline and unsupported flanks stay unidentifiable.
    opportunity = 2 * np.outer(rho[0], rho[0]) / (1 - rho[0] @ rho[0])
    assert np.isclose(opportunity[np.triu_indices(20, 1)].sum(), 1)
    assert n_region_residue_opportunities(model, [], [], [], [], threads=1).shape == (
        0,
        20,
    )
    with pytest.raises(ValueError):
        n_region_residue_opportunities(model, ["CKF"], [3], ["TOYV1*01"], ["TOYJ1*01"])


def test_conditional_exchange_fit_separates_residue_preferences():
    from scipy.special import softmax
    from vdjmatch.match.scoring_audit import fit_exchangeability

    rng = np.random.default_rng(71)
    i, j = np.triu_indices(20, 1)
    groups = np.repeat(np.arange(3), 12)
    rho = rng.dirichlet(np.ones(20) * 3, len(groups))
    baseline = rho[:, i] * rho[:, j]
    baseline /= baseline.sum(axis=1, keepdims=True)
    # Purely additive residue selection must not become pair compatibility.
    eta = rng.normal(0, 0.3, (3, 20))
    eta -= eta.mean(axis=1, keepdims=True)
    target = softmax(np.log(baseline) + eta[groups][:, i] + eta[groups][:, j], axis=1)
    beta, report = fit_exchangeability(
        1000 * target, baseline, balanced=False, prior_mass=0, residue_groups=groups
    )
    fitted_eta = np.asarray(report["residue_log_preferences"])
    predicted = softmax(
        np.log(baseline) + beta + fitted_eta[groups][:, i] + fitted_eta[groups][:, j],
        axis=1,
    )
    assert np.allclose(predicted, target, atol=2e-6)
    assert np.max(np.abs(beta)) < 2e-4
    assert np.allclose(fitted_eta, eta, atol=2e-4)

    incidence = np.zeros((190, 20))
    incidence[np.arange(190), i] = incidence[np.arange(190), j] = 1
    truth = rng.normal(0, 0.2, 190)
    truth -= incidence @ np.linalg.lstsq(incidence, truth, rcond=None)[0]
    target = softmax(
        np.log(baseline) + eta[groups][:, i] + eta[groups][:, j] + truth, axis=1
    )
    recovered, _ = fit_exchangeability(
        1000 * target, baseline, balanced=False, prior_mass=0, residue_groups=groups
    )
    assert np.allclose(recovered, truth, atol=2e-4)


def test_origin_neighbours_match_genes_and_complete_labels_without_exporting_pgen():
    from vdjtools.model.io import from_germline
    from vdjtools.model import n_region_posterior_batch
    from vdjmatch.match.scoring_audit import origin_neighbour_tables

    germline = pl.DataFrame(
        [
            dict(
                allele=a,
                gene=a.split("*")[0],
                segment=k,
                sequence=s,
                cdr3_anchor=-1,
                functional=True,
                full_germline="",
            )
            for a, k, s in [("TOYV1*01", "V", "TGT"), ("TOYJ1*01", "J", "TTT")]
        ]
    )
    model = from_germline(germline, locus="TOY", ins_max=3, palindrome_max=0)
    queries = pl.DataFrame(
        dict(
            query_id=[0, 1, 2, 3],
            cdr3=["CAF", "CCF", "CDF", "CEF"],
            v=["TOYV1*01"] * 4,
            j=["TOYJ1*01"] * 4,
        )
    )
    membership = pl.DataFrame(dict(query_id=[0, 1, 2, 3], epitope=["E", "E", "F", "G"]))
    labels = pl.concat(
        [membership, pl.DataFrame(dict(query_id=[0, 2], epitope=["TAIL", "TAIL"]))]
    )
    origins = n_region_posterior_batch(
        model, queries["cdr3"].to_list(), threads=1
    ).rename({"row": "query_id"})
    tables = origin_neighbour_tables(
        model, queries, membership, origins, labels, threads=1
    )
    observed = tables["neighbour-counts"].filter(pl.col("arm") == "all")
    same = observed.filter(pl.col("comparison") == "same")
    different = observed.filter(pl.col("comparison") == "different")
    w = (
        origins.filter(pl.col("position") == 1)
        .sort("query_id")["n_probability"]
        .to_numpy()
    )
    assert same["observed_pair_weight"].sum() == pytest.approx(w[0] * w[1])
    assert same["expected_pair_weight"].sum() == pytest.approx(w[0] * w[1])
    e_negative = w[0] * w[3] + w[1] * w[2] + w[1] * w[3]
    assert different.filter(pl.col("epitope") == "E")[
        "observed_pair_weight"
    ].sum() == pytest.approx(e_negative)
    assert different["observed_pair_weight"].sum() == pytest.approx(
        2 * (e_negative + w[2] * w[3])
    )
    assert not tables["neighbour-matrix"]["supported"].any()
    assert (
        different.filter(
            (pl.col("epitope") == "E") & (pl.col("a") == "A") & (pl.col("b") == "D")
        )["observed_pair_weight"].sum()
        == 0
    )
    assert all("pgen" not in t.columns for t in tables.values())

    one_epitope = membership.with_columns(pl.lit("E").alias("epitope"))
    sources = one_epitope.with_columns(
        pl.Series("reference_id", ["S1", "S1", "S2", "S2"])
    )
    heldout = origin_neighbour_tables(
        model,
        queries,
        one_epitope,
        origins,
        one_epitope,
        source_membership=sources,
        threads=1,
    )
    split = heldout["neighbour-source-split"].row(0, named=True)
    assert split["heldout_source"] == "S1"
    assert split["training_junctions"] == split["test_junctions"] == 2
    train_counts = heldout["neighbour-counts"].filter(pl.col("arm") == "train")
    assert train_counts["observed_pair_weight"].sum() == pytest.approx(w[2] * w[3])
    prediction = heldout["neighbour-prediction"]
    assert prediction.height == 6
    assert np.allclose(prediction["test_pair_weight"], w[0] * w[1])
    assert np.isfinite(prediction["pair_gain_bits_per_weighted_event"].to_numpy()).all()


def test_neighbour_completion_uses_all_twenty_exposed_residues():
    from vdjtools.model.io import from_germline
    from vdjtools.model import n_region_posterior_batch
    from vdjmatch.match.scoring_audit import AA, origin_neighbour_tables

    germline = pl.DataFrame([
        dict(allele=a, gene=a.split('*')[0], segment=k, sequence=s,
             cdr3_anchor=-1, functional=True, full_germline='')
        for a, k, s in [('TOYV1*01', 'V', 'TGT'), ('TOYJ1*01', 'J', 'TTT')]
    ])
    model = from_germline(germline, locus='TOY', ins_max=3, palindrome_max=0)
    queries = pl.DataFrame(dict(query_id=list(range(20)), cdr3=['C'+a+'F' for a in AA],
                                v=['TOYV1*01'] * 20, j=['TOYJ1*01'] * 20))
    membership = queries.select('query_id').with_columns(pl.lit('E').alias('epitope'))
    origins = n_region_posterior_batch(model, queries['cdr3'].to_list(), threads=1).rename({'row': 'query_id'})
    tables = origin_neighbour_tables(model, queries, membership, origins, membership, threads=1)
    assert tables['neighbour-counts'].filter(pl.col('comparison') == 'same')['observed_pair_incidences'].sum() == 190
    for _, part in tables['neighbour-matrix'].group_by('weighting', 'estimator'):
        part = part.sort('a', 'b')
        q = part['generator_rate'].to_numpy().reshape(20, 20)
        f = part['stationary_frequency'].to_numpy().reshape(20, 20)[:, 0]
        cost = part['distance_bits'].to_numpy().reshape(20, 20)
        assert part['supported'].all()
        assert np.allclose(q.sum(1), 0) and np.allclose(f @ q, 0)
        assert np.isclose(-f @ np.diag(q), 1)
        assert np.allclose(cost, cost.T) and np.array_equal(np.diag(cost), np.zeros(20))


def test_residue_fit_keeps_small_groups_when_pair_mass_is_unequal():
    from scipy.special import softmax
    from vdjmatch.match.scoring_audit import fit_exchangeability

    rng = np.random.default_rng(93)
    i, j = np.triu_indices(20, 1)
    groups = np.repeat(np.arange(3), 6)
    rho = rng.dirichlet(np.ones(20) * 3, len(groups))
    baseline = rho[:, i] * rho[:, j]
    baseline /= baseline.sum(1, keepdims=True)
    eta = rng.normal(0, .3, (3, 20))
    target = softmax(np.log(baseline) + eta[groups][:, i] + eta[groups][:, j], axis=1)
    observed = target * np.array([1e9, 1., 1e-5])[groups, None]
    _, joint = fit_exchangeability(observed, baseline, balanced=False,
                                  residue_groups=groups, fit_pairs=False)
    _, single = fit_exchangeability(observed[groups == 1], baseline[groups == 1], balanced=False,
        residue_groups=np.zeros(6, dtype=int), fit_pairs=False, prior_mass=1/3)
    assert np.allclose(np.asarray(joint['residue_log_preferences'])[1],
                       np.asarray(single['residue_log_preferences'])[0], atol=2e-3)


def test_all_site_rate_opportunity_retains_zero_neighbour_and_unknown_sites(monkeypatch):
    from vdjtools.model.io import from_germline
    from vdjtools.model import Model,n_region_posterior_batch
    import vdjtools.model.native as native
    from vdjmatch.match.scoring_audit import AA,origin_rate_tables

    germline=pl.DataFrame([dict(allele=a,gene=a.split('*')[0],segment=k,sequence=s,
        cdr3_anchor=-1,functional=True,full_germline='') for a,k,s in
        [('TOYV1*01','V','TGT'),('TOYJ1*01','J','TTT')]])
    model=from_germline(germline,locus='TOY',ins_max=6,palindrome_max={'v_3':0,'j_5':0})
    tables=dict(model.tables)
    for name in ['v_3_del','j_5_del']:
        tables[name]=tables[name].with_columns(p=(pl.col('ndel')==0).cast(pl.Float64))
    tables['vj_ins']=tables['vj_ins'].with_columns(p=(pl.col('length')==6).cast(pl.Float64))
    for name in ['vj_first_nt','vj_dinucl']:
        tables[name]=tables[name].with_columns(p=pl.lit(.25))
    model=Model(manifest=model.manifest,tables=tables,genomic=model.genomic).validate()
    queries=pl.DataFrame(dict(query_id=[0,1,2,3],cdr3=['CAAF','CCAF','CDCF','AAAA'],v=['TOYV1*01']*4,j=['TOYJ1*01']*4))
    primary=queries.select('query_id').with_columns(pl.lit('E').alias('epitope'))
    origins=n_region_posterior_batch(model,queries['cdr3'].to_list(),threads=1).rename({'row':'query_id'})
    degeneracy=dict(zip(AA,[4,2,2,2,2,4,2,3,2,6,1,2,4,2,6,6,4,4,1,2]))
    calls=[]
    real_masked = native.pgen_aa_degenerate_batch
    def masked_oracle(model,allowed,v,j,*,n_regions,threads):
        calls.append(len(allowed))
        result=[]
        for motif,(p,end) in zip(allowed,n_regions):
            assert p in [1,2] and end==p+1 and len(motif[p])==19
            other=2 if p==1 else 1
            result.append(sum(degeneracy[a] for a in motif[p])*degeneracy[motif[other]]/4096)
        import inspect
        if 'n_regions' in inspect.signature(real_masked).parameters:
            actual = real_masked(model, allowed, v=v, j=j, n_regions=n_regions, threads=threads)
            assert np.allclose(actual, result, rtol=1e-12, atol=0)
            return actual
        return result
    monkeypatch.setattr(native,'pgen_aa_degenerate_batch',masked_oracle)
    counts=origin_rate_tables(model,queries,primary,origins,threads=1)['origin-rate-counts']
    assert calls==[6]  # One native operation; fixed/unknown sites never become N queries.
    zero=counts.filter((pl.col('position')==1)&(pl.col('aa')=='D')).row(0,named=True)
    assert zero['N_site_weight']==1 and zero['endpoint_incidences']==0
    assert zero['expected_endpoint_opportunity']==pytest.approx(59/1024)
    assert counts['unknown_origin_sites'].sum()==4
    assert counts['N_endpoint_weight'].sum()==2 and counts['endpoint_incidences'].sum()==2
    assert counts.filter(pl.col('position').is_in([0,3]))['N_site_weight'].drop_nulls().sum()==0


def test_profiled_all_site_rates_recover_factorization_and_keep_zero_event_strata():
    from vdjmatch.match.scoring_audit import AA, origin_rate_models

    rng = np.random.default_rng(72)
    positions = np.repeat(np.arange(1, 4), 20)
    residues = np.tile(np.arange(20), 3)
    target = np.linspace(-.4, .4, 20)[residues] + np.array([-.3, .1, .2])[positions-1]
    baseline = rng.uniform(.5, 2., (3, len(target)))
    observed = baseline * np.exp(target) * np.array([2e5, 3e5, 0.])[:, None]
    counts = pl.DataFrame([dict(epitope=f'E{e}', v='V', j='J', L=6,
        position=int(p), aa=AA[a], expected_endpoint_opportunity=float(baseline[e, c]),
        N_endpoint_weight=float(observed[e, c]))
        for e in range(3) for c, (p, a) in enumerate(zip(positions, residues))])
    tables = origin_rate_models(counts)
    fitted = tables['origin-rate-models'].filter(pl.col('model') == 'factorized').sort('position', 'aa')
    assert np.allclose(fitted['log_relative_rate'], target-target.mean(), atol=2e-4)
    fits = tables['origin-rate-fits']
    assert fits['zero_event_strata'].to_list() == [1, 1, 1]
    assert fits['fit_strata'].to_list() == [2, 2, 2]
    assert 'fit_epitopes' not in fits.columns
    assert fits.filter(pl.col('model') == 'factorized')['working_deviance'][0] < 1e-4
    assert tables['origin-rate-intercepts'].filter(pl.col('epitope') == 'E2')['log_stratum_intensity'].to_list() == [-np.inf]*3
    assert fitted['likelihood_component'].n_unique() == 1


def test_canonical_gap_masks_partition_repeated_residue_pairs():
    from itertools import product
    from seqtree import gapblock

    prior = lambda p, k, m: 0 if 1 <= p <= m-k-1 else gapblock.UNREACHABLE
    for interior_length in range(4):
        for middle in product('AC', repeat=interior_length):
            short = 'C' + ''.join(middle) + 'F'
            for k in range(1, 4):
                represented = {}
                possible = set()
                for p in range(1, len(short)):
                    for insertion in product('AC', repeat=k):
                        long = short[:p] + ''.join(insertion) + short[p:]
                        possible.add(long)
                        if p == 1 or insertion[-1] != short[p-1]:
                            assert long not in represented
                            represented[long] = p
                assert set(represented) == possible
                for long, p in represented.items():
                    assert gapblock.gapblock_score(short, long, gap_open=0, gap_extend=0, gap_prior=prior) == (0, p)


def test_shared_immutable_codon_context_is_not_product_of_aa_frequencies():
    from arda.refbuild.translate import translate
    from vdjmatch.match.scoring_audit import AA

    frequency = np.zeros((2, 20))
    for row, fixed in enumerate(['AT', 'GT']):
        for n in 'ACGT':
            frequency[row, AA.index(translate(fixed+n))] += .25
    marginal = frequency.mean(0)
    joint = frequency.T @ frequency / 2
    covariance = (frequency-marginal).T @ (frequency-marginal) / 2
    assert np.allclose(joint-np.outer(marginal, marginal), covariance)
    i, v = AA.index('I'), AA.index('V')
    assert joint[i, v] == 0 and marginal[i]*marginal[v] == .1875
    assert joint[i, i] == .28125 and marginal[i]**2 == .140625


def test_joint_rate_preconditioning_matches_scalar_stationarity_solution():
    from scipy.special import softmax
    from vdjmatch.match.scoring_audit import fit_exchangeability

    baseline = softmax(np.linspace(-24, 0, 60))
    observed = softmax(np.log(baseline) + np.sin(np.arange(60))) * 50000
    from scipy.optimize import brentq
    # The neutral prior acts on rate effects, not on baseline-weighted probabilities.
    root = brentq(lambda z: np.sum((observed+1/60)/(1+z*baseline))-1,
                  0, (observed.sum()+1)/baseline.min())
    exact = np.log((observed+1/60)/(1+root*baseline))
    exact -= exact.mean()
    fitted, report = fit_exchangeability(observed[None], baseline[None], balanced=False, precondition=True)
    assert np.allclose(fitted, exact, atol=2e-3)
    assert report['max_scaled_gradient'] < 1e-6


def test_local_origin_pair_channels_partition_without_compulsory_conservation():
    from vdjmatch.match.scoring_audit import origin_class_pair_counts
    encoded=np.array([[0,15],[1,15]],dtype=np.uint8)
    classes=pl.DataFrame(dict(row=[0,0,1,1,0,1],position=[0,0,0,0,1,1],
        pattern=[0,1,0,1,2,2],probability=[.4,.6,.7,.3,1.,1.],
        exchangeable=[True,True,True,True,False,False]))
    channels=origin_class_pair_counts(encoded,classes)
    assert set(channels)=={'same_exchangeable','same_zero_alternative','different_class'}
    assert channels['same_exchangeable'][0,1]==pytest.approx(.46)
    assert channels['same_exchangeable'].sum()==pytest.approx(.92)
    assert channels['same_zero_alternative'][15,15]==pytest.approx(2.)
    assert channels['different_class'][0,1]==pytest.approx(.54)
    total=sum(channels.values())
    assert np.allclose(total,column_pair_counts(encoded,np.arange(2),np.ones_like(encoded))[0])
    shuffled=origin_class_pair_counts(encoded,classes.reverse())
    assert all(np.allclose(channels[k],shuffled[k]) for k in channels)
    with pytest.raises(ValueError,match='partition'):
        origin_class_pair_counts(encoded,classes.with_columns(pl.col('probability')*.5))
    with pytest.raises(ValueError,match='unique'):
        origin_class_pair_counts(encoded,pl.concat([classes,classes.head(1)]))


def test_mixed_block_owner_reuses_NNN_counts_and_keeps_unknowns():
    from vdjmatch.match.scoring_audit import mixed_origin_block_matrices
    clones,origins=_inputs()
    classes=origins.rename({'n_probability':'probability'}).with_columns(
        pl.when(pl.col('probability').is_not_null()).then(0).otherwise(-1).alias('pattern'),
        pl.lit(True).alias('exchangeable'))
    # Complete supported posteriors, retaining the original unknown receptor.
    classes=classes.with_columns(pl.when(pl.col('probability').is_not_null()).then(1.).otherwise(None).alias('probability'))
    expected=origin_block_matrices(clones,classes.rename({'probability':'n_probability'}))[0]
    tables=mixed_origin_block_matrices(clones,classes)
    got=tables['mixed-matrix'].filter(pl.col('channel')=='same_exchangeable').drop('channel')
    assert got.equals(expected)
    assert tables['mixed-support'].filter(~pl.col('origin_supported'))['query_id'].to_list()==[4]
    empty_channels=tables['mixed-matrix'].filter(pl.col('channel')!='same_exchangeable')
    assert empty_channels['pair_weight'].sum()==0
