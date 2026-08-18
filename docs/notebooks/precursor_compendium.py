"""Does the recombination model predict a measured naive precursor frequency?

A marimo notebook. Run it with::

    marimo edit docs/notebooks/precursor_compendium.py

Nothing here recomputes a `Pgen`: both tables come from the HuggingFace mirror, so the notebook
opens in a second on any machine with `vdjmatch[precursor]` installed.
"""

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # Precursor frequency: prediction against measurement

        `F(e)` estimates the fraction of a naive repertoire that can see an epitope. Other
        laboratories have measured that quantity directly — tetramer and dextramer enrichment of
        unexposed donors, umbilical cord blood, naive mice — so this is a comparison against the
        estimand itself rather than against a proxy.

        Two tables, kept deliberately apart:

        | | provenance |
        |---|---|
        | `load_compendium()` | **experimental** — what a laboratory counted |
        | `load_estimates()` | **derived** — a mass under a generative recombination model |

        Both live at `isalgo/airr_benchmark` under `vdjmatch/precursor_freq/`.
        """
    )
    return


@app.cell
def _():
    import marimo as mo
    import numpy as np
    import polars as pl

    from vdjmatch.precursor import load_compendium, load_estimates

    comp = load_compendium()      # measured, naive, one row per published record
    est = load_estimates()        # predicted, one row per (species, model, chain, epitope)
    return comp, est, load_compendium, load_estimates, mo, np, pl


@app.cell(hide_code=True)
def _(comp, est, mo):
    mo.md(
        f"""
        ## What is in the tables

        **{comp.height} measured rows**, {comp["study"].n_unique()} studies, tracing to
        **{comp["lab"].n_unique()} originating laboratories** — one of the studies is a review,
        and its rows carry the laboratory that made them in `lab`, so the corpus is countable by
        laboratory rather than by paper.

        **{est.height} predictions**, over {est["epitope_seq"].n_unique()} epitopes, scored under
        both bundled model sets where the species allows it.
        """
    )
    return


@app.cell
def _(comp, pl):
    # The measured axis, by species and subset. `freq_per_1e6_naive` is the harmonized column:
    # cells per million of the *naive* compartment, not per total CD8.
    coverage = (
        comp.group_by(["species", "subset"])
        .agg(
            pl.col("lab").n_unique().alias("labs"),
            pl.col("epitope_seq").filter(pl.col("epitope_seq") != "").n_unique().alias("epitopes"),
            pl.len().alias("rows"),
            pl.col("freq_per_1e6_naive").median().round(2).alias("median_per_1e6"),
        )
        .sort(["species", "subset"])
    )
    coverage
    return (coverage,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Agreement, per epitope

        Both axes in cells per million of the naive compartment, so the diagonal is the ideal. The
        estimator columns are repertoire *fractions*, hence the `1e6`.

        Pick the arm and the estimator:
        """
    )
    return


@app.cell
def _(est, mo):
    arm = mo.ui.dropdown(
        options={f"{s} / {m} / {c}": (s, m, c)
                 for s, m, c in est.select(["species", "source", "chain"]).unique().sort(
                     ["species", "source", "chain"]).iter_rows()},
        value="human / arda / TRB",
        label="arm",
    )
    estimator = mo.ui.dropdown(
        options=["observed_mass", "union", "retained", "occupancy"],
        value="occupancy", label="estimator",
    )
    mo.hstack([arm, estimator], justify="start", gap=2)
    return arm, estimator


@app.cell
def _(arm, est, estimator, np, pl):
    sp, md, ch = arm.value
    arm_df = est.filter(
        (pl.col("species") == sp) & (pl.col("source") == md) & (pl.col("chain") == ch)
    ).filter(pl.col(estimator.value).is_not_null() & (pl.col(estimator.value) > 0))

    pred = arm_df[estimator.value].to_numpy() * 1e6
    meas = arm_df["measured_per_1e6"].to_numpy()
    ratio = pred / meas
    return arm_df, ch, md, meas, pred, ratio, sp


@app.cell
def _(arm_df, estimator, meas, mo, np, pred, ratio):
    from scipy.stats import spearmanr

    rho = spearmanr(pred, meas)
    # A stable offset leaves rankings intact and is what a depth factor absorbs, so the *spread*
    # matters more than the median -- a median that agrees is not an epitope that agrees.
    mo.md(
        f"""
        **n = {len(pred)} epitopes** &nbsp;·&nbsp; Spearman ρ = **{rho.statistic:.3f}**
        (p = {rho.pvalue:.3g})

        median predicted/measured = **{np.median(ratio):.3g}×** &nbsp;·&nbsp;
        spread (max/min) = **{ratio.max() / ratio.min():.0f}×**

        The median is the calibration; the spread is why `{estimator.value}` supports a **ranking**
        and not an absolute frequency for a named epitope.
        """
    )
    return rho, spearmanr


@app.cell
def _(arm_df, ch, estimator, md, meas, pred, sp):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    lo = min(pred.min(), meas.min()) * 0.4
    hi = max(pred.max(), meas.max()) * 2.5
    ax.plot([lo, hi], [lo, hi], ls="--", lw=1, color="#666666", label="y = x")
    ax.scatter(pred, meas, s=34, color="#1b9e77" if ch == "TRA" else "#d95f02",
               edgecolor="white", linewidth=0.5, zorder=3)
    for x, y, e in zip(pred, meas, arm_df["epitope_seq"]):
        ax.annotate(e, (x, y), fontsize=6, alpha=0.75,
                    xytext=(3, 3), textcoords="offset points")
    ax.set(xscale="log", yscale="log", xlim=(lo, hi), ylim=(lo, hi),
           xlabel=f"predicted {estimator.value}  (per 10⁶ naive)",
           ylabel="measured naive frequency  (per 10⁶ naive)",
           title=f"{sp} · {md} · {ch}")
    ax.legend(frameon=False, loc="upper left", fontsize=8)
    ax.grid(alpha=0.25, which="major")
    fig.tight_layout()
    fig
    return ax, e, fig, hi, lo, plt, x, y


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## What else the `Pgen` distribution and the junction count give you

        The estimators above are one summary of a cognate set. The same two inputs — the `Pgen`
        spectrum and how many junctions were catalogued — answer several other questions, and they
        are not interchangeable.
        """
    )
    return


@app.cell
def _(arm_df, pl):
    # The share of the one-substitution ball's mass that the catalogued junctions already carry.
    # `union` is exact -- it never enumerates -- so this is available for every epitope at once.
    seen = (
        arm_df.select(
            "epitope_seq",
            pl.col("n_junctions").alias("catalogued"),
            (pl.col("observed_mass") / pl.col("union")).round(4).alias("seen_mass_frac"),
            (pl.col("union") * 1e6).round(3).alias("union_per_1e6"),
            (pl.col("observed_mass") * 1e6).round(3).alias("observed_per_1e6"),
        )
        .sort("seen_mass_frac")
    )
    seen
    return (seen,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        The **count** of uncatalogued junctions is a separate computation, because it needs the ball
        enumerated rather than only its mass. `union_mass(..., count_members=True)` returns it:
        """
    )
    return


@app.cell
def _(chain_pick, epitope_census, mo):
    mo.hstack([epitope_census, chain_pick], justify="start", gap=2)
    return


@app.cell
def _(mo):
    epitope_census = mo.ui.text(value="GILGFVFTL", label="epitope for the ball census")
    return (epitope_census,)


@app.cell
def _(chain_pick, epitope_census, pl):
    from vdjmatch import db as _db
    from vdjmatch.precursor import (
        check_junctions as _check,
        load_model as _load_model,
        observed_mass as _observed_mass,
        union_mass as _union_mass,
    )

    _frame = _db.load(asset="slim", species="HomoSapiens", gene=chain_pick.value)
    _js, _ = _check(_frame.filter(pl.col("epitope") == epitope_census.value)["cdr3"]
                    .unique().to_list())
    _js = sorted(set(_js))
    _model = _load_model(chain_pick.value, source="arda", organism="human")
    u = _union_mass(_model, _js, r=1, count_members=True)
    obs = _observed_mass(_model, _js)
    return _check, _db, _frame, _js, _load_model, _model, _observed_mass, _union_mass, obs, u


@app.cell(hide_code=True)
def _(chain_pick, epitope_census, mo, obs, u):
    _uncat = None if u["n_union"] is None else u["n_union"] - u["n_seqs"]
    mo.md(
        f"""
        **`{epitope_census.value}` · {chain_pick.value}** — {u["n_seqs"]:,} catalogued junctions,
        {"" if _uncat is None else f"{u['n_union']:,} in the one-substitution ball, so "}
        **{"too many to enumerate" if _uncat is None else f"{_uncat:,} uncatalogued"}** candidate
        neighbours.

        - ball mass `union` = **{u["union"]:.3g}** &nbsp;·&nbsp; catalogued mass
          `observed_mass` = **{obs:.3g}** → the catalogued junctions carry
          **{obs / u["union"]:.1%}** of the neighbourhood's mass.
        - naive per-sequence sum would have been {u["naive_sum"]:.3g}, so overlap between the balls
          double-counts **{u["overlap"]:.1%}** of it. Cognate TCRs are near-duplicates by
          construction, so their balls overlap and the *union* is the correct object — the sum is
          not.

        Two different "unseen" questions, and conflating them is the classic error:

        - **Unseen count** is a *census* of the ball: finite and exact given the ball. A
          Horvitz–Thompson richness extrapolation is not — `1/π` diverges as `π → 0` and the first
          real run returned 10¹¹–10¹⁴ unseen junctions, which is not a number.
        - **Unseen mass** *does* converge, because the Horvitz–Thompson weight `p/π → 1/N` as
          `p → 0`. That is what makes the correction work without guessing the shape of the tail.

        The missing count can be enormous while the missing mass stays small, because the unseen
        members are systematically the low-`Pgen` ones. That asymmetry is why the library ships both
        and marks a richness number unreliable when the rarest observed junction was captured too
        seldom for it to mean anything.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## The spread inside one cognate set

        There is no "typical" member of a cognate set, which is why no estimator here uses a mean.
        Scoring one epitope's junctions directly shows the range: generation probabilities inside a
        single specificity group span a median 7.8 decades on TRA and 8.9 on TRB.
        """
    )
    return


@app.cell
def _(mo):
    epitope = mo.ui.text(value="NLVPMVATV", label="epitope (VDJdb)")
    chain_pick = mo.ui.dropdown(options=["TRB", "TRA"], value="TRB", label="chain")
    mo.hstack([epitope, chain_pick], justify="start", gap=2)
    return chain_pick, epitope


@app.cell
def _(chain_pick, epitope, np, pl):
    from vdjmatch import db
    from vdjmatch.precursor import check_junctions, load_model, pgen

    frame = db.load(asset="slim", species="HomoSapiens", gene=chain_pick.value)
    js = (frame.filter(pl.col("epitope") == epitope.value)["cdr3"].unique().to_list())
    js, dropped = check_junctions(js)
    model = load_model(chain_pick.value, source="arda", organism="human")
    p = np.asarray(pgen(model, sorted(set(js))), dtype=float)
    p_nz = p[p > 0]
    return check_junctions, db, dropped, frame, js, load_model, model, p, p_nz, pgen


@app.cell
def _(chain_pick, epitope, js, np, p, p_nz, plt):
    fig2, ax2 = plt.subplots(figsize=(5.6, 3.4))
    ax2.hist(np.log10(p_nz), bins=40, color="#7570b3", edgecolor="white", linewidth=0.4)
    ax2.set(xlabel="log₁₀ Pgen", ylabel="junctions",
            title=f"{epitope.value} · {chain_pick.value} · n = {len(js):,} junctions, "
                  f"{np.ptp(np.log10(p_nz)):.1f} decades")
    ax2.grid(alpha=0.25)
    fig2.tight_layout()
    fig2
    return ax2, fig2


@app.cell(hide_code=True)
def _(js, mo, np, p, p_nz):
    order = np.sort(p_nz)[::-1]
    mo.md(
        f"""
        **{len(p_nz):,} junctions score a non-zero `Pgen`**, spanning
        **{np.ptp(np.log10(p_nz)):.1f} decades**; {int((p <= 0).sum())} score exactly zero.

        The summed mass is carried by a small minority: the top 10% of junctions hold
        **{order[: max(1, len(order) // 10)].sum() / order.sum():.1%}** of it. That concentration is
        the `Σπ²` term that makes a naive set *total* the wrong functional — under the size bias
        governing which receptors get catalogued, a set total has expectation `n·(Σπ²)/F`, which is
        proportional to how many receptors were catalogued and *inversely* proportional to the
        quantity it is meant to estimate.

        A junction that scores exactly zero contributes nothing to any mass and raises no error, so
        it is reported rather than dropped silently — under the `olga` model set 5.7% of human TRA
        junctions land there.
        """
    )
    return (order,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## From a mass to a number of cells

        `F(e)` is the probability a precursor *exists*. A detectable response needs more than one
        cell, and once `n_eff · F` is of order one the two stop being a monotone reparametrisation
        — which is the regime a rare neoepitope sits in.
        """
    )
    return


@app.cell
def _(mo):
    f_slider = mo.ui.slider(start=-7.0, stop=-3.0, step=0.1, value=-5.0,
                            label="log₁₀ F(e)", show_value=True)
    k_slider = mo.ui.slider(start=1, stop=20, step=1, value=1,
                            label="k (precursors needed)", show_value=True)
    mo.hstack([f_slider, k_slider], justify="start", gap=2)
    return f_slider, k_slider


@app.cell
def _(f_slider, k_slider, mo):
    from vdjmatch.precursor import N_T_CELLS, expected_cells, p_at_least

    f = 10 ** f_slider.value
    n_eff = 0.02 * N_T_CELLS          # naive CD8 is a fraction of the 1e11 total
    cells = expected_cells(f, n_cells=n_eff)
    prob = p_at_least(f, n_eff=n_eff, k=k_slider.value)
    mo.md(
        f"""
        At `F = {f:.3g}` and an effective naive pool of `{n_eff:.2g}` cells:

        - expected precursors: **{cells:,.0f}**
        - P(at least {k_slider.value}) = **{prob:.4g}**

        Everything here is uncalibrated: `q` defaults to 1, so `F` is a raw mass whose *ranking*, not
        scale, is the supported claim. The measured per-chain depth factors are
        `SELECTION_BY_CHAIN` (TRB 4.62, TRA 1.07) and are deliberately not applied by default.
        """
    )
    return N_T_CELLS, cells, expected_cells, f, n_eff, p_at_least, prob


if __name__ == "__main__":
    app.run()
