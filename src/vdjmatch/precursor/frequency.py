"""From a set of `Pgen` values to a precursor frequency, a cell count, and a response probability.

Everything above this module returns a *mass* under the recombination model. This is where that
becomes a number about a person: what fraction of their T cells can see the epitope, how many cells
that is, and how likely at least ``k`` of them exist.
"""
from __future__ import annotations

import math

from .mass import ALPHA_PER_EDIT, shell_profile

#: ALICE's thymic-selection factor `Q` (Pogorelyy et al., *PLoS Biol* 2019;17:e3000314), the
#: constant by which post-selection frequency exceeds the raw generation probability. Calibrated on
#: **TRB only**, against OLGA's scale, and shipped here as a named citable default rather than a
#: measured property of this library. The matched comparison against the model-free event ratio
#: puts the same factor at a median 14.8x, so the two are the same order and are not the same
#: number; pass your own ``q`` when you have fitted one for your chain and cohort.
ALICE_Q = 9.41

#: Order-of-magnitude total T-cell count for one adult human. Used only to turn a frequency into a
#: count; it is not a measured constant of this library and every reported count inherits its error.
N_T_CELLS = 1e11


def precursor_frequency(model, junctions, r: int = 1, alpha: float = ALPHA_PER_EDIT,
                        q: float = 1.0, n_cells: float = N_T_CELLS, compartment: float = 1.0,
                        n_eff: float | None = None, threads: int = 0) -> dict:
    """`F(e)` as a frequency, with the cell count and precursor probabilities that follow.

    ::

        F_hat(e) = q * sum_{k=0..r} alpha^k * mass(shell k)
        cells    = n_cells * compartment * F_hat(e)
        lambda   = n_eff * F_hat(e)        P(>=m precursors) = 1 - Poisson_cdf(m-1; lambda)

    The mass term is :func:`shell_profile`'s ``retained``: the union of the observed junctions'
    Hamming-``r`` balls, resolved into min-distance shells and down-weighted by the measured cognacy
    retention ``alpha`` per edit, so a neighbour ``k`` substitutions away contributes ``alpha^k`` of
    its mass rather than all of it. Shells partition the union, so this neither double-counts the
    overlap between two cognate junctions' balls nor assumes every neighbour is still cognate.

    ``q`` is the **selection constant** that carries a generation probability to a post-selection
    repertoire frequency. It defaults to ``1.0``, i.e. **uncalibrated** — the returned ``F`` is then
    a raw model mass and only its *ranking* is meaningful. Pass :data:`ALICE_Q` for the published
    TRB value, or a constant fitted against your own cohort. Whatever you pass is echoed back in the
    result so a reported number always carries its calibration.

    ``compartment`` is the fraction of ``n_cells`` the epitope's restriction actually addresses —
    for an MHC-I epitope, the CD8 fraction, times the naive fraction if you are after a *naive*
    precursor count rather than a realised one. Left at ``1.0`` the ``cells`` field is "cells in the
    whole T-cell pool", which is rarely the quantity a tetramer experiment measured.

    ``n_eff`` is the number of **independent rearrangements** the probabilities are about, and it is
    not the same as ``n_cells``: set it to a donor's distinct-clonotype count to ask "does this
    person have a precursor at all", or to a sample's sequencing depth to ask "will I see one in
    this tube". Left as ``None`` the Poisson fields are omitted rather than guessed.

    Returns ``{"F", "q", "alpha", "r", "union", "retained", "overlap", "n_seqs", "shells",
    "cells", "n_cells", "compartment", "lambda", "p_ge_1", "p_ge_10"}``, with the last three
    ``None`` when ``n_eff`` is not given.

    **What `F` is not.** It is the probability that a precursor *exists*, not that a response is
    *mounted*: a detectable response needs more than one cell, which is what the ``p_ge_*`` fields
    are for, and it needs priming, help and an absence of tolerance, which nothing here models.
    """
    prof = shell_profile(model, junctions, r=r, alpha=alpha, threads=threads)
    f = q * prof["retained"]
    out = {"F": f, "q": q, "alpha": alpha, "r": r,
           "union": prof["union"], "retained": prof["retained"], "overlap": prof["overlap"],
           "n_seqs": prof["n_seqs"], "shells": prof["shells"],
           "cells": expected_cells(f, n_cells=n_cells, compartment=compartment),
           "n_cells": n_cells, "compartment": compartment,
           "lambda": None, "p_ge_1": None, "p_ge_10": None}
    if n_eff is not None:
        lam = float(n_eff) * f
        out.update({"lambda": lam, "p_ge_1": p_at_least(f, n_eff, 1),
                    "p_ge_10": p_at_least(f, n_eff, 10)})
    return out


def expected_cells(f: float, n_cells: float = N_T_CELLS, compartment: float = 1.0) -> float:
    """Expected number of epitope-specific cells: ``n_cells * compartment * f``.

    ``compartment`` restricts ``n_cells`` to the subset the epitope can address (the CD8 fraction
    for MHC-I, times the naive fraction for a naive precursor count). The result inherits the error
    of both ``f`` and the assumed pool size, so quote it as an order of magnitude.
    """
    return float(n_cells) * float(compartment) * float(f)


def p_at_least(f: float, n_eff: float, k: int = 1) -> float:
    """``P(X >= k)`` for ``X ~ Poisson(n_eff * f)`` — the probability that ``k`` precursors exist.

    ``F(e)`` on its own answers "is there at least one?", i.e. ``k = 1``. A detectable assay or
    clinical response needs more, and once ``n_eff * F`` is of order one the two stop being a
    monotone reparametrisation of each other — which is exactly the regime a rare neoepitope sits
    in. That is why ``k`` is an argument rather than a fixed 1.

    ``n_eff`` is a count of independent rearrangements; see :func:`precursor_frequency`.
    """
    if k < 1:
        raise ValueError("k must be >= 1")
    lam = float(n_eff) * float(f)
    if lam <= 0:
        return 0.0
    # 1 - CDF(k-1); summing forward is stable for the small k this is used with
    term, cdf = math.exp(-lam), 0.0
    for i in range(k):
        cdf += term
        term *= lam / (i + 1)
    return max(0.0, 1.0 - cdf)


def paired_frequency(f_alpha: float, f_beta: float, kappa: float = 1.0) -> float:
    """Frequency of cells whose **both** chains are cognate: ``kappa * f_alpha * f_beta``.

    Alpha and beta pair essentially without constraint across the repertoire, so unconditionally
    ``P(A, B) = P(A) P(B)``. **Given an epitope they do not**: a cognate beta works with a restricted
    set of alphas, not with any alpha drawn from the cognate alpha set. ``kappa`` carries that
    departure — ``kappa = 1`` asserts fully permissive pairing within the cognate sets and is an
    upper bound; ``kappa < 1`` is restricted pairing. It is measurable on VDJdb's paired records and
    is not assumed here, which is why it defaults to the bound and must be passed to claim anything
    tighter.

    Both inputs must be the same kind of quantity — two calibrated ``F`` values, or two raw masses,
    never one of each.
    """
    return float(kappa) * float(f_alpha) * float(f_beta)
