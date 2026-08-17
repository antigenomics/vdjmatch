"""What the specificity database never saw: the missing mass, and the missing count.

Both come from one fitted capture curve. The textbook unseen-species estimators do not apply here
and the reason is structural, not incidental — see :func:`coverage_corrected_mass`.
"""
from __future__ import annotations

import math

from .model import pgen

#: Smallest inclusion probability among the *observed* junctions at which the Horvitz-Thompson
#: **richness** is still worth quoting. Mass and richness behave differently in the tail: the mass
#: weight ``p/pi -> 1/N`` converges as ``p -> 0``, but the richness weight ``1/pi -> 1/(N p)``
#: diverges, so a set whose rarest observed member was barely captured extrapolates to an
#: arbitrarily large unseen count. That is not a defect of this implementation — species richness
#: is simply not estimable from a sample without an assumption about the shape of the tail, which
#: is exactly the assumption this estimator was built to avoid making. Below this floor
#: ``richness_reliable`` is False and the count should be read as "many", not as a number; the
#: finite alternative is the ball census, ``union_mass(...)["n_union"] - n_observed``.
MIN_INCLUSION_FOR_RICHNESS = 0.01


def _fit_capture_rate(u, m, n):
    """MLE of ``theta`` in ``p_i = 1 - exp(-theta * u_i)`` under a zero-truncated Binomial(n, p_i).

    ``u`` are Pgen values rescaled to unit median so the search bracket is scale-free. Returns
    ``(theta, hit_boundary)``. Grid then golden-section: the likelihood is unimodal in practice but
    not provably so, and a coarse grid makes the refinement immune to a bad initial bracket.
    """
    import numpy as np

    u = np.asarray(u, dtype=float)
    m = np.asarray(m, dtype=float)

    def loglik(x):
        t = np.exp(x) * u
        # log(1 - exp(-t)) via -expm1 keeps precision for small t
        return float(np.sum(m * np.log(-np.expm1(-t)) - (n - m) * t - np.log(-np.expm1(-n * t))))

    grid = np.linspace(math.log(1e-8), math.log(1e8), 321)
    vals = np.array([loglik(x) for x in grid])
    k = int(np.argmax(vals))
    if k == 0 or k == len(grid) - 1:
        return float(np.exp(grid[k])), True

    lo, hi = grid[k - 1], grid[k + 1]
    phi = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = hi - phi * (hi - lo), lo + phi * (hi - lo)
    fa, fb = loglik(a), loglik(b)
    for _ in range(60):
        if fa < fb:
            lo, a, fa = a, b, fb
            b = lo + phi * (hi - lo)
            fb = loglik(b)
        else:
            hi, b, fb = b, a, fa
            a = hi - phi * (hi - lo)
            fa = loglik(a)
    return float(math.exp((lo + hi) / 2.0)), False


def coverage_corrected_mass(model, junctions, multiplicity, n_units: int | None = None,
                            threads: int = 0) -> dict:
    """`F(e)` with the size-biased sampling deficit put back — the estimator, not the bound.

    ``multiplicity[i]`` is the number of **independent capture units** (donors, or failing that
    studies) that **re-reported** ``junctions[i]`` for this epitope. Take it from donor/study
    counts, never from VDJdb record counts: rows duplicate the same donor's TCR across curation
    passes, so a record count is not a recapture and inflating it collapses the correction.
    ``n_units`` is the total number of capture units for the epitope (defaults to
    ``max(multiplicity)``).

    .. warning::

       ``multiplicity`` counts **recaptures of one object**, not objects. :func:`event_ratio` uses
       the opposite convention on the same-looking numbers: there each donor carrying a junction is
       a **separate recombination event**, a distinct object. Both are coherent, they are not the
       same quantity, and passing one function's counts to the other is silently wrong rather than
       an error. If you are counting how many independent rearrangements exist, you want
       :func:`event_ratio`; if you are inferring what the sampling never saw, you want this.

    **Model.** Each cognate junction is captured by one unit with probability
    ``p_i = 1 - exp(-theta * pi_i)`` — Poisson sampling at a rate proportional to `Pgen`, which is
    exactly the size-biasing that makes ``observed_mass`` biased. ``theta`` is fitted by maximising
    the zero-truncated Binomial(``n_units``, ``p_i``) likelihood over the *observed* junctions
    (unobserved ones contribute nothing, hence the truncation). The corrected mass is then
    Horvitz–Thompson:

    ``F_hat = sum_i pi_i / (1 - (1 - p_i)^n)``

    which is unbiased for the total over the **whole** cognate set, unobserved members included,
    because each observed member is upweighted by its own inclusion probability. The same
    inclusion probabilities give the richness, ``S_hat = sum_i 1 / (1 - (1 - p_i)^n)``, from which
    :func:`unseen_junctions` reads off how many members were never catalogued and how rare they are.

    Returns ``{"observed", "corrected", "coverage", "theta", "gt_coverage", "f1", "n_units",
    "n_seqs", "n_zero", "n_total", "n_unseen", "unseen_mass", "mean_unseen_pgen", "degenerate",
    "reason"}``; ``coverage = observed / corrected``.

    **Why not Good–Turing.** The textbook flat coverage ``1 - f1/N`` is *known-bad* on TCR data —
    Laydon et al., *PLoS Comput Biol* 2014;10:e1003646 (PMID 24945836) measure 61.7% median error
    on real TCR abundance data, against 43.8% for Chao1bc and 42.8% for ACE, because the
    capture-probability distribution is far too heterogeneous for the uniform-multinomial
    assumption behind it. Here the inclusion probability is **known** rather than inferred —
    `Pgen` *is* the sampling probability — and the consequential behaviour is at the rare end:
    ``p/(1 - exp(-N*p)) -> 1/N`` as ``p -> 0``, a constant. Every rare junction contributes the
    same weight however rare it is, so the estimator never has to guess the *shape* of the tail it
    cannot see, which is precisely the guess that sinks Good–Turing. The flat number is still
    returned as ``gt_coverage`` so the two can be compared, but it is a diagnostic, not the estimate.

    **Richness and mass are not interchangeable.** Because the unseen junctions are systematically
    the low-`Pgen` ones, the missing *count* can be enormous while the missing *mass* stays small.
    ``F(e)`` is made of mass; ``n_unseen`` answers a different question and should be quoted as one.

    **What it assumes.** (1) Capture is independent across units — false where two studies share
    donors or one study's TCRs were re-curated into another, which biases ``theta`` upward and the
    correction downward. (2) The capture curve is the one-parameter saturating form above; real
    capture also depends on assay sensitivity and HLA typing of the cohort, which this folds into a
    single ``theta``. (3) `Pgen` is the right size variable — it is the one the size-bias argument
    names, but clonal expansion in the source samples adds a second, unmodelled one. (4) The cognate
    set is fixed; a junction that is cognate in one donor and not another breaks the frame entirely.

    **Where it breaks, loudly.** When **every junction is a singleton** there is no recapture
    information at all, the likelihood is maximised as ``theta -> 0`` and the Horvitz–Thompson sum
    diverges: the function then sets ``degenerate=True``, names the reason and returns the
    ``observed`` bound unchanged rather than an infinity or a ZeroDivisionError. The same happens
    with fewer than two capture units, and if the fit hits the search boundary.

    **Do not compose this with the neighbourhood route.** :func:`ball_mass` puts back mass that lies
    *near* what was seen, and most of what the sampling missed *is* near what it saw, because
    cognate junctions are near-duplicates by construction. Applying both counts the same missing
    mass twice.
    """
    seqs = list(junctions)
    mult = [int(x) for x in multiplicity]
    if len(mult) != len(seqs):
        raise ValueError(f"multiplicity has {len(mult)} entries for {len(seqs)} junctions")
    if any(x < 1 for x in mult):
        raise ValueError("multiplicity entries must be >= 1: an unobserved junction is not a member")
    n = int(n_units) if n_units is not None else (max(mult) if mult else 0)
    if any(x > n for x in mult):
        raise ValueError(f"multiplicity exceeds n_units={n}: a junction cannot be captured twice "
                         "by the same unit")

    pis = pgen(model, seqs, threads=threads)
    observed = float(sum(pis))
    f1 = sum(1 for x in mult if x == 1)
    captures = sum(mult)
    out = {"observed": observed, "corrected": observed, "coverage": 1.0, "theta": None,
           "gt_coverage": (1.0 - f1 / captures) if captures else 0.0,
           "f1": f1, "n_units": n, "n_seqs": len(seqs), "n_zero": sum(1 for p in pis if p <= 0),
           "n_total": float(len(seqs)), "n_unseen": 0.0, "unseen_mass": 0.0,
           "mean_unseen_pgen": None, "min_inclusion": None, "richness_reliable": False,
           "degenerate": True, "reason": None}

    if not seqs:
        out["reason"] = "no junctions"
        return out
    if n < 2:
        out["reason"] = "fewer than two capture units: no recapture information"
        return out
    if f1 == len(mult):
        out["reason"] = ("every junction is a singleton: the capture curve is unidentified "
                         "(theta -> 0, the Horvitz-Thompson sum diverges)")
        return out

    keep = [i for i, p in enumerate(pis) if p > 0]
    if not keep:
        out["reason"] = "every junction has Pgen 0 under this model"
        return out
    kept = [pis[i] for i in keep]
    med = sorted(kept)[len(kept) // 2]
    u = [p / med for p in kept]
    theta, boundary = _fit_capture_rate(u, [mult[i] for i in keep], n)
    if boundary:
        out["reason"] = "the capture-rate fit hit the search boundary; theta is not identified"
        return out

    corrected = 0.0
    richness = 0.0
    min_incl = 1.0
    for p, ui in zip(kept, u):
        incl = -math.expm1(-n * theta * ui)
        if incl <= 0:                                # pragma: no cover - guarded by `boundary`
            out["reason"] = "an inclusion probability underflowed to 0"
            return out
        corrected += p / incl
        richness += 1.0 / incl
        min_incl = min(min_incl, incl)
    n_unseen = max(0.0, richness - len(keep))
    unseen_mass = max(0.0, corrected - observed)
    out.update(corrected=corrected, coverage=observed / corrected if corrected > 0 else 1.0,
               theta=theta / med, n_total=richness, n_unseen=n_unseen, unseen_mass=unseen_mass,
               mean_unseen_pgen=(unseen_mass / n_unseen) if n_unseen > 0 else None,
               min_inclusion=min_incl, richness_reliable=min_incl >= MIN_INCLUSION_FOR_RICHNESS,
               degenerate=False)
    return out


def unseen_junctions(model, junctions, multiplicity, n_units: int | None = None,
                     threads: int = 0) -> dict:
    """How many cognate junctions were never catalogued, and how rare each one is.

    The richness counterpart of :func:`coverage_corrected_mass`, read off the same fitted capture
    curve: ``S_hat = sum_i 1/pi_i`` is the Horvitz–Thompson estimate of the cognate set's true size,
    so ``n_unseen = S_hat - n_observed`` and the mass they carry is ``corrected - observed``.

    Returns ``{"n_observed", "n_total", "n_unseen", "observed_mass", "unseen_mass",
    "mean_unseen_pgen", "mean_observed_pgen", "rarity_ratio", "min_inclusion",
    "richness_reliable", "degenerate", "reason"}``. ``rarity_ratio`` is how many times rarer an
    average unseen junction is than an average observed one — the number that makes "the missing
    count is huge, the missing mass is not" concrete.

    **Read the mass; treat the count as an order of magnitude at best.** ``unseen_mass`` converges
    (the weight ``p/pi`` tends to a constant as ``p -> 0``) but ``n_unseen`` does not (``1/pi``
    diverges), so a group whose rarest observed junction was barely captured extrapolates to an
    arbitrarily large count. ``richness_reliable`` is False when that is happening — see
    :data:`MIN_INCLUSION_FOR_RICHNESS`. For a finite census of uncatalogued candidates, use the ball
    instead: ``union_mass(model, junctions, r)["n_union"] - n_observed``.

    Degenerates exactly where :func:`coverage_corrected_mass` does and for the same reasons; it
    reports that it cannot answer rather than returning a number.
    """
    c = coverage_corrected_mass(model, junctions, multiplicity, n_units=n_units, threads=threads)
    n_obs = c["n_seqs"] - c["n_zero"]
    mean_obs = (c["observed"] / n_obs) if n_obs else None
    mean_unseen = c["mean_unseen_pgen"]
    return {"n_observed": n_obs, "n_total": c["n_total"], "n_unseen": c["n_unseen"],
            "observed_mass": c["observed"], "unseen_mass": c["unseen_mass"],
            "mean_unseen_pgen": mean_unseen, "mean_observed_pgen": mean_obs,
            "rarity_ratio": (mean_obs / mean_unseen) if (mean_obs and mean_unseen) else None,
            "min_inclusion": c["min_inclusion"], "richness_reliable": c["richness_reliable"],
            "degenerate": c["degenerate"], "reason": c["reason"]}
