"""T-cell precursor frequency for an epitope: how much repertoire mass can see it.

The estimand is

    F(e) = sum over the cognate set C_e of pi(tau)

— the probability that a random naive-repertoire junction recognises epitope ``e``. This is the
continuous quantity behind "immunogenic": recognition is a spectrum, and `F(e)` is where on it a
given epitope sits.

Seven estimators of the same `F(e)`. One is empirical and needs no model at all; the rest ride on
the recombination model, in increasing order of commitment.

``event_ratio``
    `F(e)` counted **directly off repertoire data**, as a ratio of independent recombination events
    to independent recombination events — no `Pgen`, no model. The unit is
    ``(donor, V, J, junction_nt)``: the same nucleotide junction in two donors is *two* events, two
    rearrangements that converged. Because numerator and denominator are counted with the same key
    over the same donors, **repertoire depth divides out** and there is no coverage correction to
    make. This is the estimand itself rather than a proxy for it, which makes it the direct
    empirical check on everything below. Its own weak point is the *other* sampling — the size of
    the cognate set VDJdb holds — which it shares with every set total here; the per-junction median
    form removes it, see the function's docstring.

``observed_mass``
    Sum of `Pgen` over the junctions actually recorded for the epitope. A **strict lower bound**,
    and a biased one: VDJdb samples cognate TCRs **size-biased by Pgen** (a TCR enters the record
    roughly in proportion to its repertoire frequency), so the observed members are systematically
    the high-`Pgen` ones and the deficit does not shrink with more studies at the same depth.

``coverage_corrected_mass``
    ``observed_mass`` with the size-biased deficit put back. A capture curve increasing in `Pgen` is
    fitted to the per-junction donor/study multiplicities, then the mass is Horvitz–Thompson
    reweighted by each member's inclusion probability. Unlike the bound it is *not* monotone in
    depth by construction, and it degenerates loudly rather than silently — see the function's own
    docstring for what it assumes and where it breaks.

``unseen_junctions``
    The richness counterpart of ``coverage_corrected_mass``, from the same fitted capture curve:
    how many cognate junctions were **never** catalogued, and how rare each one is on average.
    Richness and mass are not interchangeable — because the unseen members are systematically the
    low-`Pgen` ones, the missing *count* can be enormous while the missing *mass* stays small.

``union_mass``
    Mass of the **union** of Hamming-`r` balls around the observed junctions. Cognate TCRs are
    near-duplicates by construction, so their balls overlap and the *sum* of per-sequence ball
    masses double-counts; the union is the correct object. Exact, and without enumerating the union.
    The returned ``overlap`` quantifies that double-counting directly, which is worth reporting
    rather than hiding — it is a measurement of how tight the specificity group is.
    (``ball_mass`` is the same quantity by enumeration, kept as its oracle.)

``shell_profile``
    ``ball_mass`` resolved by exact edit distance, so the empirically measured cognacy-retention
    profile `alpha_r` can be applied per shell instead of assuming every ball member is still
    cognate. `F ≈ sum_r alpha_r * mass(shell r)`.

``motif_mass``
    `Pgen` of a degenerate motif (a VDJdb cluster PWM is V/J/length-pinned, hence exactly a
    per-position residue set — see :func:`load_cluster_motifs`). This computes the mass of a **set**
    directly, so unlike the others it does not suffer the observed-sample coverage bias at all.

Which one to use
----------------

============================================  =======================================
question                                      estimator
============================================  =======================================
"what is `F(e)`, measured?"                   ``event_ratio`` (needs donor-resolved
                                              nucleotide repertoires)
"what can I defend without any assumption?"   ``observed_mass`` (report it as a bound)
"how much did the sampling miss?"             ``coverage_corrected_mass`` (needs >= 2
                                              capture units and some recaptures)
"how many TCRs were never catalogued?"        ``unseen_junctions``
"how tight is this specificity group?"        ``union_mass`` -> the ``overlap`` field
"best point estimate from observed TCRs"      ``shell_profile`` -> ``retained``
"an estimate with no sampling bias at all"    ``motif_mass`` (needs a cluster PWM)
"how much mass is still missing?"             ``cross_check`` -> ``missing_fraction``
"how many cells, in this donor?"              ``precursor_frequency`` -> ``cells``
"how many clonotypes, seen and unseen?"       ``occupancy`` -> ``S``, ``n_seen``,
                                              ``n_unseen``
============================================  =======================================

:func:`cross_check` is the scientifically load-bearing one. ``motif_mass`` and the observed-sample
estimators measure the same `F(e)` by independent routes with *different* biases, so **their
disagreement is an estimate of the missing mass**. :func:`event_ratio` then adjudicates from
outside the model entirely.

Two objects, easily confused
----------------------------

``coverage_corrected_mass`` and ``event_ratio`` both take per-junction donor counts, and they mean
**opposite things**:

- in ``coverage_corrected_mass`` the count is a ``multiplicity`` — how many independent units
  *re*-reported the **same** object, i.e. recaptures, from which what was never captured is
  inferred;
- in ``event_ratio`` distinct donors are **distinct objects**, each one its own recombination event.

Feeding the same donor counts to both is silently wrong in one of them. "Seen again" is not
"happened again".

One-substitution matching is part of the estimator
--------------------------------------------------

Everything here defaults to a closed **one**-substitution ball — ``pgen(..., mismatches=1)``,
``ball_mass(r=1)``, ``event_ratio(r=1)`` — and that is a definition, not a tuned parameter. At
radius 0 the repertoire is too sparse for either route to be estimable, and the measured
consequence is large: the Pogorelyy replication scores 0.51–0.61 at radius 0 against 0.76–0.86 at
radius 1. Use radius 0 only to demonstrate that sparsity.

Requires the optional ``vdjtools`` dependency (the recombination model)::

    pip install 'vdjmatch[precursor]'

Nothing here reimplements `Pgen` — the DP, the closed Hamming-1 ball and the degenerate/masked DP
all live in vdjtools, and the neighbourhood enumeration lives in seqtree. Importing ``vdjmatch``
itself never touches vdjtools; the import happens on first use inside this subpackage.
"""
from __future__ import annotations

from .compendium import load_compendium, load_estimates
from .events import RecombinationEvent, event_ratio
from .frequency import (
    ALICE_Q,
    N_T_CELLS,
    SELECTION_BY_CHAIN,
    expected_cells,
    occupancy,
    p_at_least,
    paired_frequency,
    precursor_frequency,
)
from .mass import (
    ALPHA_PER_EDIT,
    MAX_BALL_MEMBERS,
    MAX_COMPONENT_MEMBERS,
    MOTIF_FREQ_THRESHOLD,
    ClusterMotif,
    ball_mass,
    closed_ball_mass,
    cross_check,
    load_cluster_motifs,
    motif_mass,
    observed_mass,
    shell_profile,
    union_mass,
)
from .model import check_junctions, load_model, pgen
from .report import SCHEMA, summarise, summarise_group
from .unseen import coverage_corrected_mass, unseen_junctions

__all__ = [
    "load_model", "check_junctions", "pgen",
    "observed_mass", "coverage_corrected_mass", "unseen_junctions",
    "closed_ball_mass", "ball_mass", "union_mass", "shell_profile",
    "RecombinationEvent", "event_ratio",
    "ClusterMotif", "load_cluster_motifs", "motif_mass",
    "cross_check",
    "load_compendium", "load_estimates",
    "precursor_frequency", "occupancy", "expected_cells", "p_at_least", "paired_frequency",
    "summarise", "summarise_group", "SCHEMA",
    "ALPHA_PER_EDIT", "SELECTION_BY_CHAIN", "MOTIF_FREQ_THRESHOLD", "MAX_BALL_MEMBERS", "MAX_COMPONENT_MEMBERS", "ALICE_Q", "N_T_CELLS",
]
