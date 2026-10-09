Understand matching evidence
============================

A query junction can match several reference observations and several peptide–MHC
labels. The software keeps the sequence search, reference observations, candidate
ranking and background calibration as distinct layers.

Identity and observation counts
-------------------------------

Search uses a unique junction index for each locus and expands its matches back to
reference observations. Identical junctions with different V/J calls remain
separate clonotypes. Independent observations retain their record, study and
pairing metadata; duplicate metadata does not create an additional sequence key.

A candidate is identified by peptide, MHC allele fields and MHC class. The same
peptide with two restrictions produces two candidates. ``annotate`` preserves every
input row and input order, including invalid and unmatched rows. ``candidates``
returns a long table of candidate evidence indexed by query identity.

.. list-table:: Evidence counting units
   :header-rows: 1
   :widths: 30 70

   * - Field
     - Quantity
   * - ``n_records``
     - Matched reference observations for the candidate.
   * - ``n_clonotypes``
     - Distinct matched reference junction/V/J keys for the candidate.
   * - ``n_studies``
     - Distinct nonempty reference identifiers supporting the candidate.
   * - ``n_competing_clonotypes``
     - Distinct matched clonotypes supporting other candidate labels, including
       shared clonotypes when they also support another label.
   * - ``n_reference``
     - Distinct junctions in the candidate's reference set under the stated predicate.
   * - ``n_control``
     - Distinct control junction matches in the search scope.
   * - ``n_hits``
     - Matched reference observations across the query's single-chain candidates;
       paired output identifies its paired-reference support separately.
   * - Query ``count``
     - Supplied abundance in the producer's unit; retain whether it means reads,
       cells, molecules or another count. It is not a reference evidence count.

Candidate ranking
-----------------

NED-v1 ranks candidates by the sum of control-adjusted neighbour weights. Each
reference clonotype contributes an exponentially decaying weight based on its
search score, divided by one plus the number of control neighbours at or within
that hit's edit distance. The existing default score temperature is 400; the
existing soft-V weighting gives same-family matches weight 1 and known cross-family
matches 0.25 times their germline CDR1/CDR2 similarity. Missing or unknown V calls
provide neutral evidence in the production interface.

Without a background, the control count is zero and ranking is explicitly
``uncalibrated``. ``ned_score`` is a ranking statistic. It is distinct from edit
distance, the search penalty ``score``, the database evidence score ``vdjdb_score``,
and background expectation ``E``.

``rank`` starts at 1 within a query. ``competing_score`` is the strongest alternative
candidate's score, and ``score_margin`` is the candidate score minus that alternative.
The highest-scoring candidate is called only when its ranking is unique; an exact
tie produces ``ambiguous``. All alternatives remain available in the candidate table.

Fixed-ball calibration
----------------------

For one query and one candidate reference set, let ``N`` be the number of distinct
reference junctions, ``M`` the number of distinct background junctions and ``n_control``
the background matches under the same search predicate. The expected background
match count is::

   E = N * n_control / M

The observed target count is the number of distinct reference junction matches in
that candidate set. The Poisson-tail ``p_enrichment`` asks how often the background
expectation would produce at least that observed count. ``p_any`` is the probability
of at least one match under that expectation.

A finite background with zero matches does not establish zero background probability.
The existing rule-of-three path uses ``E = 3 * N / M`` and reports
``rule_of_three=True``. NED's ``+1`` denominator is a separate ranking convention;
it does not replace this finite-control calculation.

Target counts, reference size and background size must have the same counting unit.
The target and background search predicates must also agree. The current sequence-only
background cannot represent V/J restrictions, so hard-V/J calibrated requests fail
explicitly. Requested calibration failures propagate rather than silently switching
to an uncalibrated call.

Neither NED, ``E`` nor ``p_enrichment`` is a posterior probability that the winning
peptide is the receptor's specificity. Competing labels, cross-reactivity and incomplete
reference coverage remain visible in the candidate evidence.

Paired evidence
---------------

A paired candidate requires alpha and beta matches in the same reference complex
and under the same peptide–MHC label. Alpha/beta agreement without a shared complex
is chain evidence fusion, a different claim. Repeated observations of a paired
receptor remain available as support metadata; ranking uses distinct paired receptor
keys rather than treating repeated metadata as new sequence evidence.

Paired ranking uses the sum of
``exp(-(score_alpha + score_beta) / score_scale)`` over
distinct paired junction keys. ``n_hits`` counts those keys, ``n_records`` counts
distinct supporting reference complexes, and ``n_reference`` counts distinct
junction pairs for the candidate peptide–MHC set. Repeated observations add complex
support without increasing the paired similarity score. The estimator is named
``paired-reference-v1`` separately from single-chain NED-v1. The API default
``score_scale`` is 400 and must be finite and positive; unit-cost CLI scoring uses
scale 1. Detailed paired hits retain the source observations separately, with
``alpha_``/``beta_`` metadata and optional alignment fields. Request them together
with candidates using ``return_hits=True`` without repeating reference searches.

Paired calibration requires both chain backgrounds. It reports
``n_control_alpha``, ``n_control_beta``, ``E``, ``p_enrichment`` and
``rule_of_three`` with ``calibration='paired_independent_fixed_ball'``.
The product of the chain background probabilities uses a stated conditional
independence assumption. It describes the chance of the joint search event under
those backgrounds; it is not a posterior probability of specificity and does not
establish independence in a biological repertoire. Dependence and sample-level
co-occurrence models are outside this calibration. Paired CLI matching rejects
``--match-v`` and ``--match-j`` until predicates can be applied separately to each
chain; these flags are not silently ignored.

The separate historical ``evalue.paired.paired_scan`` evaluates a nearest paired
radius, defined as the maximum alpha/beta cost. Default two-tuple hits count
observation associations. Its optional ``include_identity=True`` emits three-tuples
carrying the two junctions; ``pvalue`` counts distinct pairs in that mode and its
reference size ``N`` must count distinct pairs as well. This adaptive-radius helper
is a separate contract from the Annotator's fixed-ball calibration.

Statuses and sample summaries
-----------------------------

``matched`` means the ranking has a unique top candidate. ``ambiguous`` means an
exact tie. ``no_hit`` means a retained query has no reference neighbour in scope;
``no_joint_hit`` means no same-complex paired candidate. ``no_reference`` means the
required locus is absent from the loaded reference. ``invalid_query`` means the
query cannot supply a valid matching sequence/locus.

CLI summary tables report distinct junction/V/J/locus clonotypes (``unique``),
matched input rows (``n_query_rows``), and their summed abundances (``reads``).
Paired summaries count matched query pairs (``n_query_pairs``).
They do not implement repertoire-level statistical enrichment, donor comparisons
or a Pgen-matched naive repertoire analysis. Those require a separately defined
sample-level estimand and background. The optional precursor APIs remain separate
from query annotation and background-match calibration.

For paired finite controls, each zero-neighbour chain uses the existing rule-of-three
upper bound ``min(1, 3/M)``. Capping at one matters for tiny controls and prevents
pseudo-counts larger than the control size. Both paired implementations share this bound.

Global fixed-K statistics
-------------------------

For a fixed query, score rule and predeclared K, let T be the Kth-smallest
distinct-reference score and N the original reference exposure. Under an IID
background score CDF F, the null CDF of this order statistic is
``BinomialSF(K-1, N, F(T))``. This remains conservative for discrete score ties.
An independent control of size M supplies C scores at or below T. The one-sided
Clopper–Pearson upper bound U uses failure budget δ=10⁻⁶; the reported test is
``min(1, δ + BinomialSF(K-1, N, U))``. Conditioning on target-selected T leaves
the independent control unchanged, so a pointwise bound suffices. This is the
confidence-set/failure-budget construction applied to a global order statistic,
not the fixed-ball plug-in Poisson statistic.

For paired score ``max(score_A, score_B)``, the same construction uses marginal
control bounds for the Cartesian intersection. Exact/exact exclusion reuses the
four-category bound described above; without exclusion two marginal CP bounds
share δ equally. An additive paired score threshold cannot use this rectangle.
Global enrichment is separate from per-pMHC ranking; it is not the probability
that the winning epitope assignment is correct.
