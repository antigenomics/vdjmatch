Two directions of control calibration
=====================================

Select one raw control repertoire. Calibration is computed on demand after
reference filtering; no fitted calibration file is needed. The
``database-evidence --control`` command implements the two directions below
for full-junction scoped nearest distance or distinct-neighbour count.
For complete-score single, mixed and genuine paired inputs, use the opt-in
``annotate --control`` route below. The scoped ``database-evidence`` statistics
remain available separately.

Complete-score annotation
-------------------------

.. code-block:: bash

   vdjmatch annotate sample.tsv --vdjdb reference.zip --control control.tsv --threads 4 -o calibrated

With default ``--background vdjdb-other``, the control opts into the combined
reference-weighted score and both on-demand tails. ``--mock-replicates`` defaults
to 199; ``--seed`` defaults to zero. No fitted calibration file is required.
Reference epitope/MHC, ``--evidence-type`` and ``--min-score`` filters apply first.
Ordinary annotation includes every valid selected reference observation; explicit
``--reference-population curated-benchmark`` retains benchmark curation.

The control must supply full junctions and unambiguous V/J annotations. Counts
are validated before filtering; false/missing productive tags and invalid or
missing amino-acid junctions are excluded with separate counts. Without a
productive column, valid amino-acid junctions remain eligible. Query rows are
retained even when invalid or from another species. Paired controls require TRA
and TRB, randomized independently. Incomplete genuine pairs retain the missing
partner row and calibrate each eligible chain as a marginal query.

Outputs are four flat TSVs: ``.queries.tsv`` preserves every raw chain row and
calibration identity/status; ``.database.tsv`` tests the maximum across labels
and available axes; ``.epitopes.tsv`` reports fixed full-pMHC tails and score
components; ``.calibration.tsv`` records filters, populations, model/native hashes
and seeds. Original query identifiers remain intact; generated calibration
identifiers are unique across input files and populations. No JSON input or
output is required by this route.

The combined statistic uses an independent unconditional normalization bank
and recomputes reference-derived weights on every mock model. Final query-null
draws remain matched to the observed Pgen strata in both chains.
Paired score axes include linked evidence and available unlinked marginals; the
axis maximum is recomputed within every draw. The original component prior uses
the matching-chain control for singles and the TRB control for paired scoring.
This candidate one-control context differs from accepted historical benchmark
model assets. It remains opt-in pending whole-task nonregression checks; the
no-control stem4 default and older real/generative control routes are unchanged.
Nondefault geometry and target-from-sample selection are rejected before mocks.

The event being tested
----------------------

Let :math:`q` be a query clonotype and :math:`R_i` the reference for one complete
peptide--MHC target. Let :math:`R_0` be the whole selected reference. Species
and loci define separate populations. Write :math:`T_A(q,R_i)` for the final
higher-is-better statistic under procedure :math:`A`. This procedure includes
scope, distance geometry, V/J restrictions, exact-match policy, score formula,
weights, reference filtering and any target or paired-axis selection.
For nearest distance, :math:`T_A=-d_{\min}`; for count, it is the number of
distinct matching junctions within the declared scope. A no-hit is weaker than
every finite supported score. Exact gene-clonotype presence remains a separate
deterministic lookup.

Null reference populations preserve each group's source-record count and
rounded log2-Pgen histogram. The current fixed bin width is one bit, with zero
model mass in a separate bin. These are approximate same-Pgen strata, rather
than an assertion of equal exact Pgen. Source records and distinct search keys
remain separate quantities: repeated observations contribute to the mock
population's requested size, whereas the scoped count statistic deduplicates
matching junctions.

The control is collapsed at amino-acid junction level using rearrangement-event
weights. Read counts do not define these weights. Pooled ``occurrences`` counts
are used when supplied; otherwise distinct nucleotide/V/J events within samples
are counted. Pooled nucleotide keys can carry their supplied sample incidence.
When only amino-acid convergence is supplied, that available event unit is used;
sample incidence times convergence is never assumed to recover an unknown
number of events. The control's supplied population metadata selects the
declared species/locus. No anchors are fabricated.

Random reference, fixed query
-----------------------------

Draw :math:`R_i^{*(b)}`, :math:`b=1,\ldots,B_R`, from the control, matching
the Pgen strata and source-record count of :math:`R_i`. Within a group,
sampling draws integer rearrangement-event slots without replacement, then
collapses to amino-acid search keys. Distinct nucleotide events producing the
same amino-acid junction can therefore occur together in a mock reference.
Each replicate is generated and reduced in memory. An insufficient stratum
raises an error instead of silently changing the requested null.

For observed score :math:`t=T_A(q,R_i)`, define

.. math::

   p_R(q,i)=\frac{1+\sum_{b=1}^{B_R}
       \mathbf1\{T_A(q,R_i^{*(b)})\ge t\}}{B_R+1}.

This asks whether the observed reference supports the query more strongly than
a composition-matched random reference. Chance exact overlaps remain eligible
in the mocks. Removing every real reference junction from the control would
make a chance exact mock match impossible for an exact database query.

Random query, fixed reference
-----------------------------

Draw :math:`q^{*(b)}`, :math:`b=1,\ldots,B_Q`, from the query's Pgen stratum
in the same control, with replacement and event weights. Keep :math:`R_i`
fixed and evaluate the same complete procedure:

.. math::

   p_Q(q,i)=\frac{1+\sum_{b=1}^{B_Q}
       \mathbf1\{T_A(q^{*(b)},R_i)\ge t\}}{B_Q+1}.

This asks whether the query's support is unusual among clonotypes with similar
generation probability. It uses the actual reference motif structure and its
overlapping neighbourhoods.

A related quantity answers the literal match-probability question:

.. math::

   \widehat\pi_Q(q,i)=\frac1{B_Q}\sum_b
      \mathbf1\{q^{*(b)}\text{ has any match to }R_i\text{ under }A\}.

``null_query_match_rate`` reports this Monte Carlo estimate, which can be zero.
It differs from ``p_query``: a score-tail test compares against the observed
score, whereas the match rate counts every supported match within scope.
``null_reference_match_rate`` analogously estimates any-match frequency for
the fixed query against mock references. These rates do not estimate the
probability that an observed specificity assignment is false.

Combining the directions
------------------------

Report both directional P-values and

.. math::

   p_{\mathrm{both}}=\max(p_R,p_Q).

Let :math:`H_R` and :math:`H_Q` denote the respective conditional randomization
nulls. The combined null is :math:`H_R\cup H_Q`: at least one direction is
explainable by its control law. If either null, say :math:`H_R`, holds and its
P-value is valid, then, for any :math:`a\in[0,1]`,

.. math::

   \Pr\{p_{\mathrm{both}}\le a\}
   \le \Pr\{p_R\le a\}\le a.

Thus the maximum tests evidence in both directions without independence.
This is the intersection--union construction. It cannot rescue an invalid
directional null. Multiplication or an independent Fisher reference
distribution would count overlapping information without a dependence model.
For example, :math:`p_R=0.002,p_Q=0.3` gives :math:`p_{\mathrm{both}}=0.3`;
the reference looks unusual, but many control clonotypes have comparable
support. Two small directional tails support unusual resemblance in both
directions; neither tail is a posterior binding probability.

The plus-one construction includes the observed statistic in its rank and
retains ties and no-hit draws. Its finite-sample validity requires exchangeable
observed/null statistics under the specified conditional model; a fitted
Pgen approximation alone does not establish that assumption. See
`Phipson and Smyth <https://doi.org/10.2202/1544-6115.1585>`_ for Monte Carlo
rank tails and the
`intersection--union discussion <https://pmc.ncbi.nlm.nih.gov/articles/PMC4783210/>`_
for conjunction inference. The maximum proof above applies directly here.

Analytic checks for scoped matches
----------------------------------

For a fixed query, suppose stratum :math:`h` contains :math:`W_h` control
events, of which :math:`c_h(q)` match, and a mock reference draws
:math:`n_{ih}` event slots without replacement. Independent stratum sampling gives

.. math::

   \Pr\{\text{any match}\}
   =1-\prod_h\frac{\binom{W_h-c_h(q)}{n_{ih}}}
                         {\binom{W_h}{n_{ih}}},\qquad
   \mathbb E[K_i^{\mathrm{events}}(q)]
   =\sum_h n_{ih}\frac{c_h(q)}{W_h}.

The empty-match factor is zero when too few nonmatching keys remain.
The expectation counts matched events, not distinct amino-acid keys or a P-value.
If matching junction :math:`r` has :math:`w_r` events in stratum :math:`h(r)`,
the expected number of distinct matching junctions is instead

.. math::

   \mathbb E[K_i^{\mathrm{keys}}(q)]
   =\sum_{r:\,r\text{ matches }q}
     \left[1-\frac{\binom{W_{h(r)}-w_r}{n_{i,h(r)}}}
                        {\binom{W_{h(r)}}{n_{i,h(r)}}}\right].

Under independent
with-replacement event draws with match probability :math:`u_h(q)`, the
any-match probability instead is
:math:`1-\prod_h(1-u_h(q))^{n_{ih}}`. The Poisson approximation
:math:`1-\exp[-\sum_hn_{ih}u_h(q)]` requires sufficiently rare individual
matches. Integer event-slot sampling has this hypergeometric law even when
different amino-acid keys represent different event counts. Sampling weighted
unique amino-acid keys without replacement would have a different law and
would suppress convergent rearrangements.

For random queries, the exact model probability is the mass of the union of
all matching reference balls within the query's Pgen stratum. A native batch
against a generated repertoire estimates this union directly and counts a
query once, even when several balls overlap. Summing per-reference ball masses
double counts intersections. Formula and sampling checks agree only when
they share the event unit, conditional generative/selected law, scope and
reference. A raw generative model need not reproduce a selected empirical
repertoire exactly.

Complete-score and selection requirements
-----------------------------------------

For reference-derived weights :math:`w(R)`, the null statistic is
:math:`T_A(q,R^*)` with :math:`w(R^*)` recomputed. Freezing weights chosen from
the observed reference omits part of the procedure. Component normalization
must also treat observed and null scores symmetrically, or use a separate
internal fitting split from the same selected control. One user control does
not require reuse of the same draws for fitting and testing.

For a fixed pMHC, these are fixed-target tests. A winning label, maximum over
targets, paired-axis choice, radius selection or combined-score selection
must be repeated inside each null replicate. A whole-reference scoped
nearest/count statistic already includes all reference junctions, but is not
the null distribution of a maximum VDJmatch score. Independently sampled
per-target mocks also do not define a coherent mock for selected labels;
that event must preserve overlapping memberships and genuine reference pairing.

The paired null draws TRA and TRB independently, separately matching the two
query Pgen strata: :math:`Q_{\alpha\beta}=Q_\alpha Q_\beta`.
``iter_mock_paired_references`` implements this law with separate random streams
and independent within-target pairing permutations, preserving source-pair
counts and both marginal Pgen histograms. It accepts one raw control table
containing both locus populations. Independently sorted lanes must not simply
be zipped: sorting would introduce an artificial association.

Even under this product null, the statistic must match both chains to the same
linked reference pair and retain all unlinked same-target marginal support.
The random-query any-call probability is

.. math::

   \pi_{Q,\alpha\beta}(q,i)=
   \sum_{x,y}Q_\alpha(x\mid h(q_\alpha))
              Q_\beta(y\mid h(q_\beta))
              \mathbf1\{A(x,y,R_i)\text{ has any linked or marginal call}\},

where :math:`h` denotes the chain-specific Pgen stratum. Multiplying two
chainwise *any-database-match* rates would permit the two matches to come from
different linked reference pairs and would omit the declared marginal/axis
selection. A sum-distance ball against one fixed reference pair instead uses
the convolution :math:`\int F_\beta(r-u)\,dF_\alpha(u)`, where :math:`r` is
the total radius and :math:`F_\alpha,F_\beta` are independent chain-distance
distributions to that pair. The union across reference pairs still needs its
overlap handled by the complete native batch. The ``annotate --control``
adapter repeats these full-score and selection operations. Whole-cohort runtime
and scientific nonregression remain separate acceptance gates.

Changing evidence filters or supplying a custom labelled table changes
:math:`R`, its size and its target compositions. Recalibrate both directions
on that selected table with the same procedure. ``p_both`` then describes
unusual matching under the declared controls. It does not prove binding,
specificity exclusivity or a biological false-discovery rate.
