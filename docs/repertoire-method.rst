Repertoire enrichment: algorithm and parameter choices
======================================================

Implementation status
---------------------

The four established comparison families are unit, V-unit, weighted and
V-weighted. They cross unit/weighted edit-bounded search with absence/presence
of hard V equality. Their NED score also retains the existing soft V weight
when hard equality is disabled. The implementation is owned by
``Annotator.graded_candidates`` and the shared candidate reducer. The M3
count adapter reuses the existing repertoire-indexed native batch and shared
gene normalization, omitting NED rankings and per-mock Annotator construction.

M3 now exposes those four geometries through the same matcher names. Supply
``--tra-radius-factor`` or ``--trb-radius-factor`` as an explicit upper bound
on an individual alignment's dimensionless distance: raw cost divided by 1
for unit families and 400 for weighted families, using the existing graded
score scales. ``--mismatches`` is the total edit ceiling (1–3), with insertion
and deletion ceilings min(radius, 2). Weighted families use the existing VDJAM
matrix and gap costs; they are not the normalized BLOSUM presets.

A receptor contributes one hit if at least one alignment passes the bound,
including equality. Its abundance contributes once. This is a nearest-distance
endpoint: soft-V NED weights and summed weak evidence do not enter acceptance.
The ``v-`` prefix enforces gene equality; without it V does not filter hits.
Optional J equality remains separate. Real and mock panels use the same rule,
with exact hits retained and each panel's complete model preserved.

The bulk ``repertoire-exact`` command accepts the same four names through
``--matcher`` and the bound through ``--radius-factor``. It preserves its original
panel sampler and source-row/read units, including repeated source rows. Exact J
is retained in every bulk family; V equality applies only to ``v-`` families.
For these graded families, substitution costs are symmetric and insertion and
deletion budgets/costs are equal, so exchanging the two alignment rows preserves
the minimum cost. The bulk count path therefore indexes the repertoire once and
batches the smaller real/mock panel junctions as queries. Resolved V/J predicates
are joined to the sparse hits before deduplicating source observations. This
reuses the existing native search for large-repertoire/small-panel workloads;
receptor-level NED rankings keep their own query-conditioned interpretation.

Software parity is checked against the existing graded API. The fourteen-setting
bulk comparison and the separate pooled/sample-specific conditioning assessment
are complete; their observed operating choices are reported below. The paired
four-family comparison is still running. Software parity alone does not select
a biological cutoff. Continuous-score summation remains deferred. Normalized
BLOSUM results retain their separate presets and producing settings.

Start with the biological question
----------------------------------

Receptor annotation asks which reference targets support a sequence. Repertoire
enrichment asks whether a fixed sample supports a target more strongly than
matched mock references. Exposure summaries ask whether that evidence is
consistent with prior exposure. These are three different outputs.

For repertoire work, retain junctions including conserved anchors, V/J calls,
sample identities, receptor pairing and abundance. Keep peptide and MHC together.
A pair must refer to the same receptor; separate alpha and beta matches anywhere
in the sample are not a paired match. Count each receptor once per target and
once in an across-target union. Read, cell and receptor counts are different
endpoints. Technical replicas are not independent donors.

Recommended choices from completed M3 comparisons
-------------------------------------------------

.. list-table:: Current recorded operating choices
   :header-rows: 1
   :widths: 20 40 40

   * - Property
     - Bulk longitudinal cohort
     - Paired single-cell cohort
   * - Shipping command
     - repertoire-exact
     - repertoire-enrichment --null reference
   * - Current cohort choice
     - positional-blosum62, radius factor 0.30
     - Hamming <=1 substitution on each chain
   * - V / J predicates
     - Gene equality for both V and J
     - Gene equality for V; J equality off
   * - Count units
     - Source rows and reads, once per panel
     - Unique sample-receptors and cell abundance, once per target
   * - Pgen / composition
     - Exact joint olga-legacy, Q=1; pooled J|V correction
     - Exact marginal isalgo; fixed pgen-matched reference null
   * - Mock settings
     - 200 replicates, seed 42
     - 200 replicates, seed 42
   * - Reported decision
     - Original conditional tail, donor/day BH q<0.10 and positive effect
     - Chain nominal P<0.05 with positive excess; cohort BH retained
   * - Observed coverage
     - 5/7 day15 libraries by rows, 7/7 by reads; 0/12 day0
     - 5/6 day14 donors, 6/6 ever post; 0/6 baseline, paired descriptive heuristic

These settings reproduce the currently strongest recorded cohort endpoints.
They are not promoted software defaults; the paired four-family comparison
must finish before its final recommendation is selected.

There is no single matching setting demonstrated optimal across both cohorts.
For the bulk longitudinal cohort, positional BLOSUM62 radius factor 0.30 improved
significant day15 clonotype libraries from 4/7 to 5/7 while retaining 7/7 read
calls and 0/12 day0 calls on both endpoints. This is the selected observed
tradeoff among tested settings, not an independently validated universal optimum.

For paired single-cell exposure coverage, retain the demonstrated Hamming<=1
reference: 0/6 prechallenge donors, 5/6 day14 donors and 6/6 donors at some
postchallenge timepoint. All six radius-0.30 substitution arms retained the
14 supplied positive assignment rows, but none improved that coverage. Those
14 rows represent eight sample-receptors and three target tests in one sample,
and their sequences are present in the source-publication reference.

Use exact Pgen (``--pgen-mismatches 0``) as recorded for the accepted flexible
reference-null comparisons. Matching radius and Pgen neighbourhood radius are
separate parameters. Increasing mock replication improves tail resolution; it
does not repair a matcher that also admits much more null support. At 200 mocks
the minimum empirical P-value is 1/201. Retain unavailable-null rows as missing
inference, not negative evidence. Preserve nominal and BH results and name the
original correction family.

Large repertoire, small target panel
------------------------------------

With N=10^7 repertoire observations and M=100 reference keys, enumerating every
comparison means N*M=10^9 sequence pairs. The counting task only needs qualifying
memberships. Exact membership uses the existing columnar key join; bounded graded
matching uses the existing native repertoire index and batches the smaller panel
as queries. For B=200 mock panels of M keys, deduplicating the real and mock
query keys gives at most (B+1)*M=20,100 queries for this one target, often fewer.
Search direction can be exchanged for symmetric graded costs with equal
insertion/deletion budgets, while preserving V/J predicates and source-row/read
weights. It cannot be exchanged on that argument alone for position-dependent
or query-conditioned scores. The current positional implementation groups
panel queries by length and uses its producing geometry.

This reduces redundant traversal and avoids a dense pair matrix; it does not
prove constant search time or biological optimality. Index construction, emitted
hit count, joins, mock preparation and all targets still contribute to runtime.
The current measured cohorts do not establish ten-million-observation scaling.
Nor does a small panel make weak per-receptor hits harmless: if b is the null
probability that one observation has any qualifying target match, linearity
gives E[K]=N*b, even without independence. At b=10^-6 and N=10^7 that is ten
background observations. Matching geometry therefore remains a signal-versus-null
choice evaluated through the declared count endpoint; changing the search
direction accelerates that same test rather than selecting its threshold.

Two algorithms, two null hypotheses
-----------------------------------

``repertoire-exact`` reproduces the historical bulk-panel calculation, including
source-row/read units, joint probability strata and its recorded sampler.
``repertoire-enrichment --null reference`` fixes the sample and generates mock
target references using marginal full-junction Pgen. Do not exchange their
P-values, denominators or significance conventions.

For the reference-null route:

1. Load and normalize the raw sample, curated reference and raw controls.
2. Compute the selected model's full-junction Pgen and rounded log2 strata.
3. Retain the full reference model and target inventory. Prepare control-bin index arrays once, then draw the required number
   of control event slots with replacement for each reference bin. Reference
   values below the smallest control bin use that bin; repeated clonotypes
   remain eligible and observed receptors are counted once. Optional V/J-conditional policies change this sampling law.
4. Apply exactly the same locus-specific matching and reducer to real and mock
   panels, in native batches. Keep zero-hit observations.
5. Export observed support, mean mock support, excess, empirical upper tails,
   BH results, identities and complete hits. Keep chain inference separate from
   downstream descriptive exposure summaries.

For observed support K and mock supports K_b, with B mock panels, the exported
upper tail is (1 + number of K_b >= K)/(B+1). Abundance uses the same formula
with each receptor contributing its original count once. This arithmetic is
exact for saved mock draws. Its interpretation as a calibrated test requires
that the conditional mock law adequately represents the biological null;
matching Pgen alone is not a proof of exchangeability.

Does the construction force Hamming 0 or 1 to win?
--------------------------------------------------

No. A nested matcher increases both observed and null counts; there is no
monotonicity theorem for their tail probability. Here is a constructive
counterexample. Under exact matching let observed support be 1 and every mock
support be 2, giving P=1. Under a broader rule add three real neighbours and no
mock neighbours: observed support is 4 and every mock remains 2, giving
P=1/(B+1). Both count arrays obey nesting. Conversely, start with observed 2
and mocks 0, then add three mock neighbours and no real neighbours: broadening
changes P from 1/(B+1) to 1. Distinct receptor identities can realize these
counts around different real/mock references. Thus neither exact nor one-edit
matching is optimal by construction; no theorem singles out radius one either.

A useful conditional criterion is the signal-to-null ratio of the added shell.
For a binary membership feature X, let a=P(X=1|specific) and b=P(X=1|null).
The hit likelihood ratio is a/b when b>0. Enlarging the set by disjoint masses
u and v improves this ratio exactly when ub>av, by cross multiplication of
(a+u)/(b+v)>a/b. This criterion concerns per-receptor discrimination, not a
proof of enrichment-test power. The latter also depends on sample size,
abundance, reference diversity, null variance and tail discreteness.

Potential design bias still deserves attention: an exact-Pgen stratum does not
necessarily match neighbourhood mass for broad weighted distances. The earlier distinct-key comparisons excluded known target keys; the
author-requested binned control draw now retains incidental matches and clips
values below the control minimum. Earlier results retain their producing law. These observations motivate a comparison; they
do not establish that the existing result is invalid or guarantee an alternative
will improve it. No new null model or threshold sweep is introduced here.

What is optimal under a stated model?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fix the reference, matching rule and threshold before testing a sample of N
independent receptors. Let X_i indicate whether receptor i has at least one
qualifying match, and K=sum(X_i). Suppose every X_i is Bernoulli(b) under the
background hypothesis and Bernoulli(a) under the response hypothesis, with
0<b<a<1. The likelihood ratio is

L(K)=(a/b)^K ((1-a)/(1-b))^(N-K).

Increasing K by one multiplies L by a(1-b)/(b(1-a)), which exceeds one because
a>b. For a fixed background rejection probability, assigning rejection to
outcomes with larger likelihood ratios maximizes response probability: swapping
equal background mass from a smaller ratio to a larger one increases the
response mass. Therefore an upper-tail test on K is the most powerful test at a fixed
false-positive rate for these two simple Bernoulli hypotheses (randomizing at
the boundary if needed). This proves optimality of the count ordering under
that model, not optimality of a particular matching geometry, Hamming radius
or V predicate. Those choices change a and b and must be compared empirically.
Read-weighted support has a different sampling model; this proof does not
make reads independent receptor observations. M3 retains its matched mock tails
rather than replacing them with this illustrative binomial model.

Training determines a scoring rule; calibration determines an operating point
for a specified event and reference/null population; evaluation measures its
performance on the declared test population. The earlier selected-kernel
calibration is therefore useful evidence for its receptor-level event. It does
not select the most powerful sample-level count test across these four
alignment families. M3 can keep the trained scoring parameters fixed while
comparing individual-distance cutoffs using the recorded sample endpoints.

The completed evidence addresses three different endpoints:

* Earlier sequence-level comparisons associate shared V usage with greater
  neighbour co-specificity and favour local, position-aware evidence. This
  motivates retaining V information without proving hard equality optimal.
* The selected-kernel generated-null check fits six curated species/chain
  contexts at a nominal 1% withheld-call target, then evaluates 100,000 generated
  queries per context. It retains no-score queries and paired axis selection.
  Those cutoffs belong to that score, reference and selected-call event.
* M3 compares thresholded membership followed by sample-level counts under its
  fixed mock-reference law. A broader matcher can admit both more response
  receptors and more mock support. Its detection endpoint therefore requires
  the matched cohort comparison, even with scoring parameters held fixed.

Thus the previous training and calibration remain useful: they establish their
own scoring and operating-point evidence. They do not analytically settle M3's
choice of distance threshold or V predicate.

A low weighted bound can equal exact matching by construction. In the bundled
VDJAM matrix, the smallest nonidentity substitution penalty is 76; gap opening
and extension cost 100. All costs are nonnegative. At weighted bound 0.1 the
budget is ``400 * 0.1 = 40``. Every nonidentity alignment exceeds that budget,
whereas an identical junction has cost zero. Consequently weighted bound 0.1
and unit bound zero have identical membership when their V/J predicates and
inputs agree. This equivalence concerns thresholded counts, not NED scores or
the biological optimality of exact matching. Higher bounds can admit genuine
substitutions or indels and need their own comparison.

V and J matching are part of the hypothesis
-------------------------------------------

Exact V matching keeps germline context aligned and narrows accidental support;
it can lose real responses across related genes or uncertain calls. Allele-pooled
gene equality and allele equality must be distinguished. No-V matching can
recover missing/cross-gene support but expands the null neighbourhood. A fuzzy V
rule based on fixed germline-loop geometry is biologically plausible, but must
be defined in comparable distance units rather than an arbitrary gene-name edit.
J matching should be stated separately from V and from junction distance.

The strict Hamming and normalized BLOSUM routes retain the documented V
predicate, with optional J equality. Kernel routes use their own production
V-loop geometry and score semantics; this is not a controlled V-only comparison.
Do not attribute differences between those routes exclusively to V matching.
The completed bulk exact-junction count comparison holds reference, sample
cohort, Pgen model, sampler, J predicate and aggregation fixed. Hard V equality
retains 4/7 day15 clonotype and 7/7 read calls, compared with 2/7 and 5/7 without
the hard V filter; both have 0/12 day0 calls. The strict-V route reproduces the
accepted exact matched counts and all call flags across 42 libraries. These
results favour enforcing V for this exact-junction bulk endpoint; they do not
assume that cross-V recognition is impossible. The separate day-1 read call
remains in the complete record.

The fourteen-setting bulk comparison is complete under its fixed pooled null.
Hard-V weighted bounds 0.1 and 0.2 preserve every exact-baseline call; none of
these settings improves on the recorded positional BLOSUM62 bound-0.30 result.
An isolated fuzzy-V contrast remains a proposal: keep the junction rule and
other settings fixed, and apply the V rule identically to real and mock panels.
Preserve raw missing gene status and imputation; imputed equality is not an
observed gene call.

V matching and V-usage correction are separate decisions. The historical bulk
correction is ``P_sample(J|V) / P_model(J|V)``, with add-one frequencies on the
declared V/J support. It adjusts J usage within V. A joint-frequency correction
would multiply this ratio by ``P_sample(V) / P_model(V)``. This identity follows
from ``P(V,J) = P(J|V) P(V)``; it does not establish which null best represents
the biological background. Pooling samples versus fitting each sample also
changes those frequencies and hence the Pgen strata and mock law. Compare those
choices with matching held fixed, preserving observed counts and common final
FDR families. The completed 42-library comparison retains 4/7 day15 clonotype
and 7/7 read calls for hard-V exact under either correction. Positional bound
0.30 gives 5/7 clonotype calls pooled versus 4/7 sample-specific, retaining 7/7
read calls and 0/12 day0 calls on both endpoints in either case. Pooled
correction therefore retains the strongest observed bulk result on this
endpoint. This comparison does not introduce or select a marginal-V correction.

When can hard V equality be optimal?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let H1 mean target-specific recognition and H0 mean matched background. For
junction comparison d and V relationship g, the full likelihood ratio is
L(d,g)=P(d,g|H1)/P(d,g|H0), wherever the denominator is positive. A sufficient condition for a hard V
filter to be lossless is the strong assumption
P(g is unequal|H1)=0 (and correctly observed V calls). Then every unequal-V
case has zero numerator and can be rejected by an optimal likelihood-ratio
rule. This assumption is not established by the present results. If unequal V
occurs under H1, hard equality can discard real signal. At a fixed likelihood-ratio cutoff c, the weaker condition that every
excluded unequal-V case has L(d,g)<=c is already sufficient (with ties handled
consistently). Conversely, dropping V is lossless across operating points if
L(d,g) is independent of g conditional on the retained features. Neither
sufficient condition follows from Pgen matching.

Thus exact V, fuzzy loop geometry and no V correspond to different assumptions
about germline contribution. Test their matched-case effects before declaring
one optimal; do not justify hard equality simply because a database contains
many same-gene neighbours. Gene imputation can further violate the correctly
observed-V assumption and should be audited separately.

Combining chains and targets
----------------------------

Use separate TRA/TRB tests as the inferential unit. Same-receptor intersection
supports a stronger paired biological statement than union; union preserves
single-chain coverage and orphan receptors. Neither is a calibrated exposure
P-value in the current manuscript interpretation. Sum counts only after
receptor deduplication: one receptor matching several epitopes must not become
several observed cells in an exposure denominator.

The actual exposure exporter first restricts each sample to compatible COVID
pMHC tests with available chain inference. For peptide i and compatible HLA
entry h, let p_ih be its empirical chain tail and e_ih its observed excess.
It sets u_ih=p_ih when e_ih>0, and u_ih=1 otherwise, then uses
p_i=min_h(u_ih) over available HLA entries for that peptide. A peptide with
only nonpositive excess remains present with p_i=1 and contributes zero to
the sum. The exported k counts all retained peptides, including those terms;
unavailable tests remain explicit. This minimum is a selection within each
peptide, so independence between peptide labels alone cannot give p_i a
uniform null distribution. No HLA-minimum correction is applied or implied.

For evidence from k distinct epitopes, Fisher's score is
T=-2 sum(log(p_i)), where p_i are the specified chain-level P-values. Under
independent continuous uniform null P-values, -2 log(p_i) has a chi-square
law with 2 degrees of freedom, and the sum has 2k degrees of freedom. This
follows from P(-2 log U >= t)=exp(-t/2) for uniform U. The continuous-uniform
assumption is not exact for discrete 200-mock tails; sharing receptors and
reference motifs may also induce dependence. Distinct epitopes are treated as
independent in the manuscript framework, but that is an assumption, not a
consequence of different peptide labels.

Keep T as an explicitly labelled evidence score under the present contract.
Do not call it a calibrated exposure P-value. Do not combine TRA and TRB as
independent replicates or combine a chain P-value with a union/intersection
P-value built from that chain. Fix the eligible target set before combination,
retain non-significant tests, and report k; selecting only significant P-values
would change the null distribution. Missing tests must remain explicit.
The statistical assumptions are documented in the `SciPy combine_pvalues
reference <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.combine_pvalues.html>`_.

Refinement priorities
---------------------

* First compare V rules while holding junction matching fixed; this separates
  germline restriction from the substitution metric. Success means better
  matched-case coverage without extra prechallenge calls, with every loss kept.
* Keep separate chain evidence, same-receptor paired support and cross-target
  exposure unions visible. This improves interpretation without fitting a model.
* Report Fisher score, component tests and k together. Any new calibrated
  combination requires its own explicit null specification and authorization.
* Check whether matched neighbourhood mass explains broad-rule dilution before
  proposing a changed sampler. Preserve exact-Pgen accepted results throughout.

These are proposals with distinct verification requirements, not measured gains.
No new experiment, fitted rule, cache, null model or software default is added.

Runtime records for later optimization
--------------------------------------

Bulk CLI manifests include ``stage_seconds`` and ``sample_timings``. Stage
timings separate input/identity checks, reference/probe loading, Pgen/selection,
the initial source V/J pass, strata, mock generation, sample loading/matching,
count assembly, statistics/FDR and export. Per-library records retain source
row counts and loading versus matching/staging time; divide by those counts
when comparing differently sized libraries. For stem4, sample staging and the
subsequent cohort kernel stage are reported separately.

The per-library timings break down the aggregate stage and must not be added
to it again. Manifest timings cover the CLI body; the external whole-process
profile includes startup and records sampled CPU, threads and summed RSS.
Use both records to identify expensive stages before changing algorithms.
Frozen existing runs retain their original whole-run resource records; missing
internal stage times cannot be reconstructed exactly from untimestamped logs.
