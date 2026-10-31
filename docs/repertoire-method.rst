Repertoire enrichment: algorithm and parameter choices
======================================================

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
3. Retain the full reference model and target inventory. Sample unique control
   junction/V/J keys without replacement within target/bin, excluding known
   target keys. Optional V/J-conditional policies change this sampling law.
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
necessarily match neighbourhood mass for broad weighted distances. Known-target
exclusion also means the sampler is a conditional reference law rather than
unconditional random labeling. These observations motivate a comparison; they
do not establish that the existing result is invalid or guarantee an alternative
will improve it. No new null model or threshold sweep is introduced here.

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
A future comparison should hold junction rule, reference, cohorts, Pgen model,
null sampler and aggregation fixed, then vary exact V, fixed fuzzy V and no V.
Apply each V rule identically to real and mock panels. Preserve raw missing-gene
status and imputation; do not silently treat imputed equality as observed truth.
This comparison is proposed, not a new recommendation or executed experiment.

When can hard V equality be optimal?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let H1 mean target-specific recognition and H0 mean matched background. For
junction comparison d and V relationship g, the full likelihood ratio is
L(d,g)=P(d,g|H1)/P(d,g|H0), wherever the denominator is positive. A hard V
filter is justified without loss only under the strong assumption
P(g is unequal|H1)=0 (and correctly observed V calls). Then every unequal-V
case has zero numerator and can be rejected by an optimal likelihood-ratio
rule. This assumption is not established by the present results. If unequal V
occurs under H1, hard equality can discard real signal. Conversely, dropping V
is lossless only if L(d,g) is independent of g conditional on the retained
features. Neither condition follows from Pgen matching.

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
