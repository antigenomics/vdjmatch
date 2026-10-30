How VDJmatch scores a receptor
==============================

All routes start with anchored junctions, chain identities and V genes. They
retain reference peptide–MHC labels and query identity. Ordinary annotation
keeps exact evidence; a benchmark filters exact query hits while keeping the reference model fixed.

Fast reference classification
-----------------------------

Index preparation classifies reference structure, separately by species/chain
and pooled epitope or observed full pMHC. It uses reference labels, never query
assay outcomes. Eligibility is at least 30 valid source observations; counts of
distinct junction/V keys remain separate. All observations and pair links stay
in the saved model. These descriptive classes do not establish binding or select
a scoring default.

Let the population contain :math:`K` distinct full-junction/allele-pooled-V keys.
Each key retains its complete binary target-membership vector, including shared
labels. For a target with :math:`N` keys, let :math:`n_v` count its V-gene keys and
:math:`q_v=K_v/K` be the population V distribution. The gene statistic, in bits
per target key, is

.. math::

   G=\sum_{v:n_v>0}\frac{n_v}{N}\log_2\frac{n_v/N}{q_v}.

The background includes the target, as in the fixed-cohort label-shuffle
comparison. Every observed category consequently has population support; no
pseudocount or zero-frequency extrapolation is needed. Unresolved V calls
retain their mass and prevent a biological class assignment.

Tokenize the junction after trimming three residues at each end, retaining each
exact amino-acid 3-mer once per key. For motif :math:`t` within V gene :math:`v`,
let :math:`a_{vt}` count target presence and :math:`A_{vt}` population presence.
Write :math:`p_{vt}=a_{vt}/n_v` and :math:`q_{vt}=A_{vt}/K_v`. The sequence
statistic is the strongest enriched binary projection:

.. math::

   M=\max_{v,t:p_{vt}>q_{vt}}
      \frac{n_v}{N}\left[
        p_{vt}\log_2\frac{p_{vt}}{q_{vt}}+
        (1-p_{vt})\log_2\frac{1-p_{vt}}{1-q_{vt}}\right].

An empty maximum is zero; :math:`0\log 0=0`. Motif depletion does not supply
positive motif evidence. The V weight expresses one projection in bits per
target key; it does not sum correlated, overlapping motifs. Stable lexical
V/motif order resolves exact ties. The selected motif and its four presence
counts remain visible.

Two independent fixed random streams (seed 0) generate 199 shuffles. For
:math:`G`, permute complete membership vectors over population keys. For
:math:`M`, permute those vectors only within V genes. Target sizes and all
shared-membership relationships remain fixed; the latter null also fixes every
target's V composition. Recompute the same maximum over motifs in each
replicate, so choosing the best motif is part of the null statistic.

For either observed statistic :math:`S` and shuffled values :math:`S_b`, report

.. math::

   \Delta S=S-\frac{1}{199}\sum_{b=1}^{199}S_b,\qquad
   p_S=\frac{1+\sum_{b=1}^{199}\mathbf1\{S_b\ge S\}}{200}.

Support requires :math:`\Delta S>0` and :math:`p_S\le0.005`; with 199 draws,
this requires exceeding every shuffled value. Ties count against support.
Under the respective exchangeability nulls these Monte Carlo tails are
conservative; allocating 0.005 to each of two axes bounds the per-profile
probability of falsely supporting either true-null axis by 0.01, without an
independence assumption. This does not control errors across every database
profile, remove study/ascertainment bias or calibrate selected query calls.
Effects and null maxima are shown even when support is absent.

The two support flags define V-and-motif, V-only, motif-only and no-supported-axis
leaves. Missing V is assigned reproducibly from the same population's resolved
gene frequencies before label shuffling, independently of target membership;
original key identities and missingness counts remain fixed. With no resolved
V, one neutral missing category carries no V-usage signal. Missing background
and inadequate distinct support remain unclassified. A small effect measured precisely can pass its noise test while
remaining weak for prediction. Neighbour coverage and edit-one pair-density
enrichment are separate structural descriptors; pair counts share vertices and
are not independent binomial observations.

Implementation reuses the core tokenizer, fixed sparse key/target membership and
key/V-motif presence matrices. Identical target memberships are reduced once.
Shuffles perform count reductions only; no sequence search, SW, Pgen calculation
or reference-model refitting occurs inside them. The completed class table and
native search index are durable prepared-reference assets; query predictions
are not reused.

Current kernel evidence
-----------------------

.. graphviz::

   digraph kernel {
     graph [rankdir=TB, bgcolor="white"];
     node [shape=box, style=rounded, fontname="Helvetica", fontsize=11, color="#0072B2"];
     edge [color="#555555"];
     input [label="Query TRA, TRB or linked alpha/beta\n+ labelled VDJdb reference"];
     distance [label="BLOSUM62 junction and V-loop distance\n+ joint terminal stem (depth 4)"];
     counts [label="One batched search\nTarget / other-peptide counts by distance"];
     score [label="Distance-decaying reference support\nSame-V weight 1; cross-V weight 0.25\nDivide by background counts and target size"];
     select [label="Maximum restriction per peptide\nRank peptides; retain exact ties"];
     pair [label="Paired: linked score + alpha + beta\nKeep unlinked same-target marginal support"];
     output [label="kernel_score and competing score\nSelected labels and support status"];
     input -> distance -> counts -> score -> select -> output;
     select -> pair -> output;
   }

Larger ``kernel_score`` means stronger evidence. A linked score compares both
chains against the same linked reference receptor. Separate alpha and beta
scores also use unlinked reference observations. The score is not itself a
P-value. :doc:`scoring` gives the mathematical definitions of retained variants.

Conditional-rank evidence
-------------------------

.. graphviz::

   digraph rank {
     graph [rankdir=TB, bgcolor="white"];
     node [shape=box, style=rounded, fontname="Helvetica", fontsize=11, color="#009E73"];
     query [label="Query and target / background receptors\nSame distance and exclusion rule"];
     ranks [label="Rank target distances against background\nAverage tied reciprocal ranks"];
     censor [label="Retain all reference receptors\nRepresent beyond-radius ranks by their mean"];
     nearest [label="If target has no within-radius support:\nrefine with exhaustive nearest-reference rank"];
     result [label="conditional_rank_bayes_factor\nHigher means stronger relative evidence"];
     query -> ranks -> censor -> nearest -> result;
   }

This retained alternative measures distance ranks relative to the declared
background. Component probability bounds belong to their original statistics;
they do not calibrate the conditional refinement or a selected antigen label.

Legacy scoring routes
---------------------

.. graphviz::

   digraph legacy {
     graph [rankdir=TB, bgcolor="white"];
     node [shape=box, style=rounded, fontname="Helvetica", fontsize=11, color="#D55E00"];
     input [label="Query junction + V/J\nReference receptors + declared background"];
     search [label="Legacy sequence-neighbour search"];
     original [label="historical-density\nPosition-weighted substitution penalty\nV similarity and background-normalised density"];
     originalout [label="Historical density / log ranking score"];
     match [label="match\nAligned neighbour evidence\nCandidate aggregation and ambiguity"];
     matchout [label="Legacy candidate ranking\nFixed-ball expected count and tail (when available)"];
     input -> search;
     search -> original -> originalout;
     search -> match -> matchout;
   }

``historical-density`` reproduces the retained density model. ``match`` is a
separate older annotation interface. Neither should be described as the current
``annotate`` default, and their probability fields cannot be transferred to the
kernel score.

From ranking to a match decision
--------------------------------

.. graphviz::

   digraph calibration {
     graph [rankdir=TB, bgcolor="white"];
     node [shape=box, style=rounded, fontname="Helvetica", fontsize=11, color="#555555"];
     nulls [label="Matched null queries\nSame reference, scorer and complete label selection"];
     freeze [label="calibrate-calls\nFreeze operating cutoff + producer identity"];
     query [label="Independent query selected scores\nSame complete pipeline"];
     apply [label="apply-call-calibration --null-cases\nTie-aware empirical tail + frozen cutoff"];
     out [label="Ranking score + p_null\ncall_status: accepted / no_supported_match"];
     nulls -> freeze -> apply;
     query -> apply -> out;
   }

The null asks whether an unrelated query would achieve at least this selected
score. All searched labels, ties and unsupported null queries belong in that
calculation. Paired any-axis calls require a null for the maximum over linked,
alpha and beta evidence. See :doc:`explanation` for the output fields.

Reference structure for scoring weights
---------------------------------------

.. graphviz::

   digraph structure {
     graph [rankdir=TB, bgcolor="white"];
     node [shape=box, style=rounded, fontname="Helvetica", fontsize=11, color="#0072B2"];
     input [label="Reference full junctions, V genes and pMHC\nDeclared species, chain and population policy"];
     search [label="One native Hamming-one batch\nUnique full junctions; no indels"];
     target [label="Partition edges by target membership\nPreserve source records and V annotations"];
     connectivity [label="Components and articulation fragments\nExact full-junction removal profiles"];
     usage [label="Junction-length and V-usage composition\nRemove the same junction contributions"];
     tables [label="Full-reference feature tables\nGraph morphology, length and V diversity"];
     future [label="Proposed reference-to-weight rule\nGraph, length and V usage; chainwise",style="rounded,dashed",color="#777777"];
     scores [label="Original, gap-block and stem4 scores\nMatching component null populations"];
     tails [label="Conservative empirical null tails\nCommon log-evidence scale (nats)"];
     sum [label="Explicit weights: nonnegative, sum to one\nCombined ranking and channel contributions"];
     calibration [label="Required final-selection calibration\nAll candidate labels and paired axes",style="rounded,dashed",color="#777777"];
     input -> search -> target -> connectivity -> tables;
     target -> usage -> tables;
     tables -> future [style=dashed];
     future -> sum [style=dashed];
     scores -> tails -> sum;
     sum -> calibration [style=dashed];
   }

``reference-structure`` exports morphology and exact exclusion effects;
``combine-evidence --three-channel`` normalizes the components and applies explicit
weights. Solid stages are implemented. Dashed stages remain the scoring
investigation: the reference-to-weight rule and final selection calibration need
verification before adoption. Neither command reads query binding labels.

Database presence and a label-independent P-value
-------------------------------------------------

.. graphviz::

   digraph database {
     graph [rankdir=TB, bgcolor="white"];
     node [shape=box, style=rounded, fontname="Helvetica", fontsize=11, color="#0072B2"];
     reference [label="All valid junctions of one species/chain\nPool specificity labels; retain exact records"];
     rows [label="Query receptors and declared null draws\nKeep repeated draws and original query identities"];
     join [shape=point,label="",width=.04,color="#555555"];
     search [label="One native nearest-junction batch\nUnit Levenshtein distance, maximum edits fixed"];
     lookup [label="Exact junction and V/J-gene presence\nBoolean lookups; separate from significance"];
     tail [label="P_database = max source (k+1)/(N+1)\nk: null draws at least as close; N: all null draws"];
     output [label="Presence flags + nearest edits + P_database\nNo antigen label is selected"];
     reference -> join; rows -> join;
     join -> search; join -> lookup;
     search -> tail -> output; lookup -> output;
   }

``database-evidence`` tests unusual junction resemblance to any database entry
under the supplied null law. It retains exact matches and includes the full
database search in the null statistic. Its P-value ignores V/J and cannot be
transferred to epitope or paired recognition; the exact gene-clonotype flag is a
separate lookup. The chosen scoring pipeline's ``p_null`` answers a different
question about its selected evidence score.

Single-control directional calibration
--------------------------------------

The preceding ``--null-sample`` graph describes the legacy diagnostic.
The on-demand ``--control`` workflow computes both conditional directions:

.. graphviz::

   digraph directional_control {
     graph [rankdir=LR,bgcolor="white"];
     node [shape=box,style="rounded",fontname="Helvetica",fontsize=10];
     control [label="One raw control\nAA collapse retains NT event counts"];
     reference [label="Selected reference\nEvidence filters / custom table"];
     pgen [label="Fresh native Pgen\nSpecies / locus / 1-bit strata"];
     mockr [label="Mock reference event slots\nPreserve N and per-pMHC composition"];
     mockq [label="Same-Pgen query draws\nIndependent TRA/TRB for paired null"];
     search [label="Same search event A\nScope / score / selection"];
     pr [label="P_reference\nFixed query / random reference"];
     pq [label="P_query\nRandom query / fixed reference"];
     combined [label="P_both = max(P_reference, P_query)\nNo independence assumption"];
     control -> pgen; reference -> pgen;
     pgen -> mockr; pgen -> mockq;
     mockr -> search; mockq -> search; reference -> search;
     search -> pr; search -> pq; pr -> combined; pq -> combined;
   }

The scoped nearest/count adapter and independent-chain paired sampler are
implemented. Full VDJmatch-score, selected-label and paired search adapters
must repeat their complete statistic, including reference-derived weights,
inside the null calculation. See :doc:`control-calibration` for the equations,
event/unique-key distinction and dependence-safe maximum proof.

Repertoire enrichment with mock target corpora
----------------------------------------------

M3's reference-null test fixes a sample and counts its distinct receptor identities
supporting each pMHC. A receptor matches when an available chain satisfies the
specified full-junction Hamming radius and gene predicate. Chain union is evaluated
inside every real and mock calculation, so paired receptors cannot be counted twice.
The sample index is execution-owned; each real or mock target population is searched
in one native batch per locus with bounded primitive hit buffers.

For target ``e``, let ``K_e`` be the observed matched-receptor count, ``B`` the number
of fresh mock target corpora, and ``K_eb`` the count for mock ``b``. The one-sided
empirical tail is ``(1 + sum_b I[K_eb >= K_e]) / (B + 1)``. Ties and zero-hit mocks
are retained. The expected count is the mean of ``K_eb``; excess is observed minus
expected. The same calculations on summed receptor abundance are descriptive
secondary endpoints. Biological donor replication is distinct from mock replication.

The construction matches reference publicity while
removing selected target identities. It uses native marginal full-junction Pgen,
not the historical snippet's conditional probability and cohort-fitted V/J factor.
Binning is ``round(log2(Pgen + 2**-60))``. Explicit options preserve V/J composition
or require different V/junction identities; these change the null hypothesis and
must not be selected by suppressing unfavourable results. See :doc:`cli` for the
sampling laws, significance families and complete output contracts.
