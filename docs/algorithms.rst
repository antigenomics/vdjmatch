How VDJmatch scores a receptor
==============================

All routes start with anchored junctions, chain identities and V genes. They
retain reference peptide–MHC labels and query identity. Ordinary annotation
keeps exact evidence; a benchmark filters exact query hits while keeping the reference model fixed.

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
