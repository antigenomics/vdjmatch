How VDJmatch scores a receptor
==============================

All routes start with anchored junctions, chain identities and V genes. They
retain reference peptide–MHC labels and query identity. Ordinary annotation
keeps exact evidence; a benchmark explicitly removes the query clonotype.

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
     tables [label="Full and punctured feature tables\nGraph morphology, length and V diversity"];
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
