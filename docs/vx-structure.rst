.. _vx-structure:

Junction matching and repertoire structure
==========================================

Junction membership
-------------------

Use the full amino-acid junction, including the conserved anchors. For
junctions :math:`x` and :math:`y`, define

.. math::

   H(x,y)=
   \begin{cases}
   \sum_{i=1}^{|x|}\mathbf{1}\{x_i\ne y_i\},& |x|=|y|,\\
   +\infty,& |x|\ne|y|.
   \end{cases}

The inclusive rule :math:`H\le h` counts substitutions. Radius :math:`h=0` tests
exact junction equality, and :math:`h=1` adds single substitutions. Eligibility
can also require the same species and chain, a compatible peptide–MHC
target, and matching normalized V or V and J gene calls. A hard gene
mask removes ineligible comparisons before matching. Handling of
unresolved calls belongs to that mask, and gene-level agreement pools
alleles.

A paired-chain rule requires both chains of a genuine query pair to match
the two chains of one paired-chain reference receptor for the target. A
target-level chain intersection requires a TRA match and a TRB match
for the same query receptor and target, allowing the reference chains
to come from different receptors. The paired-chain rule implies the
intersection when its paired chains are included in the respective
single-chain pools and the chain predicates agree.

When a cutoff reproduces membership
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fix a finite set of query–target opportunities, the reference catalogues
and all eligibility rules. Let :math:`P` be the opportunities called by a
Hamming rule and :math:`Q` those eligible for a call with a finite score :math:`S`.
If assignment is part of the rule, eligibility includes the winning-label
and tie policy. An inclusive score cutoff :math:`\tau` gives the same calls
exactly when

.. math::

   P\subseteq Q,\qquad
   \max_{u\in Q\setminus P} S(u) < \tau \le \min_{u\in P} S(u).

Every positive must reach the cutoff and every eligible negative must
fall strictly below it. These conditions suffice because opportunities
outside :math:`Q` cannot be called. Use :math:`\max\varnothing=-\infty` and
:math:`\min\varnothing=+\infty` to cover empty sets. A tie across the boundary
prevents separation by one cutoff. Apply the extrema across all
opportunities for a shared cutoff, or within each target for
target-specific cutoffs. For an inclusive distance cutoff :math:`D\le r`,
require :math:`P\subseteq Q` and
:math:`\max_{u\in P}D(u)\le r<\min_{u\in Q\setminus P}D(u)`.

Exact and Hamming matching are recovered in an idealized geometry with
full junctions, no indels, a constant positive cost :math:`\lambda` per
substitution and the same hard gene and target masks. On eligible
comparisons, :math:`D(x,y)=\lambda H(x,y)`, so :math:`D\le\lambda h` is exactly
:math:`H\le h`. Taking the minimum distance over eligible references preserves
this equivalence for target membership, with an empty reference set
giving infinite distance. For a paired-chain receptor with the same radius on both
chains, take the maximum of the two chain distances before minimizing
over paired-chain references. Summing the chain distances bounds the total
substitution count. A two-substitution total admits two differences on
one chain, which the per-chain radius-one rule excludes.

In the :doc:`vX research score <vx-formalism>` distance and score terms
can change the ordering through trimming, unequal substitution costs, gap
costs, V-loop costs, source support, the background margin and the second
hit. Its cutoff must meet the separation condition on the stated cohort
to reproduce Hamming membership. Keep exact matches eligible when
comparing against ordinary exact matching.

Per-target degree distributions
-------------------------------

Within one reference target, let :math:`X` be the distinct valid full
junctions. Join distinct junctions at Hamming distance one. For each
:math:`x\in X`, its degree is

.. math::

   d(x)=\sum_{y\in X\setminus\{x\}}\mathbf{1}\{H(x,y)=1\}.

For a nonempty graph, the empirical degree probability mass at integer
:math:`j\ge0` is :math:`|\{x\in X:d(x)=j\}|/|X|`. The mass at zero is the orphan
fraction. Inclusive degree :math:`d(x)+1` counts the vertex itself once.
On distinct junctions, a radius-zero graph has no edges.

The structure command also reports the graph on distinct
``(junction, V gene, J gene)`` clonotypes. Let :math:`m_x` be the
number of these vertices sharing junction :math:`x`. The command pools alleles
to genes and retains unresolved calls as explicit categories. Edge
membership depends on junction distance. Exact-junction vertices form
a clique, and each one-substitution junction edge expands to a complete
bipartite set of clonotype edges. Each of the :math:`m_x` clonotypes has degree

.. math::

   d_{\mathrm{clone}}(x) = (m_x-1) +
      \sum_{y\in X\setminus\{x\}} m_y\,\mathbf{1}\{H(x,y)=1\}.

The clonotype degree distribution counts this value :math:`m_x` times and
divides by :math:`\sum_{x\in X}m_x`. At radius zero only the clique term
remains. For a hypothetical V or V/J edge mask, count the multiplicities
within the allowed gene stratum. The structure command uses the
junction-only edge rule, includes valid junctions with unresolved V
calls and reports the distribution counts and denominators. An empty
graph has no degree distribution. The implementation batches the
junction search and expands degrees without materializing clonotype edges.

For vX, a threshold defines a bipartite query–target graph. Its adjacency
:math:`B_{qe}(t)` equals one for a finite score at least :math:`t`
under the chosen call policy, and zero otherwise. The target degree
:math:`\sum_q B_{qe}(t)` counts called query opportunities for :math:`e`.
The query degree :math:`\sum_e B_{qe}(t)` counts supported targets and is at
most one under unique-winner selection. Specify whether a vertex is a
source record or a distinct receptor key before interpreting either
degree. The fixed competitor and target panels are part of this graph.

Read rarefaction and graph acquisition
--------------------------------------

For a library with :math:`N>0` eligible reads, let :math:`M` be the set of matched
source clonotypes and let :math:`a_i` be the read abundance of source
clonotype :math:`i\in M`. Hold match membership fixed and sample reads
uniformly without replacement at whole-library depth :math:`m\in\{0,\dots,N\}`.
The expected number of matched source clonotypes found is

.. math::

   R(m)=\sum_{i\in M}
   \left[1-\frac{\binom{N-a_i}{m}}{\binom{N}{m}}\right].

The ratio counts samples that avoid all :math:`a_i` reads. Take
:math:`\binom{n}{m}=0` for :math:`m>n\ge0`. Linearity of expectation gives the sum
even when discovery events depend on each other. The subset-hit
probability is the same identity used in Eq. 3 of
`Nipperess and Matsen, §4.1 <https://arxiv.org/html/1208.6552#S4.SS1>`_.
The vdjtools subset rarefaction API pools unmatched reads into one
category for the richness calculation and subtracts the expected
presence of that category. The denominator is the original eligible
library depth. Distinct source clonotypes keep their own read counts
even when they share a matching key. At depth zero the expected
richness is zero, and at depth :math:`N` it is :math:`|M|`.

For uniform sampling of distinct receptor IDs, each ID has unit
abundance. With :math:`B` called IDs among :math:`N` IDs, the expected number
retained is :math:`mB/N`, and the probability of retaining any called ID is
:math:`1-\binom{N-B}{m}/\binom{N}{m}`. Depth here is measured in receptor IDs.

Graph acquisition samples :math:`m\in\{1,\ldots,N\}` distinct vertices
uniformly without replacement from a fixed undirected graph of :math:`N>0`
observed vertices. Retain edges whose endpoints both occur in the
sample. Write :math:`E` for the edge set and :math:`d_i` for the degree of vertex
:math:`i` in the full graph. For :math:`N>1`, both endpoints of an edge are
sampled with probability :math:`m(m-1)/(N(N-1))`. The expected edge count
is :math:`|E|` times this probability, or zero for :math:`N=1`. Conditional on
retaining vertex :math:`i`, it is an orphan when none of its neighbors occur
among the other :math:`m-1` vertices. The expected sampled orphan count is

.. math::

   O(m)=\frac{m}{N}\sum_{i=1}^{N}
   \frac{\binom{N-1-d_i}{m-1}}{\binom{N-1}{m-1}}.

The expected number of sampled vertices with a neighbor is :math:`m-O(m)`,
and their expected fraction is :math:`1-O(m)/m`. A one-vertex sample is
always orphaned. At full depth the formulas return the original edge
and orphan counts. The structure command reduces these expressions
over the full-graph degree histogram.

For a fixed exact or Hamming rule, one can also rarefy reference keys.
Let :math:`K_e(q)` of the :math:`N_e` target keys be eligible neighbors of
query :math:`q`. A uniform sample of :math:`m` keys contains a hit with
probability

.. math::

   1-\frac{\binom{N_e-K_e(q)}{m}}{\binom{N_e}{m}}.

Use paired-chain receptors as the sampling units for a paired-chain rule. The expectation
of the number of recovered queries sums the probabilities over fixed
query opportunities. Under vX, dropping a target key can change both the
nearest and second-nearest evidence. Its reference-depth experiment
re-evaluates the full target score on each sampled target set, with the
other-peptide background fixed.

Epitope assignment by neighbor counts
-------------------------------------

Let :math:`\mathcal R(q)` be the eligible distinct reference keys with finite distance after
applying the chain, gene and benchmark hit-list rules. Let :math:`E(r)`
be the set of peptide–MHC labels on key :math:`r`, and let :math:`D(q,r)`
be a label-independent receptor distance, with smaller values closer.
For vX geometry, use the distance before the target-specific TRB discount. At a fixed
radius :math:`t`, use :math:`\mathcal N_t(q)=\{r\in\mathcal R(q):D(q,r)\le t\}`.
For positive integer :math:`k`, let :math:`t_k(q)` be the distance to the
:math:`\min(k,|\mathcal R(q)|)`-th eligible finite hit and use
:math:`\mathcal N_{t_k(q)}(q)`, including every hit tied at the boundary.
When there are no finite hits, the neighborhood is empty.
An optional maximum radius :math:`t_{\max}` gives the combined
threshold :math:`\min(t_k(q),t_{\max})`, so fewer than :math:`k` hits
may vote.

.. math::

   V_e(q)=\sum_{r\in\mathcal N(q)}\mathbf 1\{e\in E(r)\}.

Each key votes once for each distinct label it carries. Assign the unique
top label only when its vote count meets the declared positive minimum
:math:`v_{\min}`; abstain on a tie or insufficient support. Radius zero
recovers exact neighbors when zero distance is equivalent to equality
under the declared key and eligibility rules. These votes define a
classification rule. A cutoff calibrated for the competitive score
cannot calibrate neighbor votes.

Full-catalogue VDJdb benchmark cohorts
---------------------------------------

Build separate cohorts for human and mouse, TRA, TRB and paired-chain,
and MHC class I and II. In each cohort a target is the full recorded
peptide–MHC identity. Admit it when more than 30 valid source records
support that target in the same lane, before query exclusion. In the
paired-chain lane count genuine paired observations. Count distinct
physical query keys separately from source records and reference keys.

Report both all eligible cases and a capped query view. To construct
the capped view, sort unique physical query keys by SHA-256 of a
compact JSON array containing seed 20261128 followed by the key fields.
Use junction and V call for one chain and both junctions and V calls
for a pair. Keep a case only when each of its eligible annotations
has fewer than 300 selected cases. Retain the case's complete label
set. Fix the reference catalogue in both views and apply exact
exclusion only to the hit list. Keep unsupported cases in matched
evaluation. Other target annotations in VDJdb provide comparison
cases for this catalogue-recovery task.

The original VDJdb paper describes sequence similarity networks.
The subsequent VDJdb motif compendium used at least 30 distinct
amino acid sequences per epitope. The present benchmark uses a
strict source-record gate per full pMHC and lane.

* `Original VDJdb paper <https://pubmed.ncbi.nlm.nih.gov/28977646/>`__
* `VDJdb motif compendium <https://pubmed.ncbi.nlm.nih.gov/31588507/>`__

Paired motif ranking refinement
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let :math:`S_e` be the paired-chain score in the full compatible target
context and :math:`m_e` the smaller of the two fixed-catalogue chain motif
enrichments. For targets with more than two reference keys, add

.. math::

   \lambda_+[m_e]_+-\lambda_-[-m_e]_+.

The two reference keys in the gate correspond to the two pooled
pseudocount observations in the motif score. Keep every target at
the original maximum unchanged. Cap other scores strictly below that
maximum, at its binary64 predecessor. This preserves the winning tier
and any cutoff on the maximum score. It can change lower target ranks
and per-target AP and AUROC. The caller supplies the species weights.
Use ``vdjmatch.evalue.competitive.paired_motif_ranking_scores`` for
retrospective ranking; benchmarking reuses the complete fixed motif
model, including exact junctions in model counts.

On the capped mouse paired-chain cohorts, separately fitting the
two MHC-class motif models and setting support weight to zero and
conflict weight to four improved all 39 AP, AUROC and standardized
partial AUROC endpoints against frozen vX and native TCRdist3.
The setting was selected on the benchmark. The human class-I
full-catalogue comparison is a separate evaluation.

Ranking reactivity within an annotated assay
--------------------------------------------

An assay-local ranking starts with a supplied target annotation for every
receptor. Rank each receptor using its distances to the other receptors
with that annotation. The outcome labels are used only for evaluation.
The reusable reducer ``vdjmatch.evalue.competitive.paired_assay_local_score``
takes a square same-target distance matrix, masks self-hits and accepts an
optional exact-hit mask. It returns the negative weighted mean of the first
:math:`k` distances, with weight :math:`(k-i+1)^p` at rank :math:`i`.

The matched TCRvdb comparison used the accepted human paired-chain
contact geometry, :math:`k=\max(1,\lfloor0.124(n-1)\rfloor)`, and
:math:`p=7/3`. It led the supplied organizer scores on all 24 paired
ranking endpoints, with both self-only and paired-junction exact
exclusion. The setting was selected on the evaluation panel.
It requires the annotated assay group at scoring time and has
no fixed-reference specificity cutoff.

To reduce a prepared matrix, run:

.. code-block:: bash

   vdjmatch paired-assay-score --distances distances.npy --queries queries.tsv --output scores.tsv

Provide one annotated assay group per matrix with the query IDs in the
matrix row and column order. The ``queries.tsv`` file uses the
``query_id`` column. Supply ``--exact-hits exact.npy`` for an
aligned boolean paired-junction mask. The output gives ``query_id`` and
``score`` in input order. Insufficient finite neighbors yield
:math:`-\infty`.

Multi-epitope reciprocal nearest neighbors
------------------------------------------

The vX research benchmark builds one joint graph per species and chain.
A vertex is a distinct full junction of at least five canonical residues
with its resolved V allele. The benchmark deduplicates full source rows
and caps groups by species, chain and peptide. Once the cohort is fixed, labels enter the
evaluation after graph construction and clustering. Every pair of
vertices with equal full junctions is ineligible as a neighbor, even if
their V alleles differ.

The runner uses the weighted pair geometry from the :doc:`vX score
description <vx-formalism>` in integer units of 1/256 of a base
distance unit. For base distance :math:`D_0(i,j)`, component loop costs
:math:`D_\ell(i,j)` and configured weights :math:`w_\ell`, it computes

.. math::

   D_{ij}=256D_0(i,j)+\sum_{\ell\in\{1,2,2.5\}}
       \bigl(\operatorname{int}(256w_\ell)-256\bigr)D_\ell(i,j).

The pair geometry omits source support, the background margin and the
second-hit score. With :math:`n` placed vertices, set :math:`k=\lceil\log n\rceil`,
using the natural logarithm. Every row in the benchmark has at least
:math:`k` eligible comparisons. Let :math:`r_i` be the :math:`k`-th smallest distance
in the eligible part of row :math:`i`. Keep an arc :math:`i\to j` when :math:`j` is
eligible and :math:`D_{ij}\le r_i`. A vertex can have more than :math:`k`
outgoing arcs because all ties at the boundary are kept. For background
on directed nearest-neighbor graphs, see
`Baron and Darling <https://arxiv.org/abs/1908.07645>`_.
The benchmark evaluates exact row distances in batches.

The undirected graph keeps an edge precisely when both arcs occur.
Symmetric pair distances can give different ranks in the two rows.
For symmetric distances an edge requires
:math:`D_{ij}\le\min(r_i,r_j)` after applying the neighbor mask.
A global cutoff can reproduce those edges only if the largest distance
of a retained edge is strictly below the smallest distance of an
eligible omitted edge. The benchmark uses
unweighted Louvain clustering with resolution 1 and seed 0, then maps
communities back to source rows by junction and resolved V allele.
Each unplaced source row receives its own singleton label for evaluation.

The sensitivity analyses can use :math:`k=\lceil s\log n\rceil` with a positive
scale :math:`s` and can vary Louvain resolution. One optional edge weight
uses the distance and the two endpoint radii

.. math::

   w_{ij}^{\mathrm{local}}=
   \begin{cases}
     \exp(-D_{ij}^2/(r_i r_j)), & r_i r_j>0,\\
     1, & r_i r_j\le0.
   \end{cases}

Another optional weight multiplies each edge's current weight by
:math:`1+\eta|\Gamma_i\cap\Gamma_j|/|\Gamma_i\cup\Gamma_j|`, where
:math:`\Gamma_i` is the neighbor set in the reciprocal graph and :math:`\eta\ge0`.
The base weight is one when local weighting is off. Both options
use the reciprocal edge set.

Contact windows on an IMGT alignment
------------------------------------

For a structural contact profile, let :math:`f_j` be the observed contact
frequency at aligned column :math:`j`, with :math:`n_j` resolved strata.
Report peptide and MHC contacts separately, and count a site once in their
union. A support-weighted Gaussian with width :math:`h` in alignment columns
gives

.. math::

   \widetilde f_j(h)=
       \frac{\sum_k n_k f_k\exp[-(j-k)^2/(2h^2)]}
            {\sum_k n_k\exp[-(j-k)^2/(2h^2)]}.

Observed zeros contribute to the denominator. Unobserved columns have
:math:`n_k=0`. A value there is an interpolation.
Use valleys between stable peaks to locate contiguous windows. Retain
all observed positive scoring sites in their windows.

For the positive sites :math:`j\in R` in one window and clipped substitution
costs :math:`\delta_j`, define

.. math::

   A_R=\sum_{j\in R} f_j\delta_j,\qquad
   M_R=\sum_{j\in R}f_j,\qquad
   K_R=\frac{M_R^2}{\sum_{j\in R}f_j^2}.

The mean cost :math:`A_R/M_R` has the same mismatch scale for every window.
The contact mass :math:`A_R` also reflects the number and recurrence of
contacting sites. The effective-width cost :math:`K_RA_R/M_R` has the
units of a substitution sum. Equal site weights make :math:`K_R` the
number of sites. A concentrated profile has a smaller effective width.
Compute the window totals before scoring allele pairs, so their scale
does not depend on gaps in a particular comparison. Report score and
template-coverage checks on the same reference and benchmark populations
for each proposed substitution of a frozen loop term.

Paired receptor contacts
------------------------

For TRA site :math:`i` and TRB site :math:`j`, estimate :math:`p_{ij}`
from structures that resolve both positions. Use one interchain
contact observation per structure and apply the same grouping weights
used for the receptor contact profiles. If :math:`\delta_{\alpha i}`
and :math:`\delta_{\beta j}` are clipped mismatch costs between the query
and one paired-chain reference, the coupled contact cost is

.. math::

   C_{\alpha\beta}=
       \frac{\sum_{i,j}p_{ij}\delta_{\alpha i}\delta_{\beta j}}
            {4\sum_{i,j}p_{ij}}.

The product is zero when either endpoint is unchanged. The divisor
four is the maximum clipped single-site cost, so the coupled term has
the same units and range as one site. It supplements the two chain
distances before the shared paired-chain hit list is ordered. This research
term uses only verified canonical V-domain sites. Give missing
templates the frozen target geometry and retain the hit-list exclusions.

The class-I human derivation used 163 complexes and 54 observed V-domain
pairs across 16 sites per chain. TRA44–TRB44 contacted in 160 of 163
jointly resolved complexes and contributes 9.6% of the total pair weight.
The mouse class-I cohort yielded 31 pairs across ten TRA sites and eleven
TRB sites in 30 complexes. Twenty-nine pairs were observed in both sets.
The other 25 human pairs were jointly resolved in mouse, with zero observed
contacts. Two mouse pairs had zero contacts in the human derivation.
The 26 mouse class-II complexes supplied 25 positive V pairs, 19 of them
shared with mouse class I. Across both mouse classes there are 37
distinct positive V pairs.
Estimate the weights separately by species. The human paired benchmarks
test the human template only.

The human benchmarked extension adds :math:`C_{\alpha\beta}/32` to each
paired-chain target distance. The two target hits are reordered after the
correction. The competitor still uses the frozen native geometry. Of 16
tested amplitudes, :math:`1/32` alone improved TCRvdb ranking while
preserving all matched matrix and reference-depth endpoints. It retained
the 78 strict M1 wins over native TCRdist3.

The generated paired null retained the same cutoff and all 1,000
receptor–target calls among 100,000 queries. Five tie flags and one unique
winner changed outside the called set. The COVID replay retained the same
1,661 whole-catalogue and 314 approved-target receptor calls, with
identical visit and depth probabilities. The false-positive rate is
measured on the generated bank used to fit the cutoff. This extension
applies to the paired research score.

The same amplitude with the 31 mouse class-I contact weights was also
tested on paired-chain mouse VDJdb. At 6,978 distinct receptor cases,
the top labels and ties matched frozen vX, but 53 of 736 target AP and
AUROC endpoints decreased. The mouse scorer therefore retains frozen vX.
The mouse benchmark tests catalogue recovery. It has no matched mouse
null calibration or experimental negative binding panel.
