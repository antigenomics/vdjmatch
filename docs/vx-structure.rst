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

To recover single-chain Hamming membership from the score, compare
full junctions without indels, use a constant cost :math:`\lambda>0`
per substitution, and apply the same hard eligibility mask. Set the
loop, source, competition and second-hit contributions to zero. If
:math:`R_e` is the eligible reference set for target :math:`e`, the
score and membership rule then reduce to

.. math::

   S_e(q)=-\lambda\min_{r\in R_e}H(q,r),\qquad
   S_e(q)\geq-\lambda h\ \Longleftrightarrow\
   \min_{r\in R_e}H(q,r)\leq h.

An empty reference set gives no hit. Apply the rule to each target
independently to retain all supported labels. A winning-label policy
can further reduce those labels.

For a radius :math:`h` applied separately to each chain, let
:math:`R_{\alpha\beta,e}` be the eligible linked keys and
:math:`R_{\alpha,e},R_{\beta,e}` the eligible chain keys of target
:math:`e`. The two matching rules use

.. math::

   \begin{aligned}
   D_{\rm pair}(q,e)&=\min_{r\in R_{\alpha\beta,e}}
      \max\{H(q_\alpha,r_\alpha),H(q_\beta,r_\beta)\},\\
   D_{\rm int}(q,e)&=\max\left\{
      \min_{r\in R_{\alpha,e}}H(q_\alpha,r),
      \min_{r\in R_{\beta,e}}H(q_\beta,r)\right\}.
   \end{aligned}

Apply :math:`D\le h` to either distance. The minimum of an empty
set is :math:`+\infty`. When the chain pools include the corresponding
linked chains under the same masks, :math:`D_{\rm int}\le D_{\rm pair}`.
Summing the chain distances sets a total
substitution budget. A budget of two admits two differences on
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

With :math:`N>0` clonotype vertices, uniform vertex mass and
self-inclusive similarity give :math:`z_i=(d_i+1)/N`. The
reported similarity-aware Hill orders zero, one and two are

.. math::

   H_0=\frac1N\sum_i z_i^{-1},\qquad
   H_1=\exp\!\left[-\frac1N\sum_i\log z_i\right],\qquad
   H_2=\frac{N}{1+\bar d}.

Here :math:`\bar d` is the mean ordinary degree. All three return
:math:`N` for an edgeless graph and one for a complete graph.
The synthetic inclusive degree mass is :math:`\sum_i(d_i+1)=N+2|E|`,
where :math:`E` is the undirected edge set. The degree-based Chao
fields use these synthetic counts. Occurrence-based extrapolation
uses source counts and answers a different question.

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

For :math:`N\geq0` eligible items, a specified subset of
:math:`a\in\{0,\ldots,N\}` items, and a uniform sample of
:math:`m\in\{0,\ldots,N\}` without replacement, define

.. math::

   A(N,a;m)=\frac{\binom{N-a}{m}}{\binom{N}{m}}
           =\frac{\binom{N-m}{a}}{\binom{N}{a}}.

This is the probability that the sample avoids the subset. Use
:math:`\binom{n}{m}=0` for :math:`m>n\geq0`. Its complement
is the probability of finding at least one subset item, as in
`Nipperess and Matsen, Eq. 3 <https://arxiv.org/html/1208.6552#S4.SS1>`_.
At depth zero the avoidance probability is one. At full depth it is
one for an empty subset and zero otherwise.

For a library with :math:`N>0` eligible reads, let :math:`M` be the set of matched
source clonotypes and let :math:`a_i` be the read abundance of source
clonotype :math:`i\in M`. Hold match membership fixed and sample reads
uniformly without replacement at whole-library depth :math:`m\in\{0,\dots,N\}`.
The expected number of matched source clonotypes found is

.. math::

   R(m)=\sum_{i\in M}\bigl[1-A(N,a_i;m)\bigr].

Linearity of expectation gives this sum even when discovery events
depend on each other.
The vdjtools subset rarefaction API pools unmatched reads into one
category for the richness calculation and subtracts the expected
presence of that category. The denominator is the original eligible
library depth. Distinct source clonotypes keep their own read counts
even when they share a matching key. At depth zero the expected
richness is zero. At full depth it is :math:`|M|` because each
eligible source clonotype has positive integer read abundance. The
common-depth CLI grid ends at the smallest eligible library depth.

For uniform sampling of distinct receptor IDs, each ID has unit
abundance. With :math:`B` called IDs among :math:`N` IDs, the expected number
retained is :math:`mB/N`, and the probability of retaining any called ID is
:math:`1-A(N,B;m)`. Depth here is measured in receptor IDs.

Graph acquisition samples :math:`m\in\{1,\ldots,N\}` distinct vertices
uniformly without replacement from a fixed undirected graph of :math:`N>0`
observed vertices. The structure command uses the junction/V/J clonotype
graph defined above. Retain edges whose endpoints both occur in the
sample. Write :math:`E` for the edge set and :math:`d_i` for the degree of vertex
:math:`i` in the full graph. For :math:`N>1`, both endpoints of an edge are
sampled with probability :math:`m(m-1)/(N(N-1))`. The expected edge count
is :math:`|E|` times this probability, or zero for :math:`N=1`. Conditional on
retaining vertex :math:`i`, it is an orphan when none of its neighbors occur
among the other :math:`m-1` vertices. The expected sampled orphan count is

.. math::

   O(m)=\frac{m}{N}\sum_{i=1}^{N}
   A(N-1,d_i;m-1).

The expected number of sampled vertices with a neighbor is :math:`m-O(m)`,
and their expected fraction is :math:`1-O(m)/m`. A one-vertex sample is
always orphaned. At full depth the formulas return the original edge
and orphan counts. The structure command reduces these expressions
over the full-graph degree histogram. If the proportion of
vertices of degree :math:`j` is :math:`p_j`, the expected
orphan fraction simplifies to :math:`O(m)/m=\sum_j p_j A(N-1,j;m-1)`.
This samples the observed graph with its original distances and edges.

For a fixed exact or Hamming rule, one can also rarefy reference keys.
Let :math:`K_e(q)` of the :math:`N_e` target keys be eligible neighbors of
query :math:`q`. A uniform sample of :math:`m` keys, with
:math:`0\leq m\leq N_e`, contains a hit with
probability

.. math::

   1-A(N_e,K_e(q);m).

The probability is zero at depth zero. At full depth it is one
if at least one eligible key exists and zero otherwise.

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
the capped view, retain unique physical query keys with at least one
admitted label. Put those keys in canonical order, then shuffle with
Python's ``random.Random(20261128)``.
Use junction and V call for one chain and both junctions and V calls
for a pair. Keep a case only when each of its eligible annotations
has fewer than 300 selected cases. Retain the case's complete label
set. Fix the reference catalogue in both views and apply exact
exclusion only to the hit list. Keep unsupported cases in matched
evaluation. Other target annotations in VDJdb provide comparison
cases for this catalogue-recovery task.

Apply this eligibility cutoff and cap only to full-VDJdb benchmarking.
Use TCRvdb, NLV and the LLL–LLW–ELA assay as supplied. Epitope
structure classification, motif counts, rarefaction, neighbourhood
graphs and mismatch-aware functional diversity use the full depth.

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
enrichments. Replace a missing chain motif by zero before taking the
minimum. This prevents a positive bonus from one chain alone and
still lets a negative motif penalize the pair. For targets with more
than two reference keys, add

.. math::

   \Delta_e=\lambda_+[m_e]_+-\lambda_-[-m_e]_+.

Here :math:`[x]_+=\max(x,0)`. A negative enrichment means that
the observed sequence windows occur less often in that target
than in the pooled reference. With support weight :math:`\lambda_+=0`,
the rule subtracts :math:`\lambda_-|m_e|` in that case and otherwise
adds zero. The conflict weight :math:`\lambda_-` converts one unit of
natural-log enrichment into distance score units. The reported
weights were selected on their stated comparison panels.

The gate requires more observed keys than the two prior observations
used for motif shrinkage at either window length. Set
:math:`\Delta_e=0` when the gate fails. Keep every target at the
original maximum unchanged, as well as every original no-hit target.
For other targets let :math:`L_e` be the supplied finite lower-tier
baseline, falling back to :math:`S_e` when it is unavailable. The
refined value is :math:`\min(L_e+\Delta_e,S_*^-)`, where
:math:`S_*^-` is the binary64 predecessor of the query's original
maximum. The cap preserves the winning tier and any cutoff on that
maximum. A separate lower-tier baseline can still change a rank when
:math:`\Delta_e=0`. The caller supplies the species weights.
Use ``vdjmatch.evalue.competitive.paired_motif_ranking_scores`` for
retrospective ranking; benchmarking reuses the complete fixed motif
model, including exact junctions in model counts.

For a single observed chain, use its motif enrichment directly.
The optional compatible-sibling prior uses the same peptide under
compatible alternative recorded MHC restrictions. Keep distinct
physical keys and remove keys shared with the focal target from its
sibling pool. For motif token :math:`t` and target :math:`e`, use

.. math::

   p_{et} = \frac{2}{N_s+2}\frac{1+c_t}{2+N}
          + \frac{N_s}{N_s+2}\frac{1+c_{st}}{2+N_s}.

Here :math:`N` and :math:`N_s` count keys in the fixed global and
sibling pools, and :math:`c_t` and :math:`c_{st}` count keys carrying
the token in each. Both pool frequencies include a Beta(1,1) prior.
The outer coefficients mix the two smoothed frequencies, with
the sibling weight increasing with its key count. The enrichment for a token is
:math:`\log[(c_{et}+2p_{et})/((N_e+2)p_{et})]`, where :math:`N_e`
counts focal reference keys and :math:`c_{et}` counts those carrying
the token. Average over the query's available anchored tokens. With
zero sibling keys the prior reduces to the global estimate. This
borrowing preserves distinct peptide--MHC labels.

On the capped mouse paired-chain cohorts, separately fitting the
two MHC-class motif models and setting support weight to zero and
conflict weight to four improved all 39 AP, AUROC and standardized
partial AUROC endpoints against frozen vX and the native distance comparator.
The setting was selected on the benchmark. The human class-I
full-catalogue comparison is a separate evaluation.

Replay the verified mouse refinement from prepared scores:

.. code-block:: bash

   vdjmatch paired-motif-rank --scores context.npz --queries queries.tsv \
       --targets targets.tsv --support-weight 0 --conflict-weight 4 --output ranks.tsv

The NumPy archive contains the original paired scores ``base``,
the two chain motif score matrices ``tra`` and ``trb``, and the
one-dimensional ``target_key_counts`` array for the same whole
context. The TSV columns ``query_id`` and ``target`` identify the
matrix rows and columns. If the motif preparation has already taken
the minimum across chains, pass that matrix in both chain slots.
Keep all targets that compete for a call in the score matrix, then
select reporting targets after refinement. The mouse ranking
benchmark used one MHC class per run. A search across MHC classes
uses the combined target axis to determine the winning tier.

Optionally include a ``lower_tier_base`` matrix in the archive.
It supplies the baseline for the lower ranks only. A missing entry
falls back to the corresponding ``base`` score, and original no-hits
are preserved. Use the direct addition of the motif adjustment to
this matrix. Reconstructing the same value as a difference from
``base`` can alter floating-point ties.

Using the unweighted nearest-hit matrix and conflict weight four
improved all 117 TRB and all 57 paired-chain human MHC-II endpoints
in the full, capped and common-support scopes. This was a separate
benchmark-selected setting. It did not pass the human MHC-II TRA
or the capped mouse paired comparison.

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
contact geometry, :math:`k=\max(1,\lfloor31(n-1)/250\rfloor)`, and
:math:`p=7/3`. It led the supplied organizer scores on all 24 paired
ranking endpoints, with both self-only and paired-junction exact
exclusion. The setting was selected on the evaluation panel.
It requires the annotated assay group at scoring time and has
no fixed-reference specificity cutoff.
With square weights and the same width, four endpoints fall below
the organizer scores and 18 fall below the current local score.
Taking one eighth of the available neighbours with square weights
raises the organizer losses to five.
Exponential and logarithmic tapers also lower priority endpoints.
At the original width, :math:`r^2\log(1+r/k)`, with descending rank
weight :math:`r=k-i+1`, lowers six endpoints relative to the current
score and one against the organizer under both exclusion policies.

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
At :math:`k=1`, the graph retains all tied mutual nearest
neighbors. An exact strictly increasing transform of distance preserves
the unweighted reciprocal graph because it preserves ranks and ties.

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

Radius graphs at the reciprocal edge budget
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let :math:`B` be the number of reciprocal edges and let :math:`d_*`
be the first observed distance at which the inclusive count of eligible
unordered pairs exceeds :math:`B`. Write :math:`m` for the number of
pairs closer than :math:`d_*` and :math:`t` for the number at that
distance. The hard-radius graph keeps the :math:`m` closer pairs.
For the fractional boundary graph, each eligible pair receives weight

.. math::

   w_{ij} =
   \begin{cases}
     1, & D_{ij}<d_*,\\
     (B-m)/t, & D_{ij}=d_*,\\
     0, & D_{ij}>d_*.
   \end{cases}

Thus all boundary ties have equal weight and the graph has total
undirected edge weight :math:`B`. The archived benchmark requires
an over-budget distance. If the full eligible graph fits the budget,
the mathematical rule reduces to unit weights on all eligible pairs.
Apply weighted Louvain with the same seed
and resolution. An isolate variant of the hard-radius graph adds the
nearest tied reciprocal neighbors of vertices isolated by that radius.

On the locked human and mouse single-chain cohorts, the hard-radius
and fractional graphs raised adjusted Rand and V-measure over the
original native reciprocal graph in all four lanes. Relative to the
hard-radius vX graph, fractional boundary weights increased same-label
pair recall in all four lanes. Comparisons against native geometry
using the same graph rule showed a mix of gains and losses. These are
clustering comparisons using the historical peptide labels.

The graph can combine reciprocal adjacency :math:`R` with fractional
radius weights :math:`F` through :math:`W_\theta=(1-\theta)R+\theta F`.
Both graphs have undirected weight :math:`B`, so the blend has that
same total. The radius contribution is equal for all boundary ties.
An edge's reciprocal membership can change its final weight. All
three tested fractions :math:`\theta\in\{1/4,1/2,3/4\}` increased
same-label pair recall relative to fractional radius alone across all
four frozen vX lanes. Every blend lowered at least one adjusted Rand
or V-measure endpoint relative to the original native reciprocal graph.

Contact windows on an IMGT alignment
------------------------------------

For a structural contact profile, let :math:`f_j` be the observed contact
frequency at aligned column :math:`j`, with support :math:`n_j`
measured in resolved strata.
Report peptide and MHC contacts separately, and count a site once in their
union. A support-weighted Gaussian with width :math:`h` in alignment columns
gives

.. math::

   \widetilde f_j(h)=
       \frac{\sum_k n_k f_k\exp[-(j-k)^2/(2h^2)]}
            {\sum_k n_k\exp[-(j-k)^2/(2h^2)]}.

Observed zeros contribute to the denominator. Unobserved columns have
:math:`n_k=0`. A value there is an interpolation. The smoothed frequency
is undefined if there is no resolved support anywhere in the profile.
Use valleys between stable peaks to locate contiguous windows. Retain
all observed positive scoring sites in their windows.

For the positive sites :math:`j\in R` in one window and clipped substitution
costs :math:`\delta_j`, define

.. math::

   A_R=\sum_{j\in R} f_j\delta_j,\qquad
   M_R=\sum_{j\in R}f_j,\qquad
   K_R=\frac{M_R^2}{\sum_{j\in R}f_j^2}.

For :math:`M_R>0`, the mean cost :math:`A_R/M_R` has the same
mismatch scale for every window. An empty window has no region
contribution to estimate.
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
four is the maximum clipped single-site cost, so
:math:`0\le C_{\alpha\beta}\le 4` for positive total contact weight.
The coupled term has the same units and range as one site. It supplements the two chain
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

The frozen distance has zero contact amplitude. The human benchmarked
extension adds :math:`C_{\alpha\beta}/32` to each
paired-chain target distance, at most :math:`1/8` distance unit per key.
The two target hits are reordered after the
correction. The competitor still uses the frozen native geometry. Of 16
tested amplitudes, :math:`1/32` alone improved TCRvdb ranking while
preserving all matched matrix and reference-depth endpoints. It retained
the 78 strict M1 wins over the native distance comparator.

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
