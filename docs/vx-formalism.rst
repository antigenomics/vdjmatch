Frozen vX research score
========================

This research candidate uses a fixed receptor catalogue. Let :math:`q` be a
TRA receptor, a TRB receptor, or a paired-chain TRA–TRB receptor and let :math:`e` be a
candidate peptide–MHC. A chain is represented by its anchor-inclusive amino
acid junction and mapped V allele. Search keys are distinct encodings. A
key can have several source observations. Paired keys represent observed
links and both chain distances are evaluated against the same key.

For one chain, let :math:`B(a,b)` be the BLOSUM62 substitution score.
The native comparison charges :math:`\delta(a,b)=\min(4,[4-B(a,b)]_+)`
for canonical residues, with :math:`[x]_+=\max(0,x)`. It discards three
junction residues at the N terminus and two at the C terminus from the
substitution sum. For junctions of different lengths, the native alignment
places one contiguous gap block and charges four units per gap residue.
The junction's weight of three makes the maximum cost of one core
substitution and the cost of one gap residue both 12 native distance units.
Let :math:`d` be the full junction length difference and
:math:`\Pi(q,r)` the permitted alignments of their retained residues.
At equal lengths there is one ungapped alignment. With unequal lengths,
the native routine tries internal positions for the length difference as
one block, with the position window adjusted for short junctions. For
the aligned index pairs :math:`(i,j)\in\pi`, the junction cost is

.. math::

   D_J(q,r)=4d+\min_{\pi\in\Pi(q,r)}
      \sum_{(i,j)\in\pi}\delta(q_i,r_j).

Let :math:`D_1,D_2,D_{2.5}` be the aligned V-loop costs. The frozen per-chain
geometry is

.. math::

   D_w = 3D_J + D_1 + \frac34 D_2 + \frac32 D_{2.5}
             + \frac{D_1+D_2+D_{2.5}}{256}.

CDR2.5 is the selected
short segment within FR3, and full FR1–3 contributions are absent. In
code the native distance has unit loop weights, and vX adds the three
weight differences to that native distance in the original order.
The fractions encode the selected geometry exactly. The simpler base weights
:math:`(1,3/4,3/2)` lowered 23 of the 72 frozen ranking endpoints on the matched comparison.

The frozen geometry minimizes a sum of substitution and gap costs.
The clipped BLOSUM62 cost :math:`\delta` is defined above. A native
kernel that accepts similarities can represent this cost as
:math:`s(a,b)=-\delta(a,b)`. Its diagonal is zero, so the kernel's
conversion :math:`s(a,a)+s(b,b)-2s(a,b)` yields :math:`2\delta(a,b)`.
The frozen junction code divides out that factor before applying the
junction weight of three. Gap costs must follow the same unit conversion.

The experimental Smith–Waterman and Needleman–Wunsch arms instead
align using raw BLOSUM62 similarities and affine gap penalties. Let
:math:`H(x,y)` be the resulting optimal similarity for a segment,
using the same mode and penalties for self-alignment and cross-alignment.
Convert it once to distance units through

.. math::

   D_H(x,y)=\frac{H(x,x)+H(y,y)-2H(x,y)}{2}.

The native alignment replay returns the numerator. The two self-scores
set the similarity scale for the two sequences. For a fixed query and
reference, they are constant across candidate alignments. Therefore
maximizing :math:`H(x,y)` minimizes :math:`D_H(x,y)` for that pair.
Applying a clipped residue cost before alignment defines a different
objective and can select different alignments. Local and global modes also
optimize different alignment sets. Treat the resulting distance as a
dissimilarity unless its metric properties are established for the
chosen mode. Apply segment weights after the conversion to distance.

For target key :math:`r`, order hits by
:math:`d(q,r;e)=D_w(q,r)-a\log_2\max(1,n_{er})`, where :math:`n_{er}` counts
distinct reference IDs supporting that target key. Only TRB uses a nonzero
:math:`a`, a discount of :math:`9/32` units per doubling of support. The
replay takes the maximum distinct-ID count among compatible raw keys
mapping to a native key. The generated-null catalogue counts distinct IDs
across the compatible rows for a native key. The benchmark uses the native
distance when the loop correction is unavailable. The code accumulates
loop corrections before adding them to the native distance.

When both chain corrections are available, for a pair :math:`q=(q_\alpha,q_\beta)` and a paired-chain key
:math:`r=(r_\alpha,r_\beta)`, use

.. math::

   d(q,r;e)=D_w(q_\alpha,r_\alpha)+D_w(q_\beta,r_\beta).

Order these joint distances across paired-chain keys for the same target. Compute
the native competitor distances by summing the two native distances to each
paired-chain competitor. A missing correction on either chain leaves the native
pair distance in the replay. Apply the pair margin and gap adjustment once, after
ordering. For example, two paired-chain keys with chain distances :math:`(2,14)`
and :math:`(9,3)` give joint hit distances 12 and 16. Taking the minimum of
each chain separately would give 5, a distance to neither paired-chain key.


Let :math:`d_1\le d_2` be the two smallest target hit distances, including
equal distances from distinct keys. A singleton repeats :math:`d_1`.
Let :math:`c` be the :math:`k`-th smallest **native** distance in the pooled
background of compatible other-peptide keys. Use the last position if the
fixed background has fewer than :math:`k` keys, and :math:`+\infty` for an
empty background. Under benchmark exact exclusion, read masked hit distances
as :math:`+\infty` without changing model counts or rank positions. The exact
mask compares anchored junctions, requiring equality on both chains for a
pair, and applies to target and background hits. The generated-null run
retains exact hits.

Let :math:`L` count observed chains and let :math:`C=3\cdot4=12`
be the maximum weighted junction cost of a single substitution or
gap residue. The gap cap is :math:`g_\star=2LC`. Define

.. math::

   \begin{aligned}
   X &= [d_1+Lm-c]_+,\\
   G &= \min\{g_\star,[(d_1+d_2)/2-d_1]_+\},\\
   S(q,e) &= -d_1-\gamma X+bG.
   \end{aligned}

The leading term :math:`-d_1` rewards the best hit. The per-chain margin
:math:`m` sets its desired lead over the background, and :math:`\gamma`
sets the cost of falling short. The order :math:`k` sets the background
comparison rank. The gap weight :math:`b` penalizes separation for TRA and
paired-chain receptors and rewards it slightly for TRB. All distance terms have native
distance units. The weights :math:`\gamma,b` are dimensionless.

The gap contribution is in :math:`[-2.4,0]` for TRA, :math:`[0,0.15]` for
TRB and :math:`[-4.8,0]` for paired-chain receptors. In real arithmetic the midpoint gap
is half the separation :math:`d_2-d_1`, so saturation begins at
:math:`d_2-d_1=4LC`. The cap uses four maximum-cost substitutions or
gap residues per chain as the separation scale. For example, an
equal-length junction of length :math:`\ell` can contribute up to
:math:`12[\ell-5]_+` before V-loop terms, and other terms or unequal lengths
change the range further. In the matched catalogue ranking queries
(2,636 records per axis), median full junction lengths were 13 for TRA
and 14 for TRB. Their median core lengths of eight and nine give loose
substitution-only midpoint scales of 48 and 54 units. A length-based cap
would set :math:`g_\star(q)=6\sum_i[\ell_i-5]_+`. Here :math:`\ell_i`
is the full junction length of observed chain :math:`i`. Neither this
capacity nor the cost matrix determines how far apart the two best
reference hits are in a real panel. The original support-gap sweep fixed the
cap across its variants.
Its value remains a design choice. In a matched 78-endpoint ranking replay,
halving the cap reduced 32 endpoints. Doubling it or using the query-core
capacity each reduced six. The frozen cap preserved all 78 endpoints.
Binary64 evaluation forms the midpoint
before subtracting :math:`d_1`.

An infinite :math:`c` gives :math:`X=0`, and a target with no finite hit has
:math:`S=-\infty`. The nonnegative competitor penalty and capped midpoint
gap give the exact all-target search bound
:math:`S\le -d_1+g_\star\max(b,0)`. A target can be skipped once its bound
is strictly below the best score already found. Equality still needs to
be checked for an exact tie. The settings below were fixed by the matched
tuning study.

.. list-table::
   :header-rows: 1

   * - Receptor
     - :math:`k`
     - :math:`\gamma`
     - :math:`m`
     - :math:`b`
     - :math:`a`
   * - TRA
     - 8
     - :math:`1/10`
     - 12
     - :math:`-1/10`
     - 0
   * - TRB
     - 1
     - :math:`1/4`
     - :math:`17/4`
     - :math:`1/160`
     - :math:`9/32`
   * - Paired-chain
     - 1
     - :math:`1/40`
     - 6
     - :math:`-1/10`
     - 0

Single-chain scores
-------------------

Set :math:`g=\min(24,[(d_1+d_2)/2-d_1]_+)` for the two ordered hits
from one chain. Use the target-hit distances for the appropriate locus
and the native competitor distance :math:`c_k` at order :math:`k`.
Then

.. math::

   \begin{aligned}
   S_{\mathrm{TRA}} &= -d_1-\tfrac1{10}[d_1+12-c_8]_+-\tfrac1{10}g,\\
   S_{\mathrm{TRB}} &= -d_1-\tfrac14[d_1+\tfrac{17}{4}-c_1]_++\tfrac1{160}g.
   \end{aligned}

The TRB support discount is already included in the ordered hit distances.
For example, if a TRB target has only one hit, the repeated distance gives
zero gap contribution. With an empty competitor panel, its score is
then simply :math:`-d_1`.

Assignment and complete-query calibration
-----------------------------------------

Fixed-catalogue motif assignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a TRB junction, trim the same three N-terminal and two C-terminal
positions used by vX. Read up to four overlapping windows of
three consecutive amino acids from each end. The optional replay uses
four consecutive amino acids per window. These sequence windows are
called 3-mers and 4-mers, respectively. For example, the retained
sequence ``RSAGT`` has three-residue windows ``RSA``, ``SAG`` and
``AGT``, and four-residue windows ``RSAG`` and ``SAGT``.
Each token :math:`t` records the window sequence, the end from which
it was counted and its offset from that end. The two lengths test
how much adjacent sequence helps distinguish a target from the
pooled reference.
For :math:`N_e` distinct reference keys for target :math:`e`, let
:math:`c_{et}` count keys with token :math:`t`. Across the full model let
:math:`N=\sum_eN_e` and :math:`c_t=\sum_ec_{et}`. Use

.. math::

   u_t=\frac{1+c_t}{2+N},\qquad
   M_e(q)=\frac{1}{|T(q)|}\sum_{t\in T(q)}
       \log\frac{c_{et}+2u_t}{(N_e+2)u_t}.

The pooled frequency uses one occurrence and one nonoccurrence
pseudocount. The target estimate borrows two reference observations
at that frequency at both motif lengths. The mean is an enrichment score for overlapping
windows. A unique maximum assigns an epitope in the closed-set assay.
Empty query token lists and absent targets are unsupported.

``vdjmatch.evalue.competitive.anchored_motif_scores`` implements the
fixed-catalogue rule. On the 1,218 matched M2 cases it assigns 1,116
correctly with zero ties, including 68 of 87 LLW cases. The full
motif catalogue includes the assay junctions, and the parameters
were selected on this panel. The score is available for fixed-catalogue
assignment diagnostics. It has no calibrated open-set cutoff.
The four-residue replay improves all three class counts to 870 ELA,
232 LLL and 86 LLW, with zero ties. It corrects 72 additional
cases and preserves all correct three-residue calls. This comparison
uses the same full, overlapping catalogue and the supplied assay.

To replay from a prepared model, run:

.. code-block:: bash

   vdjmatch anchored-motif-score --reference reference.tsv --queries queries.tsv --output scores.tsv

The reference columns are ``target``, ``junction_aa``, ``v_call`` and
``j_call``. The command deduplicates this four-column key. Encode the
full peptide–MHC identity in ``target``. The query columns are
``query_id`` and ``junction_aa``; IDs are unique and junctions include
their anchors. The output retains all query–target scores and flags
unique winners and tied maxima.
Add ``--motif-length 4`` to replay the four-residue setting.

Competitive score diagnostics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the score-only assignment diagnostic, write :math:`x_e=\log(1+N_e)`
and maximize :math:`A_\lambda(q,e)=S(q,e)+\lambda x_e` over the three
declared labels. Here :math:`N_e` is the fixed count of native-eligible
target keys before exact-hit masking. When :math:`x_i\ne x_j`, the two label
scores cross at

.. math::

   \lambda_{ij}(q)=\frac{S(q,j)-S(q,i)}{x_i-x_j}.

This partitions the slope axis into intervals with fixed winners. The
original fit selected the midpoint of an interval improving all three assay
recalls, with boundaries approximately :math:`0.00201` and
:math:`0.00365`. The latest vX audit found that the retained midpoint
attains the maximum of :math:`\min_e[R_e(\lambda)-R_e^{\rm native}]`
across the crossing intervals. Here :math:`R_e` is the fraction of assay
records of class :math:`e` uniquely assigned to that class. The gain is flat across multiple intervals. The midpoint maximizes the distance
to the nearer endpoint of the original interval. Use the simpler
:math:`\lambda=1/400` inside the same assignment plateau. It changes the
3,654 saved scores by at most 0.00337 native distance units and preserves
all ties, within-target ranks and within-query label orderings. The term
adds :math:`(\log 2)/400` native distance units per doubling of
:math:`1+N_e`. Its value comes from the three-label assay replay.

The implementation uses ``log1p``. An exact tie or no finite maximum
abstains. The comparison uses queries with finite scores for all three
labels and applies no calibrated score cutoff.

The within-assay prevalence diagnostic uses a second assignment statistic.
Let :math:`p_K(q,e)` be the archived target-specific kernel tail. With
:math:`O_e` the distinct prepared other-peptide reference identities across
HLA restrictions and :math:`K(q,e)` kernel evidence for :math:`e`, this is
:math:`p_K=(1+\#\{r\in O_e:K(r,e)\ge K(q,e)\})/(1+|O_e|)`.
Use zero kernel evidence for no kernel hit. Let
:math:`\widetilde S_e` be the median raw vX score for target :math:`e` across
the scored assay queries and :math:`\pi_e` the observed assay class fraction.
The diagnostic maximizes

.. math::

   H(q,e)=-\log p_K(q,e)+\frac{S(q,e)-\widetilde S_e}{40}
          +0.95\log\pi_e.

Adding the same size correction to every query for a fixed target shifts its
median by the same amount, so the size slope cancels. The evaluation divides
the centered score by ten native distance units and multiplies it by a
quarter. The :math:`0.95` weight on log class frequency and the quarter
weight on centered evidence were selected for the three-class recall
tradeoff in this assay. This rule uses the fixed complete three-label assay,
whose median and class fractions belong to that cohort. Exact ties abstain.

For the generated-receptor study, search the declared compatible pMHC
catalogue and use :math:`S`, including the search across competing targets.
Let :math:`W(q)` be the unique finite winning score, undefined on ties and
unsupported queries. At a higher-is-better cutoff :math:`t`, define
:math:`C_t(q)=\mathbf{1}\{W(q)\text{ defined and }W(q)\ge t\}`. For null bank
:math:`B_s` from source :math:`s`, the empirical false-positive rate is

.. math::

   \widehat f_s(t)=\frac{1}{|B_s|}\sum_{q\in B_s} C_t(q).

The denominator includes ties, no-hits and unsupported queries. The cutoff
fit uses the most permissive observed tied threshold at which all
source-specific rates meet the cap, with an abstain-all option. The generated
single-chain banks come from rearrangement models and the pair bank uses
independent TRA and TRB draws.

Expected counts and repertoire enrichment
-----------------------------------------

After winning-label selection, the empirical expected number of false calls
per query at a cutoff is :math:`\max_s\widehat f_s(t)` over the supplied null
sources. For :math:`Q` query opportunities at that cutoff it gives
:math:`Q\max_s\widehat f_s(t)` expected false calls under the corresponding
empirical null rate. This quantity has units of calls and uses the whole
decision procedure. For a scored query the optional rank tail takes the
maximum across sources of :math:`(1+j_s)/(1+|B_s|)`, with :math:`j_s` equal
to the number of null final scores at least as strong as the query score.
A query with no score has empirical expected calls zero and rank tail one.
Rank-test validity requires exchangeability under the same final-score
pipeline.
For distance-ball E-values with reference-key exposure, see
:doc:`control-calibration`.

For repertoire inference, apply the same membership rule to observed and
mock data. Write :math:`T_0` for observed matched source richness, matched
read count, or the declared finite-effort statistic, and :math:`T_1,\ldots,T_R`
for the same statistic under :math:`R` matched mock draws. The inclusive upper
tail and fold enrichment are

.. math::

   \widehat P_{\rm sample}=\frac{1+\#\{j:T_j\ge T_0\}}{R+1},
   \qquad
   F=\frac{T_0}{R^{-1}\sum_{j=1}^R T_j}.

The ratio is undefined for a zero mock mean. The accepted yellow-fever
curve-mean decision restricts the tail numerator and denominator to mocks
with at least one matched source row, and leaves the tail undefined for a
library with no matched source rows. A call requires a positive standardized
excess and a BH-adjusted value below 0.10 within the declared target,
donor, day and endpoint family. Observed-control mock draws use the
declared Pgen strata. Apply the sample rule at the declared query depth
and target breadth. A per-receptor threshold defines membership. The mock
distribution calibrates the sample statistic.

See :doc:`vx-structure` for structural degree, rarefaction and reciprocal k-NN.
