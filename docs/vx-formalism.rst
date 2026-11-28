Frozen vX research score
========================

This research candidate uses a fixed receptor catalogue. Let :math:`q` be a
TRA receptor, a TRB receptor, or a paired-chain TRA–TRB receptor and let :math:`e` be a
candidate peptide–MHC. A chain is represented by its anchor-inclusive amino
acid junction and mapped V allele. Search keys are distinct encodings. A
key can have several source observations. Paired keys represent observed
links and both chain distances are evaluated against the same key.

The vX contribution is the receptor-to-database evidence score defined
here. It combines a weighted receptor distance with a target's first
two eligible hits and an explicit comparison against HLA-compatible
other-peptide references. The same construction applies to single
chains and genuine paired-chain keys, with chain-specific settings.
TCRdist3 supplies the native distance baseline. Fixed-site gaps, movable
single-block gaps and Smith–Waterman are existing alignment choices.
The novel vX contribution is how distance and reference evidence
are combined into a target score, together with the separately
verified neighborhood rules.


For one chain, let :math:`B(a,b)` be the BLOSUM62 substitution score.
The native comparison charges zero for identical residues. For distinct
canonical residues it charges :math:`\delta(a,b)=\min(4,[4-B(a,b)]_+)`,
with :math:`[x]_+=\max(0,x)`. This is a clipped cost derived from BLOSUM62
with a four-unit ceiling, which is also the unweighted cost of one gap
residue. It discards three
junction residues at the N terminus and two at the C terminus from the
substitution sum. For junctions of different lengths, the native alignment
places one contiguous gap block and charges four units per gap residue.
The junction's weight of three makes the maximum cost of one core
substitution and the cost of one gap residue both 12 native distance units.
Let :math:`d` be the full junction length difference and
:math:`\Pi(q,r)` the permitted alignments of their retained residues.
The accepted vX geometry reproduces the TCRdist3 default junction
alignment. At equal lengths there is one ungapped alignment. With
unequal lengths, one block spans the full length difference. If
the shorter full junction has :math:`s\geq 10` residues, the block
can follow :math:`5,6,\ldots,s-5` residues counted from its N terminus,
giving :math:`s-9` candidate sites. For example, a shorter junction of
length 13 gives four sites. For shorter junctions, widen the initial
interval symmetrically until it contains a site, then intersect it
with legal boundaries of the trimmed core. The trim always removes
three N-terminal and two C-terminal residues from the substitution sum.
A 3-N/3-C research variant trims one more C-terminal residue.

TCRdist3 also allows a fixed position for the single gap and permits
custom metrics including Smith–Waterman local alignment. The
native-default comparator uses the movable gap. The vX local and
global alignment experiments are separate alternatives to that default.
See the `TCRdist3 methods chapter <https://pmc.ncbi.nlm.nih.gov/articles/PMC9719034/>`_
and `distance documentation <https://tcrdist3.readthedocs.io/en/latest/tcrdistances.html>`_. For
the aligned index pairs :math:`(i,j)\in\pi`, the junction cost is

.. math::

   D_J(q,r)=4d+\min_{\pi\in\Pi(q,r)}
      \sum_{(i,j)\in\pi}\delta(q_i,r_j).

Let :math:`D_1,D_2,D_{\rm FR3}` be the aligned costs for CDR1, CDR2 and
the short FR3 segment historically indexed as CDR2.5. The frozen
per-chain geometry is

.. math::

   D_w = 3D_J + D_1 + \frac34 D_2 + \frac32 D_{\rm FR3}
             + \frac{D_1+D_2+D_{\rm FR3}}{256}.

The frozen distance scores the selected FR3 segment and the two CDR loops. In
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
Unconstrained local alignment on nonpositive similarities can always
choose the empty alignment at score zero. Use raw similarities for the
local arm, or require residue coverage in a cost-minimizing alignment.

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
Align each contiguous framework segment independently. An isolated
TRA–TRB contact joins positions on two chains and is evaluated as a
pair of sites under the interface term. It does not define a continuous
sequence segment for alignment.

The human FR3 research candidate blends two costs for the same short
germline segment. With :math:`w_{\rm FR3}=3/2+1/256` and
:math:`\eta=1/512`, replace this segment's contribution by

.. math::

   w_{\rm FR3}\bigl[(1-\eta)D_{\rm FR3}+\eta D_{\rm SW}\bigr].

Here :math:`D_{\rm SW}` uses the self-normalized Smith–Waterman
dissimilarity above with raw BLOSUM62 and gap opening and extension
parameters both set to two. The compiled correction uses the equivalent
residual :math:`\eta w_{\rm FR3}(D_{\rm SW}-D_{\rm FR3})` for each chain.
The replacement covers the existing FR3 segment only. CDR1, CDR2,
the 3-N/2-C trimmed one-gap junction score, and any separately
configured paired interface term retain their baseline definitions.
The same native fallback applies when the frozen correction is
unavailable. Full-junction global and local alternatives and the
separate 3-N/4-C flank split caused regressions in the organizer
replay.

The fraction :math:`1/512` was selected by matched benchmark replay.
At :math:`1/256` the local variant lowered two fixed-reference paired
TCRvdb endpoints. The smaller local fraction improved five organizer
and four fixed-reference endpoints with zero losses, keeping all six
measured NLV endpoints and all three LLL–LLW–ELA correct counts.
This is an explicit human research component. Whole-search
calibration and broader verification remain pending.

Expanded gap placement candidate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The vX candidate expands the set of legal gap placements by one
boundary at each end and keeps the gap length fixed at the
difference in junction lengths. For a shorter junction of length 13,
the native boundaries follow residues 5--8 and the expanded boundaries
follow residues 4--9. Both sets are limited to legal boundaries in the
trimmed core. Equal-length junctions retain their ungapped distance.

Let :math:`D_J` be the native-window cost and :math:`D_J^+` the cost
minimized over the expanded window, using the same residue and gap
costs. The proposed CDR3 contribution per chain is

.. math::

   3\,\frac{511D_J+D_J^+}{512}.

The expanded alignment can lower the cost only when a newly eligible
placement scores better. The small fraction was selected on the
exact-masked human paired organizer panel. It improved all three
TTDPSFLGRY endpoints, with the other 21 unchanged. Transfer to the
fixed-reference priority panels is required before using the change
in the accepted scorer.

Substitution matrix choice
~~~~~~~~~~~~~~~~~~~~~~~~~~

With the same target keys, exact-hit mask, loop weights and reducer,
the matched matrix replay replaced only the target distance. For
the unclipped BLOSUM62 arm, the residue cost was
:math:`B(a,a)+B(b,b)-2B(a,b)`. Each alternative distance was scaled
so that its largest substitution cost was four. The native
other-peptide comparison stayed clipped. This is a comparison of
costs in the same single-gap geometry. The Smith–Waterman FR3
comparison above uses a separate alignment objective.

Clipped costs beat native TCRdist3 on all 78 matched human ranking
endpoints, compared with 65 for VDJAM, 63 for unclipped BLOSUM62, 53
for PAM250 and 47 for MJ. In a direct comparison on those 78
endpoints, clipped won 33, tied 2 and lost 43 against unclipped
BLOSUM62. The unclipped arm lowered 30 of the 72 TCRvdb endpoints
relative to frozen vX and lowered measured NLV TRB average precision
by 0.02108. Clipped costs remain the reference-ranked choice under
the priority nonregression requirement. The direct endpoint tally
supports a trade-off between gains and losses across tasks.

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
gap residue. Limit the total difference between the first and second
hit distances to :math:`4LC`, equivalent to four maximum-cost residue
edits per observed chain. The midpoint cap is :math:`g_\star=2LC`.
In real arithmetic,

.. math::

   \begin{aligned}
   X &= [d_1+Lm-c]_+,\\
   G &= \tfrac12\min\{d_2-d_1,4LC\},\\
   S(q,e) &= -d_1-\gamma X+bG.
   \end{aligned}

The leading term :math:`-d_1` rewards the best hit. The per-chain margin
:math:`m` sets its desired lead over the background, and :math:`\gamma`
sets the cost of falling short. The order :math:`k` sets the background
comparison rank. The gap weight :math:`b` penalizes separation for TRA and
paired-chain receptors and rewards it slightly for TRB. All distance terms have native
distance units. The weights :math:`\gamma,b` are dimensionless.

The gap contribution is in :math:`[-2.4,0]` for TRA, :math:`[0,0.15]` for
TRB and :math:`[-4.8,0]` for paired-chain receptors. An equal-length
junction of length :math:`\ell` can contribute up to
:math:`12[\ell-5]_+` before V-loop terms, and other terms or unequal lengths
change the range further. In the matched TCRvdb ranking queries
(2,636 records per axis), median full junction lengths were 13 for TRA
and 14 for TRB. Their median core lengths of eight and nine give loose
substitution-only midpoint scales of 48 and 54 units. A length-based cap
would set :math:`g_\star(q)=6\sum_i[\ell_i-5]_+`. Here :math:`\ell_i`
is the full junction length of observed chain :math:`i`. Neither this
capacity nor the cost matrix determines how far apart the two best
reference hits are in a real panel. The original support-gap sweep fixed the
cap across its variants.
Its value remains a design choice. In the matched ranking replay,
the midpoint gap exceeded the frozen cap in 69 of 8,295 finite
query–target scores across TCRvdb and measured NLV. Halving the cap
changed 374 scores and lowered 32 of 78 endpoints.
Doubling it or using the query-core capacity each changed 69 scores and
lowered six endpoints. All cap options beat native TCRdist3 on the 78
endpoints. The frozen cap also preserved all frozen vX endpoints.
Binary64 evaluation forms the midpoint before subtracting :math:`d_1`,
then limits the result to :math:`[0,g_\star]`.

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

Paired-chain score and neighborhood decisions
------------------------------------------------

For a genuine paired-chain reference, first add the two chain distances.
Let :math:`d_1,d_2` be the first two ordered distances to the same target
and :math:`c_1` the nearest compatible native other-peptide distance.
The paired score applies the competitive terms once to those joint
distances:

.. math::

   S_{\alpha\beta}=-d_1-\frac{[d_1+12-c_1]_+}{40}
     -\frac{\min\{48,[(d_1+d_2)/2-d_1]_+\}}{10}.

The margin 12 is twice the fitted six units per chain. The cap 48
is twice the single-chain cap. Separate best hits on the two chains
cannot supply :math:`d_1`.

On the paired organizer panel, the vX replay evaluates within-target
neighborhood density.
For a group of :math:`n` constructs, mask the query itself and any
additional exact paired-junction hits specified by the protocol. Order
the remaining distances :math:`d_{(1)}\leq\cdots\leq d_{(k)}` and use

.. math::

   S_{\text{assay}}=-\frac{\sum_{j=1}^k (k-j+1)^{7/3}d_{(j)}}
                              {\sum_{j=1}^k(k-j+1)^{7/3}},\qquad
   k=\max\{1,\lfloor 31(n-1)/250\rfloor\}.

The decreasing rank weights emphasize the closest neighbors. The
organizer assay selected the fraction and exponent. All eight observed
neighbor counts stay fixed for fractions in
:math:`[15/121,16/129)`. The endpoints have cross-product difference
one, so their mediant :math:`31/250` has the smallest denominator among
fractions strictly between them. It equals 12.4%. One eighth would add
one neighbor in four groups. If fewer than :math:`k` distances are
finite, the case is unscored. Equal distances
at the boundary have equal effect on this distance average. This score
operates within an already specified target group.

For cross-target classification, a count-based neighborhood can instead
use the same receptor distances with a radius :math:`r`. For the finite
eligible reference set :math:`H(q)` and the labels :math:`L(h)` of reference
key :math:`h`, define

.. math::

   V_e(q;r)=\sum_{h\in H(q)}
     \mathbf{1}\{D(q,h)\leq r\}\mathbf{1}\{e\in L(h)\}.

One fixed-radius rule chooses a unique maximum positive vote, subject to
a declared minimum support. A kNN version sets :math:`r` to the
:math:`k`-th finite distance, using the last distance if fewer than
:math:`k` eligible hits exist. It includes all ties at the boundary.
Count each distinct physical key once per supported label and
abstain for no eligible hit, insufficient support, or tied target votes.
This is a classification rule for evaluation. The accepted organizer
density score and the fixed-reference evidence score above retain
their respective definitions and cutoff scopes.

Epitope enrichment across a reference catalogue
------------------------------------------------

In the fixed-reference vX score, :math:`d_1,d_2` are taken separately
for each candidate pMHC. A nearer key for a different peptide can
influence the background comparison, and another candidate can win.
For a cross-target neighborhood, count each distinct physical key
once. If it supports several pMHC labels, count it once for each
supported label while keeping one key in the neighborhood denominator.

Let :math:`N` be the size of the fixed catalogue, :math:`N_e` the
number of keys supporting target :math:`e`, :math:`m` the number
of finite neighbors including ties, and :math:`x_e` the number of
those neighbors supporting :math:`e`. The local fraction, background
fraction and observed enrichment ratio are

.. math::

   f_e=x_e/m,\qquad \pi_e=N_e/N,\qquad R_e=f_e/\pi_e.

For a neighborhood defined by distance alone, a finite-catalogue
reference model draws :math:`m` keys without replacement under
exchangeable target membership. Rank target support using the upper
hypergeometric tail

.. math::

   p_e = \Pr\{X_e\geq x_e\},\qquad
   X_e\sim\operatorname{Hypergeom}(N,N_e,m).

A binomial with parameters :math:`m,\pi_e` approximates this tail
when the neighborhood is small relative to the catalogue. Exact-hit
exclusion only changes the benchmark hit list. The fixed reference
fractions remain unchanged. A target can be ranked by :math:`-\log p_e`
when :math:`x_e>0`, with abstention at a tied maximum or without support.

This uses the observed-versus-expected neighborhood principle of
`TCRNET <https://pmc.ncbi.nlm.nih.gov/articles/PMC6156623/>`_ for a
different question. TCRNET tests within-sample neighborhood degree
against a control repertoire. TCRdist3 offers neighborhood tests
for covariates, including a documented chi-squared example with
epitope as the category. It can also choose radii that control
neighbor discovery in a background repertoire. The vX design
combines label enrichment from a fixed multi-label catalogue with
its receptor distance and a query-level open-set decision. The
one-sided hypergeometric tail is the Fisher exact enrichment test
for one label inside versus outside the neighborhood. It remains
defined for sparse counts.
See the `TCRdist3 enrichment workflow <https://tcrdist3.readthedocs.io/en/latest/#neighborhood-enrichment>`_
and `radius procedure <https://tcrdist3.readthedocs.io/en/latest/radius.html>`_.

The exchangeability model is a ranking reference. Study composition,
HLA restriction and repeated biological observations can give targets
different exposure to the same kinds of receptors. A stronger call
procedure prespecifies the compatible catalogue and neighborhood
rule. If the search requires a V match, use the same fixed V stratum
in both the neighborhood and the reference fractions. Then calibrate
the maximum target evidence on held-out null
queries after rerunning the entire label search. If the rule searches
several neighborhood sizes, the null replay repeats that search.
Test generated and withheld-specificity nulls separately. The resulting
empirical any-call rate applies to that decision and null population.

Peptide and restriction ambiguity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Keep all labels attaining the exact maximum when auditing an abstention.
If several winning pMHC labels share one peptide, the peptide is resolved
at the score maximum and the MHC restriction remains ambiguous.
If the maxima span distinct peptides, peptide assignment also abstains.

On the full paired TCRvdb catalogue replay, 1,360 of 1,404 tied maxima
among resolved assay cases involved a single peptide. Reporting that
peptide at a tie raises correct peptide calls from 376 to 1,040 among
1,325 assay positives. It also raises calls of the assayed peptide
from 48 to 291 among 1,032 assay negatives before a score cutoff.
This is a peptide-level diagnostic. Applying it as a call rule
requires calibrating its own decision, including the newly admitted ties.

Sample-level neighborhood association
--------------------------------------

Fix the vX distance, neighborhood centers and radii before testing
sample covariates. In sample :math:`i`, let :math:`X_i(q,r)` be the
number of distinct eligible receptor keys within radius :math:`r`
of center :math:`q`, and let :math:`N_i` be the total eligible keys
in the corresponding sample or prespecified stratum. For
:math:`N_i>0`, the sample-level neighborhood fraction is

.. math::

   Y_i(q,r)=X_i(q,r)/N_i.

Handle zero-exposure samples by the prespecified design and retain
their counts in cohort accounting. A paired comparison with
:math:`n` donors uses the mean within-donor change,

.. math::

   T(q,r)=\frac{1}{n}\sum_{i=1}^n
          [Y_{i,1}(q,r)-Y_{i,0}(q,r)].

The two indices denote the paired conditions. One donor contributes
one contrast. For a randomized paired design, swap the condition
labels within each donor only as allowed by the randomization. The
exact null is that assignment has no effect on any donor's outcome.
An observational design requires joint swap invariance under its
conditional null. Use the same swaps for all overlapping
neighborhoods. For a global test of any association, use the maximum
absolute statistic across the tested centers and radii. Let
:math:`M_{\rm obs}` be that maximum with the observed labels and
:math:`M_\sigma` its value under an allowed rearrangement
:math:`\sigma\in\mathcal G`, including the identity. The exact
global tail for equally likely assignments is

.. math::

   p_{\rm global} =
     \frac{\#\{\sigma\in\mathcal G:M_\sigma\geq M_{\rm obs}\}}{|\mathcal G|}.

Repeat any label-dependent selection inside each rearrangement.
This formulation targets sample association. It does not turn
receptor keys from the same donor into independent replicates.

Background-controlled radius selection
--------------------------------------

Let :math:`A_r(x)` indicate that the entire annotation rule calls
background query :math:`x` at radius :math:`r`, including target
selection and abstention. Choose the radius using a training
background and freeze the decision before evaluation. On an
independent evaluation null of :math:`B>0` queries, count
:math:`C=\sum_{j=1}^B A_r(x_j)` calls. The observed any-call rate is
:math:`C/B`. For independent queries with a common false-call
probability, a one-sided exact upper confidence limit at confidence
:math:`1-\delta` is

.. math::

   U=\operatorname{Beta}^{-1}(1-\delta;\,C+1,B-C),\qquad 0\leq C<B.

Here :math:`\operatorname{Beta}^{-1}` denotes the beta quantile.
Use :math:`U=1` if :math:`C=B`. A bound for a fixed rule does
not justify selecting a radius on the same evaluation queries.
For clustered queries, use an interval that respects the sampling
design. The expected number of background calls across :math:`Q`
query opportunities at the same marginal rate is estimated by
:math:`QC/B`. Reference-hit counts inside a ball are a separate
quantity.

These are proposed vX procedures and specify a matched comparison
with the TCRdist3 neighborhood tests and empirical background-radius
workflow. Compare the same samples, truth criteria and background
sources. Measure power or correct-call recall at matched
false-call control before claiming a performance advantage.

Calibration from an ordered null sample
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a single cutoff on a fixed whole-query score, choose the allowed
calibration count before examining the null scores. Let :math:`B` be
the number of independent draws from the declared null, :math:`k<B`
the allowed count, and :math:`Z_{(B-k)}` the ascending order statistic
of the calibration scores. Give abstentions score :math:`-\infty`.
Call only when the new score is strictly above :math:`Z_{(B-k)}`.
This excludes the boundary tie in full. If :math:`f` is the resulting
population any-call probability and :math:`\alpha` its intended cap,

.. math::

   \Pr(f>\alpha)\leq\Pr\{\operatorname{Binomial}(B,\alpha)\leq k\}.

The probability on the left is over calibration samples. The one-sided
upper bound at confidence :math:`1-\delta` is the beta quantile with
shape parameters :math:`k+1` and :math:`B-k`. Exceeding the cap
requires the boundary score to fall below the null
:math:`(1-\alpha)` quantile. In that case at most :math:`k` calibration
scores reach a quantile whose upper tail has mass at least
:math:`\alpha`. The binomial bound is tight for a continuous score
and conservative with ties. See `Tong, Feng and Li, Proposition 1
<https://arxiv.org/abs/1608.03109>`_.

Choose the largest :math:`k` whose binomial tail is at most :math:`\delta`.
For 100,000 nulls, a 1% cap and 95% confidence give :math:`k=948`.
Allowing 1,000 calls instead gives an upper bound of 1.0533%. Use an
abstain-all rule when no count qualifies. Freeze the catalogue, score,
and label decision independently of this sample. A radius family uses
the same argument when its call sets are nested and admit one scalar
score. A family in which adding neighbors can create an abstention
uses an independent evaluation background after radius selection.

Use ``vdjmatch calibrate-ordered-calls`` to fit this boundary and
``vdjmatch apply-call-calibration`` to apply the saved recipe.

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

Paired-chain re-ranking by motif conflict
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each target, let :math:`S_e` be the base score and
:math:`M_{\alpha e}, M_{\beta e}` the fixed reference motif enrichment
on the two chains. Treat missing motif values as zero and let
:math:`M_e=\min(M_{\alpha e},M_{\beta e})`. The selected mouse
rule penalizes conflict only. For a supported target strictly below the
query maximum :math:`S_*` and with more than two distinct reference
keys, the refined score is

.. math::

   S'_e = S_e+4\min(0,M_{\alpha e}, M_{\beta e}).

Every other target keeps its base value. The general reducer caps
changed lower-tier scores at the greatest representable value below
:math:`S_*`. The conflict-only correction can only lower a score,
so that cap is already satisfied by the displayed rule.
The motif score is independent of gap placement and is applied after
pairwise distances have been reduced to target scores. The mouse
ranking rule does not add a k-mer term to the pairwise distance.
The penalty weight four was selected on the mouse catalogue benchmark.
The two-key gate applies the penalty only to targets with more than two
reference keys. On the
matched class-I and class-II mouse paired catalogue it improved all 39
ranking endpoints versus both frozen vX and native TCRdist3.

For the same base scores and full target context, the maximum and the
set of targets attaining it are exactly invariant. Consequently any final
decision that uses only a unique maximum and a fixed cutoff has identical
calls, including on null queries. This holds pointwise, so null call counts
and the FPR curve are identical at every cutoff on a fixed null bank.
The refinement provides no improvement to that detector's FPR. A change
to the base geometry or any downstream decision based on lower ranks
requires its own calibration check.

Replay with ``vdjmatch paired-motif-rank`` using the complete target
axis and ``--support-weight 0 --conflict-weight 4``. The command accepts
either aligned base and chain motif matrices with target key counts,
or an archive with the base plane and a prepared fixed reference:

.. code-block:: bash

   vdjmatch paired-motif-rank --scores base.npz --queries queries.tsv \
     --targets targets.tsv --reference memberships.tsv \
     --support-weight 0 --conflict-weight 4 --output ranks.tsv

In fixed-reference mode, the archive contains ``base`` in the supplied
query and target order. The query TSV contains ``query_id``,
``junction_aa_alpha`` and ``junction_aa_beta``. The target TSV contains
``target`` and ``mhc_class`` (``MHCI`` or ``MHCII``). Each reference row
contains ``reference_key``, ``target``, and both full junctions. Use the
same key for a physical pair across its target memberships and distinct
keys for distinct reference pairs even when their junctions coincide.
The motif fit pools memberships separately by MHC class; a shared key
contributes once per target membership. The command constructs both
motif planes and derives the target key counts. The supplied base
plane must already reflect the chosen hit-list policy.

For the mouse benchmark, produce that plane and the three TSV files
from a pinned VDJdb release ZIP and a paired query table:

.. code-block:: bash

   vdjmatch paired-vx-base paired-queries.tsv --vdjdb reference.zip \
     --exclude-exact --output-dir mouse-base

Pass the four files from ``mouse-base`` to ``paired-motif-rank``.
Wide query tables use ``junction_aa_alpha/beta`` and ``v_alpha/beta``
with a pair ID. The command retains the full compatible catalogue
and scores genuine pairs. Exact-hit exclusion operates on both full
junctions and applies only to the target and background hit lists.
Omit ``--exclude-exact`` for ordinary exact-retained replay. This
command uses the verified mouse frozen base before motif re-ranking.

The output also includes ``maxima.tsv`` with one row for every query.
The ``maximum`` column retains finite tied maxima for audit, whereas
``score`` and ``target`` require a unique winning pMHC. Unsupported and
tied queries stay in the input denominator. For a declared calibration
null, label the input with ``--null-source`` and pass this file directly
to the ordered calibration command:

.. code-block:: bash

   vdjmatch paired-vx-base null-queries.tsv --vdjdb reference.zip \
     --null-source generated --output-dir mouse-null
   vdjmatch calibrate-ordered-calls mouse-null/maxima.tsv \
     --target-fpr 0.01 --confidence 0.95 --output mouse-calibration.json

Use the same reference and exact-hit policy when scoring the queries
to which the recipe will be applied. The mouse motif correction
preserves this maximum-based calling rule.

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
For numerical reproduction of the historical scores within :math:`10^{-5}`,
round the historical coefficient to :math:`0.002832`. The largest observed
score shift then is :math:`3.73\times10^{-6}` native distance units, with
the same ties, ranks and assignments in all 3,654 slots. The smallest
historical within-query label gap is :math:`0.00211` native distance units.

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

Recorded-target precision
-------------------------

Among catalogue positives, let :math:`r_t` be the fraction called
with a recorded target and :math:`u_t` the fraction called with a
different target at cutoff :math:`t`. For a mixture with positive
fraction :math:`\pi` and generated background source :math:`s`, the
fraction of all calls that recover a recorded target is

.. math::

   \operatorname{Precision}_{s,\pi}(t)=
   \frac{\pi r_t}{\pi(r_t+u_t)+(1-\pi)\widehat f_s(t)}.

Precision is undefined when the mixture produces no calls. Count
other-target calls on positives in the denominator. They receive no
credit under the recorded catalogue annotation. The estimates apply
to the stated mixture and its conditional query populations.

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
For a hit at distance :math:`d`, let :math:`R_e` be the fixed set
of eligible reference opportunities for target :math:`e`. With null
receptor :math:`Q_0` drawn under the declared background law, the
expected number of reference hits at least as close is

.. math::

   E_{\rm hits}(d,e)=\sum_{r\in R_e}\Pr[D(Q_0,r)\le d].

This sum uses linearity of expectation and permits dependent hits.
Estimating the probability of at least one hit requires the joint
distribution or a justified approximation. Empirical query-level
FPR comes from replaying the complete pipeline on null queries.
A suitable comparison between methods uses the same null queries,
target catalogue and matched positive recall or a prespecified
operating point. Report the null count and its denominator separately
for each source. For distance-ball E-values with reference-key
exposure, see :doc:`control-calibration`.

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
