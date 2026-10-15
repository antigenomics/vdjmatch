Score definitions
=================

See :doc:`selection-information` for selection factors, likelihood ratios,
information in bits, masked-four-mer weighting, geometry and comparable tails.
Benchmark exact exclusion filters hit contributions only: model populations,
priors, motifs and graphs stay fixed. Historical population-puncturing recipes
below do not define the current selection procedure.

Current VDJmatch kernel
-----------------------

The default ``annotate`` model uses the stem4 geometry and the ``v-bounded``
kernel reducer. All quantities below are computed after the declared exact-match
policy. Let :math:`N` be the number of eligible target reference keys,
:math:`t_j` the ordered distance cutoffs, :math:`n_j` and :math:`v_j` the target
and same-V target counts in shell :math:`j`, and :math:`m_j` the cumulative
other-peptide background count through that cutoff. The score is

.. math::

   S(q,X)=\frac{1}{N}\sum_j
   \frac{[0.25n_j+0.75v_j]\exp(-t_j/\tau)}{1+m_j}.

Here :math:`q` is the query, :math:`X` the target peptide–MHC and
:math:`\tau` the kernel scale in the same distance units. Same-V references
receive total weight one; other V genes receive weight 0.25. All target keys
remain in :math:`N`, including distant ones. An unavailable target or background
population produces unavailable evidence. Higher values indicate stronger
relative target support; the score is not a P-value.

The ``all-bounded`` reducer weights every reference equally; ``vonly-bounded``
uses only same-V shell mass. The retained ``rare`` variants replace
:math:`1+m_j` with :math:`\max(m_j,0.01)`. These are explicit alternatives,
not a per-query selection of whichever score performs best.

See :doc:`algorithms` for the computation and paired axes, and
:doc:`explanation` for complete-score calibration. The remaining sections
document retained historical and experimental statistics; they do not redefine
the annotation default.


Experimental count enrichment
-----------------------------

For a fixed distance radius :math:`D`, let :math:`n` be the number of
eligible nonexact epitope-reference receptors in the ball, :math:`N` the
reference population size, :math:`m` the corresponding control count and
:math:`M` its population size. Both population sizes must be positive;
an empty population is unavailable rather than zero evidence.
A direct enrichment statistic is

.. math::

   E=\frac{n/N}{m/M},\qquad
   E_J=\frac{n}{N\widetilde p_B},\qquad
   \widetilde p_B=\frac{m+1/2}{M+1}.

The last expression uses the Jeffreys posterior-predictive probability
for a binary background-ball event. It keeps no-neighbour evidence
(:math:`n=0,N>0`) at zero
and gives a finite denominator without the historical 0.01 floor.
It is an opt-in regularised ranking statistic, not a Bayes factor,
P-value or accepted replacement. Counts, eligible populations, exclusion
rules and distance must agree between reference and control.

Linear distance-weighted evidence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``--neighbour-weighting linear`` retains near-versus-far information inside
the same positive radius, without an exponent or a second scale. Define

.. math::

   K_D(t)=\max(0,1-t/D),\qquad
   W_X=\sum_{r\in R_X}K_D(d(q,r)),\qquad
   W_B=\sum_{b\in B}K_D(d(q,b)),

   \boxed{S_X(q)=\frac{W_X/N}{(W_B+1/2)/(M+1)}}.

Here :math:`R_X,B` are the eligible reference/control populations, :math:`q`
is the query and :math:`d` the identical distance on both sides. All eligible
receptors remain in :math:`N,M`, including distant zero-weight ones. Exact
exclusion precedes both weights and denominators. A boundary match contributes
to the ball count but has zero weight; no reference weight gives zero score.

The smoothed denominator is the posterior mean kernel weight under a
unit-concentration Dirichlet-process prior whose base distribution puts half
its mass at distance zero and half outside the radius. Its prior kernel mean
is therefore one half. This is an explicit regularizer, not a calibrated
specificity probability, Bayes factor or P-value. The binary kernel reduces
to the preceding ball formula. Also :math:`W_B=D^{-1}\int_0^D m(t)\,dt`, where
:math:`m(t)` is the cumulative background count: the taper averages nested
balls rather than choosing their best radius.

For nonzero reference support :math:`n`, the exact change from the ball score is

.. math::

   \frac{S_X}{E_J}
   =\frac{W_X/n}{(W_B+1/2)/(m+1/2)}.

The numerator is the mean reference weight inside the ball. The denominator
is the regularized mean control weight inside the ball. Thus a broad ball
loses score when its reference matches lie farther away than its controls;
near reference evidence gains relative weight. This identity separates the
distance effect from population-size normalization. When :math:`n=0`, both
scores are zero and their ratio is undefined; report zero support directly.

Native linear mass sums integer :math:`\max(0,200D-d_{\rm native})` in the
existing batched count traversal. Division by :math:`200D` precedes the
half-unit regularizer; raw ball counts are exported alongside the weights.
Linked-pair comparison uses :math:`d_\alpha+d_\beta` against the same linked
reference pair and other-peptide linked pairs as its joint population. It does
not multiply independently selected chain matches. Real/generated joint controls
and paired confidence calibration remain undefined.

An additive clipped BLOSUM62 geometry-like geometry uses clipped BLOSUM costs, explicit
positional weights, a linear gap charge and aligned V loops:

.. math::

   c(a,b)=\begin{cases}0,&a=b,\\
        \max(0,\min(4,4-s(a,b))),&a\ne b,\end{cases}

   d(q,r)=d_V(q,r)+\min_A\left[
       3\sum_{(i,j)\in A}w_{L,\max(i,j)}c(q_i,r_j)+12h\right].

Here :math:`s` is BLOSUM62, :math:`q,r` are full junctions,
:math:`d_V` sums aligned CDR1/CDR2/CDR2.5 costs, :math:`A` is an allowed
single-block alignment, :math:`h` is the length difference, :math:`L`
the longer length and :math:`w` the positional weight. clipped BLOSUM62 geometry uses unit
interior weights with three N-terminal and two C-terminal residues excluded.
The CLI's opt-in ``neighbours --junction-ends trim3`` excludes
three residues at each end. Exact exclusion still compares full junctions.
The existing empirical flank/core profile supplies a separate weighting
alternative; it is not an exact gene-specific NDN boundary annotation.

Generative controls measure recombination accessibility, real controls
measure observed repertoire composition, and VDJdb excluding the target
measures separation from other annotated specificities. These are different
background questions. Healthy control metadata does not establish naive-cell
provenance. Sequence Pgen is not the probability of an entire distance ball.
The historical control loader discards V/J. The experimental count CLI instead
retains supplied V calls and counts the same total junction-plus-V-loop distance
on both populations. References and real controls use distinct full-junction/
resolved-V-allele keys. Generated controls retain repeated eligible draws to
estimate generative-law mass; deduplicating them would change that law.
Ambiguous or unmodelled calls are excluded and reported. Present species/locus
metadata are filtered before gene resolution. With exact exclusion, all keys
sharing the query full junction are removed from both counts and denominators.

``neighbours --targets-from-sample --background`` writes these counts
and :math:`E_J` to a separate enrichment table. ``--position-weighting significance``
uses the shipped full-length empirical decay profile, quantised to 0.01 with a
minimum weight of 0.01; no additional terminal trim is applied. The default
nearest-distance comparator remains unweighted. Neither arm changes annotation
confidence, and scientific nonregression remains under evaluation.

V evidence already included in the distance must not be multiplied by a
second unconditional V prior. For a separately declared categorical
stratum :math:`z`, with population sizes :math:`N_z,M_z` and conditional
ball counts :math:`n_z,m_z`, a coherent unsmoothed factorisation is

.. math::

   \frac{N_z/N}{M_z/M}\frac{n_z/N_z}{m_z/M_z}
   =\frac{n_z/N}{m_z/M}.

Read the :download:`revised derivation source <../appendix/score-derivation.tex>`
for gap schemes, historical correspondence, symbols and limiting cases.
The document editor compiles its current preview; the older repository PDF
has not been re-exported for this revision.

Implemented historical and matched-kernel experiments
-----------------------------------------------------

The experimental ``matched-pssm`` score uses one positional distance for both
reference neighbours and control counts. It is a background-normalised ranking
statistic. It has not replaced the default annotation scorer or acquired a
calibrated P-value. These equations describe the implemented single-chain score;
paired joint scoring remains a separate contract.

Populations and geometry
------------------------

Let :math:`q` be an anchor-inclusive query junction, :math:`R_X` the compatible
reference for epitope :math:`X`, and :math:`B` the control population.
Reference representatives are unique full junctions, retaining the first
source-order V call; controls are unique full junctions. Write
:math:`N_X=|R_X|` and :math:`M=|B|`. Exact full-junction identity is excluded
on both sides, without subtracting it from these population denominators.
Unusable reference representatives remain in :math:`N_X` and are reported
separately. Controls are generation/repertoire controls, not the remaining
epitopes in VDJdb; the source manifest identifies the actual population.
The score does not imply that an unverified control is naive or independent.

For BLOSUM62 entries :math:`s(a,b)`, seqtree uses the nonnegative penalty
:math:`\Delta(a,b)=s(a,a)+s(b,b)-2s(a,b)` in its native integer units.
The positional weights are :math:`w_{L,j}=\max(1,\mathrm{round}(100\omega_{L,j}))`,
where :math:`\omega` is the shipped end-anchored significance profile,
normalised to mean one. It is empirical and is not a fitted parameter in
these comparisons. For aligned residue columns :math:`A(q,r)`:

.. math::

   d(q,r)=\sum_{(i,j)\in A(q,r)}w_{L,\max(i,j)}\Delta(q_i,r_j)
          +g(|\ell_q-\ell_r|),\qquad L=\max(\ell_q,\ell_r),

   g(h)=\begin{cases}0,&h=0,\\2800+1400(h-1),&h>0.\end{cases}

There is one length-difference block after six matched prefix residues,
clamped to the shorter junction length; Cys has index zero. Equal-length
junctions are aligned positionwise with no gap. This is a fixed-gap geometry,
not Smith--Waterman or a general alignment. Retrieval uses :math:`D=7000`;
therefore length differences above four cannot enter this score.

Define the historical V factor :math:`a(q,r)` as one for identical
allele-stripped V genes and :math:`0.25\,\mathrm{vsim}(V_q,V_r)` otherwise.
It is separate from the distance and is not a calibrated V-loop likelihood.
This experimental mode does not add clipped BLOSUM62 geometry CDR1/CDR2/CDR2.5 distances.

One kernel and one background predicate
---------------------------------------

The empirical nonexact background mass and kernel are

.. math::

   \widehat F_q(t)=\frac1M\sum_{b\in B}
      \mathbf1\{b\ne q,\ d(q,b)\le t\},\qquad
   K_\tau(t)=e^{-t/\tau}\mathbf1\{t\le D\}.

The implemented density is

.. math::

   \boxed{S_X(q;\tau)=\sum_{r\in R_X,\ r\ne q}
     \frac{a(q,r)K_\tau(d(q,r))}
          {\max\{N_X\widehat F_q(d(q,r)),\varepsilon\}}},
   \qquad\varepsilon=0.01,\quad\tau=400\ \text{by default}.

Equivalently, if :math:`n_X^a(t)=\sum_{r\ne q}a(q,r)\mathbf1\{d(q,r)\le t\}`,

.. math::

   S_X=\int_{[0,D]}\frac{e^{-t/\tau}}
                 {\max\{N_X\widehat F_q(t),\varepsilon\}}\,dn_X^a(t).

This is an integral over observed reference neighbours, weighting close
evidence by the expected unweighted reference count under the same ball
predicate. It combines discrete and gapped edges without a length threshold
or a top-K cap. In a single-distance shell with unit V weights, it reduces
to :math:`e^{-t/\tau}n_X(t)/\max\{N_X\widehat F_q(t),\varepsilon\}`.

Statistical derivation and limits
---------------------------------

For an independent reference drawn from the control distribution, a fixed
ball count has mean :math:`N_XF_q(t)` and, under independent draws, binomial
variance :math:`N_XF_q(t)(1-F_q(t))`. A Poisson approximation is possible
for rare balls. This motivates the expected-count denominator; it does
not make the sum a likelihood ratio or a tail probability. V weights are
applied only to the numerator, so the denominator is not the expected
V-weighted count. The background is empirical and its finite-sample
uncertainty is not integrated into this ranking statistic.

Even with known :math:`F_q` and unit V weights, the null mean is generally

.. math::

   \mathbb E_0[S_X]=N_X\int_{[0,D]}
     \frac{e^{-t/\tau}}{\max\{N_XF_q(t),\varepsilon\}}\,dF_q(t),

not one. Without the floor, its continuous part is
:math:`\int e^{-t/\tau}\,d\log F_q(t)` and can diverge near zero mass.
The floor bounds an individual contribution by :math:`a(q,r)/\varepsilon`;
it preserves a historical regularisation convention rather than an exact
finite-control calibration. A selected-distance score needs its own null
before it can support a P-value.

* No accepted neighbours gives zero density.
* With :math:`\tau\to0`, positive-penalty edges vanish; any nonexact zero-penalty
  edges, if present, retain their contribution.
* With :math:`\tau\to\infty`, exponential attenuation disappears within the
  unchanged cutoff. The result is a sum of inverse-background-mass contributions,
  not a single ball-count ratio.
* Uniform positional weights and forbidden gaps give weighted-BLOSUM Hamming
  geometry; unit substitution penalties give ordinary Hamming geometry.
* Increasing control size with a stable empirical distribution stabilises
  :math:`\widehat F_q`; it does not multiply the score by control size.
* Integer penalties, ties, changing neighbour sets and the cutoff make the
  score discrete. Smooth exponential weighting does not prove global smoothness.

Historical correspondence and output score
------------------------------------------

The original density uses equal-length neighbours with at most five
substitutions. Its positional penalty :math:`p(q,r)` enters the numerator,
but an ordinary unit-edit control ball at substitution count :math:`h(q,r)`
enters the denominator:

.. math::

   S_X^{\rm old}(q)=\sum_{\substack{r\ne q,\ \ell_r=\ell_q\\h(q,r)\le5}}
      \frac{a(q,r)e^{-p(q,r)/400}}
           {\max\{N_X\widehat F_q^{\rm edit}(h(q,r)),0.01\}}.

Matched mode changes the retrieval predicate and replaces this denominator
with the CDF of the same positional distance as the numerator. It rescores
accepted original edges; it does not add a second score to them. Numerical
equality with the old formula requires identical accepted edges, penalties,
V factors and denominator masses. An additive ``historical-pssm`` experiment
instead retains the old density and adds only edges outside its original
substitution ball; these are distinct estimators.

When the original germline prior is requested, the output preserves its
historical fusion. With :math:`G_X(q)` the original smoothed V/J[/length]
log-frequency ratio and :math:`\eta=10^{-6}`:

.. math::

   T_X(q)=\begin{cases}\log(S_X(q)+\eta)+G_X(q),&N_X<500,\\
                         S_X(q),&N_X\ge500.\end{cases}

Without that option, the output score is :math:`S_X`. The 500-reference
switch and prior smoothing are historical reproduction conventions,
not statistical consequences of the kernel derivation. Original paired
cohort rank fusion is not a stable query-level joint score. The existing
``p_enrichment`` remains the original radius-one component test and does
not calibrate :math:`S_X` or :math:`T_X`.

For prior features :math:`f` (V/J/length for beta, V/J for alpha), category
counts :math:`c_{X,f}` and :math:`c_{B_G,f}`, and raw prior population size
:math:`M_G`, the exact historical prior is

.. math::

   G_X(q)=\sum_f\log\frac{(c_{X,f}(f(q))+0.5)/(N_X+30)}
                            {(c_{B_G,f}(f(q))+0.5)/(M_G+30)}.

The prior background :math:`B_G` retains raw multiplicities and is distinct
from :math:`B`. Original alpha reproduction uses its declared original beta
prior background; this historical choice remains explicit.

If all nonexact penalties are positive, :math:`\tau\to0` gives
:math:`T_X\to\log\eta+G_X` in the sparse-reference branch. Changing
:math:`\tau` therefore changes the relative spread of density and prior;
improved density ranking need not improve the fused output ranking.
The historical 500-reference switch is also discontinuous in reference size.

Use ``--gap-geometry matched-pssm --gapped-extension --pssm-kernel-scale 400``.
Scales 200, 400 and 800 provide explicit half/default/double weighting
comparisons while holding geometry, retrieval, controls and original
comparison fixed. These choices are diagnostic, not an accuracy claim.

Distance-local rank evidence (experimental)
-------------------------------------------

``neighbours --neighbour-weighting local-rank --background ...`` uses
fixed integer distance cutoffs from zero through ``--radius``. Both populations
use the same upward distance bins, V-loop geometry and full-junction exclusion.
The score is separate from the existing ``enrichment`` count ratio.

Let :math:`N` and :math:`M` be eligible reference and control sizes, and
:math:`R_i` the rank of reference :math:`i` among the controls plus itself.
With :math:`H_k=\sum_{j=1}^k1/j`, :math:`H_0=0`, define

.. math::

   B_X(q)=\frac1N\sum_{i=1}^N\frac{M+1}{H_{M+1}R_i}.

This restores distance-local inverse background-mass evidence without a fitted
exponent, kernel scale, length switch or target-cluster requirement. Tied ranks
are averaged. If :math:`l` controls are in lower bins and :math:`t` in the same
bin, average over ranks :math:`l+1,\ldots,l+t+1`. Far references are censored,
not discarded: if :math:`m_D` controls fall inside the radius, their contribution
is the mean over ranks :math:`m_D+1,\ldots,M+1`. Full populations remain in
:math:`N,M`. A singleton target is valid; no near evidence on either side gives
neutral :math:`B_X=1`. Empty populations have unavailable evidence and retain
the declared prior, rather than implying a negative prediction.

Under exchangeable distance labels, reference ranks are uniform and the
normalization gives :math:`\mathbb E_0 B_X=1`, even with shared controls.
The alternative label law :math:`dP_1=B_XdP_0` is a specified rank-tilted model.
The CLI exports ``rank_bayes_factor`` for this model, the conservative
``p_rank_bound`` :math:`\min(1,1/B_X)`, and ``posterior_rank_signal``
:math:`\pi B_X/(1-\pi+\pi B_X)`. ``--rank-signal-prior`` declares :math:`\pi`;
equal odds (0.5) is an explicit model default, not a specificity prevalence estimate.
The P-value bound follows from Markov's inequality; it is not an exact tail.

A separate ``p_nearest_rank`` tests first-reference proximity. If :math:`m_*`
controls fall at or below the nearest reference bin, its finite-label tail is
:math:`\Pr\{K\ge1\}`, where :math:`K` is hypergeometric with population
:math:`N+M`, :math:`N` reference labels and draw size :math:`m_*+1`.
Upper-bin ties make this conservative; censored no-hit queries receive one.
It is not the P-value of the full averaged score and is not merged into it.

These are conditional model quantities, not empirically calibrated biological
recognition probabilities or false-positive risks. Generated draws versus
unique reference keys, repertoire selection and database ascertainment do not
automatically support exchangeability. ``rank_model_status`` and the manifest
retain this restriction; deployment priors cannot be inferred from VDJdb or
custom benchmark fractions. Target/route selection is not calibrated here.
This retained experiment does not replace the current kernel default.

The general expectation-one and conservative conversion are discussed by
`Vovk and Wang <https://arxiv.org/abs/1912.06116>`_. The particular reciprocal-rank
alternative and censoring construction above are this package's explicit model.

Uncensored nearest-reference rank (experimental)
------------------------------------------------

``neighbours --nearest-rank --unbounded-nearest --background ...``
adds a separate ``nearest_rank_score`` using uniform junction/V geometry.
For each query, controls are counted at its exhaustive nearest target distance,
including all distance ties. Linked receptors use the alpha+beta distance to
one linked reference pair. This recovers distance resolution when both fixed
radius balls are empty, without a radius scan or a fitted combination.

With eligible target/control populations ``N``/``M`` and ``m_nearest`` controls
at or below that distance, ``p_nearest_unbounded`` is the hypergeometric
first-target tail with population ``N+M``, ``N`` target labels and draw size
``m_nearest+1``. The ranking score is its negative; upper ties are conservative.
Empty populations remain unavailable. Native batched counts use the same V
costs, junction trimming, population keys and optional exact puncture as the
nearest search. Fixed-radius outputs and the averaged rank score are preserved.

This is a conditional exchangeable-label statistic, not biological confidence,
an E-value, a selected-label calibration or a replacement deployment cutoff.
Paired generated/real backgrounds remain undefined. This retained experiment does not replace the current kernel default; promotion
requires all-task and full-reference nonregression.

Nearest-conditioned censored rank (experimental)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
With ``--nearest-rank --neighbour-weighting local-rank``, the CLI also exports
``conditional_rank_bayes_factor``. Existing multi-reference evidence is retained
exactly when any target reference lies within the fixed radius. If all target
references are outside it, the exhaustive nearest reference refines their
null-conditional mean: the first target has reciprocal control rank ``1/(k+1)``,
and the remaining ``N-1`` targets contribute the mean over ranks ``k+1..M+1``.
Here ``N`` and ``M`` are the fixed target and control populations, and ``k``
is the number of controls at or nearer than the exhaustive nearest target.
The usual harmonic normalization is retained. Control ties are conservative.
This common rule preserves linked pairs and uses no fitted blend or task switch.
It supplies conditional model evidence, without a biological probability,
selected-label E-value or transferred probability bound.

Reference-only Kneedle radius (experimental)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The optional radius uses the sorted nonexact nearest-reference distances
:math:`d_i`, with query identity defined by the unchanged full junction
(both junctions together for a linked pair). Define
:math:`x_i=i/(n-1)` and :math:`y_i=(d_i-d_0)/(d_{n-1}-d_0)` for zero-based
indices. For this increasing convex curve, Kneedle considers local maxima of
:math:`x_i-y_i` and confirms a candidate when the difference drops by more
than :math:`1/(n-1)` before the following local minimum. Sensitivity is fixed
at one; the observed monotone curve is used without smoothing. Endpoints,
flat curves and curves shorter than three observations supply no knee.
Among confirmed candidates, the largest normalized difference wins; a minor
leading bump cannot override a stronger later bend. Ties select the earlier bend.
See `Kneedle, section III <https://raghavan.usc.edu/papers/kneedle-simplex11.pdf>`_.

A zero or near-zero knee is floored at half the median positive knee within
one species/chain geometry. With no positive knees, positive nearest distances
provide that scale. Raw and adjusted knees are retained. Multiplying all
distances by a positive constant multiplies both the raw radius and its floor
by that constant; changing substitution weights requires recomputing the curve.
Single-chain and linked-sum geometries are never pooled for this floor.
Reference-selected radii do not inherit fixed-radius probability guarantees.

Conditional J usage diagnostic
------------------------------

``neighbours --conditional-j`` requires ``--nearest-rank``,
``--neighbour-weighting local-rank`` and ``--background vdjdb-other``.
It exports ``j_conditioned_rank_score``: the existing nearest-conditioned
rank evidence multiplied by a target/control J-usage ratio conditional on
resolved V gene and full junction length. Linked pairs use joint V/length
and J tuples from the same reference pair. There is no additional V or
length marginal multiplier and no change to distance geometry.

For each population, the predictive conditional J frequency is
``(n_stratum_J + (n_J + .5)/(N + .5*K))/(n_stratum + 1)``, with one pooled
pseudo-observation and a common number ``K`` of observed population J categories plus one unseen category.
This makes predictions independent of the query batch contents.
Counts retain distinct full junction/resolved-V keys, explicitly mark conflicting
or missing J annotations, and apply benchmark exact-query exclusion consistently.
Missing query J yields a neutral ratio of one. Counts and status are exported.

This is an experimental empirical ranking residual. J usage can remain dependent
on the retained junction sequence; the product is not a calibrated Bayes factor,
P-value, E-value, posterior or biological confidence. Existing rank-model
probabilities describe their original statistic and do not transfer to this score.

Reference-defined epitope rules (experimental)
----------------------------------------------

``neighbours --epitope-rules --background vdjdb-other`` emits
``rule_v_score``, ``rule_rank_score`` and four reference-enriched motif scores.
Spatial and overlap kernels remain diagnostic output fields.
The V rule is estimated separately for each assayed target from distinct full
junction/resolved-V-allele keys. Linked pairs use a joint alpha/beta V category.
There is no benchmark-label fitting, hard-coded epitope exception, J likelihood
multiplier or length likelihood multiplier.

The rank score factors predictive V usage and multi-reference CDR3 rank conditional
on that exact gene category. The native conditional search has **zero V-loop cost**:
V evidence enters once. Length differences remain in the CDR3 distance. With
query-punctured target/control totals N/M, same-V counts n/m and kernel masses
r/b, define p=(n+m+0.5)/(N+M+0.5K), where K is the union of population V categories
plus one unseen category. V probabilities are (n+p)/(N+1) and (m+p)/(M+1).
If either n or m is zero, the V ratio is **one**, CDR3 rank uses the full target
and control populations, and motif bonuses are zero. This avoids prior-only
V evidence or erasing sequence evidence for unsupported categories.
Rank averages normalized reciprocal background ranks over **all references**,
with conservative tied-distance intervals as in the local-rank model above.
Exact distance bins cover every finite CDR3-only distance: clipped substitution
and gap costs are at most twelve per full-junction residue per chain, and distances
are multiples of three. With exact-hit filtering, final hit counts omit those hits; model population sizes stay fixed.

For the retained kernel diagnostics, the shared conditional pseudo-observation is c=(r+b+0.5)/(n+m+1); conditional
kernel means are (r+c)/(n+1) and (b+c)/(m+1). The score is the ratio of their
products. These are empirical predictive ranking scores, not calibrated
P-values, E-values, specificity posteriors or independence claims.

Sequence kernels use full junction coordinates with three residues trimmed at
both ends. The spatial kernel is max(0,1-CDR3_distance/radius). Motif kernels use
exact 3-mers, 3-mers in the eleven `IMGT physicochemical classes
<https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html>`_,
internal-wildcard 4-mers (AXAA/AAXA), and reduced-class wildcard 4-mers.
For each chain, the kernel is the fraction of its distinct query tokens present
in a reference receptor. For linked pairs these fractions are multiplied **within
the same reference pair**, then summed over distinct reference keys. Repeated
or overlapping motifs are averaged, never multiplied as independent likelihoods.
Each mixture averages the spatial and motif kernels with fixed equal weights
before conditional prediction; weights are not tuned on tasks. When a query is
too short for that motif, its mixed score retains the spatial kernel.

For ``rule_enriched_<mode>_score``, each query token looks up its target and control
presence frequencies within the supported V category. Let e be the mean of
max(0, target_frequency-control_frequency) over distinct query token tuples.
The score is ``rule_rank_score*(1+e)``. The bounded factor is a correlated
sequence bonus, not an independent likelihood or calibrated probability.
Rules use reference populations and query puncturing only; query labels never
select motifs. Too-short queries have e=0.

All motif and population counts puncture identical full junctions only when
``--exclude-exact`` is explicitly used for benchmarks. Ordinary annotation retains
exact reference evidence. This option leaves existing comparator/rank fields
unchanged and does not supply a new end-user call threshold.

Experimental original-based fixed refinement
--------------------------------------------

``combine-evidence --fixed`` reports :math:`S=L+s\log B`, where :math:`L`
is the full original log score and :math:`B` is conditional-rank evidence.
The shared sparse weight is :math:`s=1-uc[1-(1-k)e]`: :math:`u` is original
chain support, :math:`c` reference neighbour coverage, :math:`k` normalized
knee shape, and :math:`e` the knee radius divided by that radius plus the
matrix-relative reference radius. These reference features and missing
dispositions are exported separately. One missing score uses the available
component; two missing scores remain unsupported.

The original term retains unit weight. The additional term is an empirical
refinement of correlated evidence, with no independence or biological
probability claim. A fresh raw-task comparison and matching selected-score
null calibration are required before adoption.

``combine-evidence --null-evidence NULL_TABLE --structure-conditioned`` instead
places each component on its matching empirical null-tail scale and reports
:math:`S=\max(-\log p_L,-\log p_R+\log s)`, with the reference-only weight
:math:`s=1-c[1-(1-k)e]`. Query-local support does not switch routes; its
information is already present in the component scores.
Here :math:`p_L` and :math:`p_R` include conservative source tails, ties and
no-hit null queries. At :math:`s=0` only original evidence contributes;
a missing original score uses the available rank tail. Neither component
available remains unsupported. This common unfitted selection avoids raw
score-range domination while preserving the original normalized term. It is
an experimental ranking statistic, not a combined P-value or a validated
accuracy gain. Its complete selection requires separate null calibration.

Experimental ordered class alignment
------------------------------------

``neighbours --alphabet imgt`` uses the existing alignment engine with
the IMGT classes AILV, RHK, CM, ST, DE and NQ; G, F, P, W and Y remain separate.
The cost between two classes is the minimum original clipped BLOSUM62 cost
between their members, with zero within a class. Ambiguous symbols retain their
own classes. Residue order, alignment gaps and the existing positional profile
are preserved. CDR3 and V loops use the same selected matrix; gaps, radius factors
and distance bins derive from that matrix's units.

Raw junctions remain the identity keys. Benchmark exact exclusion removes raw
identical junctions, retaining distinct sequences with zero class distance.
Linked-pair distances still refer to one linked reference receptor. Reference
and background searches use identical geometry. The default ``full`` alphabet
is unchanged. Class distances and ranks are experimental components; their
component P-values do not calibrate a final combined score. Full raw-task
accuracy and selected-score calibration are required before adoption.

Experimental local KF4 chemistry
--------------------------------

``neighbours --chemistry --background vdjdb-other`` emits a separate
local-chemistry component, requiring the existing ``bench`` extra. It preserves
the original distance and rank outputs. Ordered four-residue windows come from
full anchored junctions trimmed by three residues at each end. Every reference
receptor contributes equal total mixture weight, irrespective of window count.
For linked pairs, eight-dimensional windows combine alpha and beta from the
same receptor. Query-window densities are averaged; overlapping windows are
not treated as independent observations.

``--chemistry-bandwidth-factor`` defaults to one and scales the Gaussian
bandwidth by the published VDJtools KF4 table's population standard deviation.
Reference and background use the identical scale. This is a separate chemical
coordinate system; it does not change BLOSUM substitution costs, gaps, V-loop
costs or distance radii.

The output keeps conditional chemical density ratios and V usage ratios
separate. ``chemistry_log_score`` adds the once-only V log ratio to the
conditional log-density ratio. Benchmark exact exclusion removes the entire
full junction (both junctions for a pair) across V categories before every
population denominator. No-window/category and numerical-tail dispositions
remain explicit. Reference keys without local windows are counted separately.
These are experimental ranking values, with no component P/E-value or binding
probability. Complete raw-task accuracy and selected-score calibration remain
required before adoption.
