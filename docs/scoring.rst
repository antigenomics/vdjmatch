Neighbourhood score: formula and interpretation
===============================================

The experimental ``matched-pssm`` score uses one positional distance for both
reference neighbours and control counts. It is a background-normalised ranking
statistic. It has not replaced the default annotation scorer or acquired a
calibrated P-value. These equations describe the implemented single-chain score;
paired joint scoring remains a separate contract.

Download the :download:`two-page derivation <../appendix/score-derivation.pdf>`
or its :download:`editable LaTeX source <../appendix/score-derivation.tex>`.

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
This experimental mode does not add TCRdist CDR1/CDR2/CDR2.5 distances.

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
