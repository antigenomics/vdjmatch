Rearrangement-aware substitution and indel models
=================================================

**Rearrangement origin changes the meaning of a substitution matrix.**
Two TCR junctions recognizing the same epitope need not descend from one
another. Their differences usually compare independently generated,
convergent V(D)J rearrangements. An aligned identity within a retained
germline segment is therefore not evidence that antigen recognition
selected against a mutation at that position. Under the model considered
here, retained germline is immutable and allelic mutation is excluded.
We distinguish this structural constraint from trimming and palindromic
P additions at segment boundaries, and from stochastic N additions. The
resulting matrix describes compatibility of alternative rearrangements
with recognition, rather than an evolutionary mutation process.

The original BLOSUM construction pools aligned residue pairs, including
identities, across blocks and controls sequence redundancy; PAM
additionally estimates residue-specific mutability and propagates a
transition process. Neither method literally requires equal mutability
of every residue. The problem in applying their pooled counts directly
to junctions is the absence of rearrangement-origin conditioning. The
original sources are `Henikoff and Henikoff
(1992) <https://doi.org/10.1073/pnas.89.22.10915>`__ and `Dayhoff,
Schwartz and Orcutt (1978), chapter
22 <https://compbio.berkeley.edu/class/c246/Reading/dayhoff-1978-apss.pdf>`__.
The construction below is an adaptation of these principles to
rearrangement convergence, not an assertion of receptor ancestry.

**Generative model and origin coordinates.** Let :math:`C` denote chain,
:math:`L` the amino-acid junction length, and :math:`x` an
anchor-inclusive amino-acid junction. Let :math:`r` be a complete
rearrangement history: gene choices, retained segment coordinates,
trimming/P events, N lengths and N nucleotide realizations. Write
:math:`T(r)=x` for translation and :math:`P_0(r\mid C)` for the
specified recombination model. The observed gene calls enter as
conditioning evidence with their ambiguity retained. Known retained
templates are not repaired by allowing mutations inside germline.

For each nucleotide coordinate :math:`t`, retain an origin label
:math:`Z_r(t)\in\{G,P,N\}`. The germline label :math:`G` also records
the gene and reference coordinate. An amino-acid site :math:`p` has a
three-base origin

.. math::


    Z_{r,p}=\bigl(Z_r(3p),Z_r(3p+1),Z_r(3p+2)\bigr).

Only :math:`NNN` denotes a pure-N codon. A codon spanning a segment
boundary is a mixed-origin codon, even when its amino acid differs
between two receptors. P sequence is templated by the palindromic
germline extension, whereas its presence and length vary through
boundary events. N sequence is stochastic under the insertion model;
nucleotide bias, dinucleotide dependence, reading frame and translation
make its amino-acid frequencies nonuniform. These three mechanisms
consequently require different residue and length baselines.

**Marginalize histories instead of assigning one history as truth.** For
positive model likelihood, the conditional origin probability is

.. math::


    \omega_{x,p}(z)=
    \frac{\displaystyle\sum_r P_0(r\mid C,\mathcal V,\mathcal J)
          \mathbf1\{T(r)=x,\ Z_{r,p}=z\}}
         {\displaystyle\sum_r P_0(r\mid C,\mathcal V,\mathcal J)
          \mathbf1\{T(r)=x\}},
    \qquad z\in\{G,P,N\}^{3}.

Here :math:`\mathcal V,\mathcal J` denote the gene-call evidence used by
the model. A zero denominator makes origin unresolved; it does not
establish N origin. Top-ranked histories and a most-likely inferred
nucleotide sequence are not this posterior. In particular,
model-inferred nucleotides cannot independently validate the model that
generated them.

If epitope selection depends only on the observed junction and the
conditioned receptor features, it cancels from the posterior over
histories given those features. For example,

.. math::


    P_e(r)\ \propto\ P_0(r)\,\sigma_e(T(r),\mathcal V,\mathcal J)
    \quad\Longrightarrow\quad
    P_e(r\mid x,\mathcal V,\mathcal J)=P_0(r\mid x,\mathcal V,\mathcal J).

This cancellation requires the stated conditioning. It does not hold
automatically when selection depends on unobserved paired-chain, gene or
structural features. The origin probabilities remain model-based
quantities.

**Structural conservation and recognition conservation.** At a fixed
retained gene/reference coordinate, :math:`P(\mathrm{change}\mid G)=0`
by construction in this model. Such identities must be omitted from both
the substitution numerator and the conservation denominator of a learned
N-region core. Removing germline mismatches while retaining their
identities would spuriously inflate conservation. Conversely, simply
removing conserved terminal anchors is insufficient: templated residues
can extend into the junction interior. A difference between V genes, or
between retained reference coordinates after trimming, is a
gene/scenario difference and belongs in the gene or boundary component.
It is not a substitution within a fixed germline template. Uncertain
origin is retained separately rather than silently included in N.

For independent rearrangements :math:`x,y`, a pure-N pair at aligned
positions :math:`p,q` has weight

.. math::


    u_{xy,pq}=\omega_{x,p}(NNN)\,\omega_{y,q}(NNN),

conditional on the independence and model assumptions above. Multiplying
this by the declared record, clonotype, redundancy and epitope weights
gives an origin-weighted pair contribution. Pair definitions, exact
matches, contradictory labels and source identities are retained in the
audit. Within each declared comparison, ordered counts are used
consistently for identities and exchanges.

**BLOSUM-style convergence counts.** Let :math:`h` index aligned site
pairs in the chosen blocks, and let :math:`w_h` denote their sampling
weights. The pure-N counts for epitope :math:`e` are

.. math::


    F^{N,e}_{ab}=
    \sum_{h=(x,y,p,q)} w_h\,\mathbf1\{e_x=e_y=e\}\,
    u_{xy,pq}\,\mathbf1\{x_p=a,\ y_q=b\}.

Only distinct observations allowed by the declared redundancy convention
contribute. Identity counts :math:`a=b` remain essential, but only on
the same N-origin opportunity as exchanges. The block alignment, cluster
definition and origin weights are separate operations. An aligner
proposes correspondence; it does not establish the rearrangement
history.

Normalize the weighted counts to :math:`q^+_{ab}`. A baseline
:math:`q^0_{ab}` must use the same chain, length, position,
gene/scenario, origin and pair-ascertainment conditions. It may come
from the specified recombination model or from a declared, matched
different-epitope comparison. Different-epitope repertoires are
themselves selected populations and must not be called an unselected
generative background without justification. The recognition score and
the scorer-compatible cost are

.. math::


    S_{ab}=\log_2\frac{q^+_{ab}}{q^0_{ab}},\qquad
    D_{ab}=S_{aa}+S_{bb}-2S_{ab},\qquad
    d_{ab}=\max(0,D_{ab}).

Consequently,

.. math::


    D_{ab}=
    \log_2\frac{q^+_{aa}q^+_{bb}}{(q^+_{ab})^2}
    -
    \log_2\frac{q^0_{aa}q^0_{bb}}{(q^0_{ab})^2}.

A factorized background :math:`q^0_{ab}=p_ap_b` cancels from this
distance. A structured origin-conditioned background need not cancel.
This distinguishes correcting pooled residue frequencies from correcting
rearrangement opportunity. The diagonal is zero after conversion, but
the self-pair counts used to derive it cannot be omitted. Missing
origin-conditioned support is reported; pseudocounts cannot create
evidence for structurally forbidden events. A BLOSUM prior, if used for
shrinkage, is declared as reference supervision and is not evidence of
recovering BLOSUM from VDJdb.

**The bridge from codons to a genetic-code baseline.** Let
:math:`u_s(c)` be the productive pure-N codon distribution in a fully
specified stratum :math:`s`, and let :math:`\tau(c)` translate codon
:math:`c`. For independent codon draws in that stratum,

.. math::


    f_{s,a}=\sum_{c:\tau(c)=a}u_s(c),\qquad
    q^0_{s,ab}=f_{s,a}f_{s,b}.

This factorization is a conditional construction; it does not assert
independence of junction positions. If the ascertainment rule also fixes
shared flanking context :math:`h`, the pooled opportunity is instead a
mixture :math:`\sum_h w_h f_{s,h,a}f_{s,h,b}`, where :math:`w_h` are
normalized context weights under that rule. Such context mixing must be
accounted for before assigning its conservation to recognition.

To obtain a one-nucleotide genetic-code proposal, explicitly define a
codon coupling :math:`J_s(c,d)`, for example one supported on sense
codons separated by one nucleotide, and project it through translation:

.. math::


    K^0_{s,ab}=\sum_{c,d}J_s(c,d)
    \mathbf1\{\tau(c)=a,\ \tau(d)=b\}.

Here :math:`J_s` is a specified joint proposal mass, and :math:`K^0_s`
is its amino-acid pushforward; synonymous proposals contribute to the
diagonal. This is a compatibility hypothesis for independently
rearranged receptors, with no implied mutational ancestry. It differs
from independent draws and must be declared when interpreting a
genetic-core comparison. For a recognition multiplier
:math:`R_{s,ab}\ge0`, write

.. math::


    q^+_{s,ab}=\frac{K^0_{s,ab}R_{s,ab}}{Z_s},\qquad
    Z_s=\sum_{a,b}K^0_{s,ab}R_{s,ab}.

Then log-odds against the normalized proposal recover
:math:`\log_2R_{s,ab}` up to a stratum constant, which cancels in the
diagonal-zero distance. This explicitly separates the codon proposal,
its amino-acid frequencies and recognition compatibility. A strict
one-nucleotide proposal cannot represent observed exchanges outside its
support; any finite-time codon walk or added proposal component used to
expand that support must be specified before fitting. Model-inferred
nucleotides remain unsuitable as independent evidence for selecting this
coupling.

**One-substitution neighbours require an ascertainment correction.**
When the input contains only Hamming-one neighbours, observing a change
at :math:`p` requires every other compared site to match. Under an
explicitly assumed independent-position baseline with residue
frequencies :math:`f_{e,p,a}`, the opportunity for unordered pair
:math:`a\ne b` is

.. math::


    B_{eg,p,ab}=
    \bigl(f_{e,p,a}f_{g,p,b}+f_{e,p,b}f_{g,p,a}\bigr)
    \prod_{r\ne p}\sum_c f_{e,r,c}f_{g,r,c}.

For within-epitope independent draws this is multiplied by
:math:`\binom{n_e}{2}`; for distinct epitopes by :math:`n_en_g`. These
factors refer to the declared sampling unit, not automatically to source
records. In the origin-aware model, replace this approximation by the
corresponding scenario-conditioned opportunity, including the N-origin
requirement. A generative baseline and observed counts must undergo the
identical neighbour selection. One-position amino-acid frequencies alone
do not correct for matching at every other position. Position
independence is an approximation, not a property of V(D)J rearrangement.

For same/different-epitope counts in a fixed stratum :math:`s`, a
conditional comparison is

.. math::


    O^+_{s,ab}\mid O^+_{s,ab}+O^-_{s,ab}
    \sim\operatorname{Binomial}(O^+_{s,ab}+O^-_{s,ab},\pi_{s,ab}),

.. math::


    \operatorname{logit}\pi_{s,ab}
    =\log\frac{E^+_{s,ab}}{E^-_{s,ab}}+\alpha_s+\beta_{ab}.

The exposures :math:`E^+,E^-` include compartment and neighbour
ascertainment. Here :math:`\beta_{ab}` measures shared recognition
enrichment and :math:`\alpha_s` absorbs the declared stratum intercept.
A common proposal factor cancels only within a correctly matched
stratum. This density ratio alone cannot identify an unrestricted
baseline exchange process. Pairs sharing receptors are dependent;
independent-pair confidence intervals do not follow from this working
likelihood.

**Residue and position terms.** The requested distribution is explicitly

.. math::


    P(a,p\mid\text{exchange endpoint},L,C,NNN),

and is distinct from the probability that a site has at least one
neighbour. A position-independent residue factor is justified only if
this joint table factorizes. For a source-site exposure
:math:`v_p f_{pa}` and relative escape intensity :math:`h_{pa}`, the
endpoint distribution is proportional to

.. math::


    v_p f_{pa}h_{pa}.

It factorizes into residue and position terms exactly when this
nonnegative table has rank one. Even stochastic N can violate this
condition through sequence context, codon phase and selection. A
neighbour incidence per exposed site may exceed one and is an intensity,
not a per-site probability.

**An identifiable exchange-intensity fit.** Let :math:`s` denote the
declared chain, length, position, gene and origin stratum; let
:math:`E^0_{s,ab}` be its baseline pair opportunity after applying
exactly the observed neighbour rule. A working count model is

.. math::


    O_{s,ab}\sim\operatorname{Poisson}(\mu_{s,ab}),\qquad
    \mu_{s,ab}=E^0_{s,ab}\exp(\alpha_s+\beta_{ab}).

Here :math:`O` counts observed pairs, :math:`\alpha_s` sets stratum
density and :math:`\beta_{ab}=\beta_{ba}` describes relative exchange
enrichment. Fix the additive gauge, for example by
:math:`\sum_{a<b}\beta_{ab}=0` over the declared supported cells. With a
single unrestricted cell multiplier and no shared intercept, the
likelihood equation is :math:`\widehat{\exp\beta}_{ab}=O_{ab}/E^0_{ab}`;
shared effects instead solve :math:`\sum_s(O_{s,ab}-\mu_{s,ab})=0`. Thus
raw mismatch frequency is not the fitted compatibility: it must be
divided by its opportunity. A structural zero opportunity is excluded
from this likelihood and retained in the support table; a positive
observation there diagnoses a model, origin or alignment mismatch.
Fractional origin-weighted counts use this as a weighted estimating
equation, not as literal independent Poisson observations. Neither this
working likelihood nor fractional counts justify independent-pair
uncertainty, because receptors and studies are reused. Position-only or
residue-only factors are restrictions on these effects, to be compared
with the joint effect rather than assumed from raw endpoint frequencies.

**PAM-style compatibility process.** Origin-weighted N-region exchange
counts can define a reversible compatibility walk. Let
:math:`A_{ab}=A_{ba}` be off-diagonal exchange counts and :math:`f_a>0`
the residue exposure from the same eligible N-origin sites, with
identical weighting and coordinate scope. Set

.. math::


    Q_{ab}=c\,A_{ab}/f_a\quad(a\ne b),\qquad
    Q_{aa}=-\sum_{b\ne a}Q_{ab}.

Choose :math:`c` to state a mean rate unit, and verify
:math:`f_aQ_{ab}=f_bQ_{ba}`, nonnegative off-diagonals and zero row
sums. When :math:`\sum_a f_a=1`, setting
:math:`c=(\sum_{a\ne b}A_{ab})^{-1}` gives :math:`-\sum_a f_aQ_{aa}=1`:
one expected compatibility jump per unit walk time. The relative escape
intensity is :math:`-Q_{aa}=c\sum_{b\ne a}A_{ab}/f_a`, which separates
residue abundance from its exchange count. Reversibility follows from
symmetrizing unordered convergence pairs; it is imposed by this
construction, not evidence of an evolutionary equilibrium. Then
:math:`P(t)=\exp(tQ)` gives :math:`q_{ab}(t)=f_aP_{ab}(t)`, from which
the same log-odds distance follows. Neither the walk nor :math:`t` is an
ancestral mutation process or biological elapsed time. Hamming-one
counts alone do not identify holding probabilities or a unique
finite-time distance. Reference comparisons therefore state a fixed
time/information convention and distinguish logarithmic distances from
linear genetic-code propensities.

**Indels are changes of rearrangement histories.** Let :math:`G_k(x,y)`
indicate that the existing alignment convention accepts exactly one
contiguous amino-acid gap of size :math:`k`, with shorter length
:math:`L`. For an origin channel :math:`z`, the baseline opportunity is

.. math::


    \begin{aligned}
    E^0_{z,L,C,k} &\propto
    \sum_{r,r'}P_0(r\mid C)P_0(r'\mid C)
    \mathbf1\{|T(r)|=L,\ |T(r')|=L+k\}\\
    &\qquad{}\times
    \mathbf1\{\operatorname{channel}(r,r')=z\}
    G_k(T(r),T(r')).
    \end{aligned}

The proportionality factor is the declared sampling opportunity.
Separate N-length changes, trimming/P boundary changes and gene-choice
changes. Internal deletion of a fixed retained germline block is
forbidden in this model; trimming at a permitted gene end is a different
event. An amino-acid length change :math:`k` requires a net productive
nucleotide change :math:`3k`, but individual trimming/P/N events need
not themselves be multiples of three. The support of :math:`k` depends
on the scenario and gene context, not just on shorter length and chain.
Zero observed counts do not establish a structural prohibition. Where
several gap placements describe one pair, count the pair once using the
owning aligner convention or an explicit placement posterior. Summing
placements as independent observations overcounts repetitive sequence.
Fit recognition multipliers against these opportunities, with the same
exposure-offset logic as substitutions.

To define the requested gap-size probability, normalize the
opportunities only after specifying the conditioning:

.. math::


    P_0(k\mid L,C,\text{one-gap neighbour})=
    \frac{\sum_z E^0_{z,L,C,k}}
         {\sum_{j\ge1}\sum_z E^0_{z,L,C,j}}.

Here :math:`j` indexes allowed positive gap sizes, not a J gene; the
denominator must be positive. This conditional size distribution does
not identify how often one-gap neighbours occur among all eligible
pairs. The latter requires the total pair exposure, including pairs
without a gap event. For observed counts :math:`I_{z,L,C,k}`, fit
:math:`\nu_{z,L,C,k}=E^0_{z,L,C,k}\exp(\gamma_{L,C}+\delta_{z,k})`, with
an explicit gauge for :math:`\delta`, and normalize :math:`\sum_z\nu` to
obtain the recognition-conditioned size distribution. The structural
allowed set is :math:`\{k:\sum_zE^0_{z,L,C,k}>0\}` under the specified
model; distinguish this from finite sampling and floating-point
underflow. There is no universal shortest-length-only list of allowed
gaps when V/J choice, trimming bounds and insertion-length support
differ.

**Generalization across epitopes.** Define the target epitope mixture
before estimating a general matrix. Pair pooling estimates the sampled
event mixture; equal-epitope weighting estimates an equal mixture of
contributing epitopes. Neither automatically removes
sequence-composition or study-selection bias. A weighted matrix is
unbiased only relative to a stated target and a justified
observation/exposure model with adequate support. Epitope-specific
residual effects are assessed after origin, composition and
ascertainment adjustment, with independent source or epitope validation.
Separate target-free estimation, target-guided representational fits and
annotation accuracy. A high correlation created by reference-fitted
weights is not independent recovery of a general substitution law.

**Algorithm and implementation contract.**

1.  Fix the cohort, chain/species, junction coordinates,
    record/clonotype units, gene evidence, ambiguity handling and target
    epitope mixture. Preserve unsupported records and complete
    provenance.
2.  Use ARDA for germline alignment and candidate evidence, and the
    reconciled vdjtools model for trimming/P/N scenarios. Do not create
    an analysis-side alignment or scenario engine.
3.  Compute origin probabilities by the existing sum-product
    recurrences. For the N core, evaluate the exact joint mass requiring
    the entire queried codon to avoid all retained/P-extended templates,
    then divide by the total compatible AA mass. Include single- and
    tandem-D branches.
4.  Keep immutable templates, trimming/P boundaries, pure N, mixed
    codons and unresolved origin in separate accounting. An N-only API
    leaves its complement unresolved among non-N categories.
5.  Run one existing native neighbour/alignment batch per population.
    Attach paired origin weights to substitution and gap events, without
    changing receptor identity, gene calls or epitope labels.
6.  Build origin-matched exposure tables under the identical geometry.
    Include N-origin self-pairs for block matrices and account for
    Hamming/gap ascertainment for neighbour estimators.
7.  Fit the declared shared exchange/position and gap-size terms, with
    explicit identifiability constraints and reference-free priors. Keep
    baseline opportunity distinct from recognition enrichment.
8.  Convert similarities to the actual scorer’s distance units, zero
    diagonal, rounding and clipping. Explicit cost tables are loaded
    directly, without a second Gram transform.
9.  Validate origin normalization, fixed-template exclusion,
    P/mixed-codon handling, ambiguous histories, tandem-D completeness
    and serial/parallel equality. Validate independent source/epitope
    transfer and all 190 off-diagonal distances; report every
    unsupported cell.
10. Evaluate the core and the aa/position/length plus gap/length
    additions on the same fixed accuracy tasks. Promote complexity only
    from added accuracy under the declared comparison, then regenerate
    authoritative tables, figures and manuscript summaries.

**Implementation status and limits.** The exact pure-N constrained sum
and public ``n_region_posterior_batch`` interface are provided by
vdjtools 5.7.0, reusing the native sum-product recurrences with single-
and tandem-D histories. Its release checks cover P-containing
boundaries, ambiguous histories, fixed-template exclusion,
zero-likelihood observations and serial/parallel equality. The VDJmatch
``scoring-audit --task origin-matrix`` command implements
posterior-weighted same-epitope/V/J/length block counts, pooled and
equal-epitope matrices, and source/support/exposure tables. This first
estimator uses unique junction/V/J clonotypes as unit clusters; it does
not claim BLOSUM62 redundancy clustering. The interface does not assign
a full separate G/P/mixed posterior; the full origin-conditioned indel
extension and independent empirical validation remain open. Previously
reported interior/all-residue fits remain descriptive composition
diagnostics, and the accepted core/bundle accuracy comparison retains
its original producing scope. No completed independent high-:math:`R^2`
recovery is claimed for the origin-aware estimator.
