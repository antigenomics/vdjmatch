.. _`sec:selection-information`:

Selection evidence as information
=================================

**Scope and status.**

The objective is to identify which receptor features distinguish a pMHC-specific repertoire from the declared control, and use those features to choose the scoring rule. The equations below define this objective and its calibration contract; they do not promote a new default or establish superiority over other predictors. The V-conditioned masked-four-mer score and the reference/control gene–length information decomposition are implemented diagnostics. A complete information-directed selector and its whole-cohort calibration remain to be established. The original/general score is retained when the available features do not support a more specific rule.

Likelihood ratios, selection factors and bits
---------------------------------------------

Fix species, locus and full peptide–MHC restriction, denoted by context :math:`c`; write :math:`e` for the target within that context and suppress :math:`c` in the equations. Let :math:`x=(L,V,J,Z)` describe an anchor-inclusive junction, its length :math:`L`, gene calls :math:`V,J`, and remaining sequence features :math:`Z`. Define :math:`P_e` as the distribution in the observed target repertoire and :math:`P_0` as the declared background distribution, using the same sampling unit on both sides. For example, distinct junction/gene keys on both sides measure repertoire composition; reference keys divided by control event frequencies instead compare two explicitly different sampling schemes. The empirical selection factor and pointwise information are

.. math::
    Q_e(x)=\frac{P_e(x)}{P_0(x)},\qquad
    i_e(x)=\log_2 Q_e(x).

A positive :math:`i_e(x)` denotes enrichment, a negative value depletion, and one bit doubles the likelihood ratio. For known distributions, :math:`Q_e` is the likelihood ratio for the simple hypotheses :math:`P_e` versus :math:`P_0`. If the target prior in an explicitly defined two-population mixture is :math:`\pi_e`, Bayes’ rule gives

.. math::

   \frac{\Pr(e\mid x)}{\Pr(0\mid x)}
    =\frac{\pi_e}{1-\pi_e}\,Q_e(x).

Without that sampling model and prior, neither an empirical ratio nor a calibrated tail is a binding posterior. Observed repertoire enrichment includes biological selection, sampling and ascertainment; it does not by itself isolate an antigen-specific causal selection factor. In particular, a selected AIRR control is not the recombination distribution :math:`P_{\rm gen}`.

The mean information per sampled receptor is the Kullback–Leibler divergence,

.. math::
    I_e=\sum_x P_e(x)\log_2\frac{P_e(x)}{P_0(x)}
    =D_{\rm KL}(P_e\Vert P_0)=H(P_e,P_0)-H(P_e),

where :math:`H(P_e)=-\sum_xP_e(x)\log_2P_e(x)` and :math:`H(P_e,P_0)=-\sum_xP_e(x)\log_2P_0(x)`. All these quantities are in bits per receptor under the declared sampling unit. A low target entropy alone does not show antigen-specific selection: a concentrated feature can be equally concentrated in the control. These identities assume :math:`P_e` has no positive mass where :math:`P_0` is zero. A zero empirical control count is unresolved support, not evidence that the population probability is exactly zero. Report the affected counts and leave the comparison unresolved; silently dropping those categories, renormalising them away or treating them as infinite evidence changes the question.

**Sampling likelihood and mutual information.**


For :math:`n` independent observations categorised by a fixed feature map, let :math:`n_y` be the category counts and :math:`\widehat p_y=n_y/n`. With a known background probability :math:`p_{0y}`, the multinomial log-likelihood ratio of the fitted categorical distribution to that background is

.. math::

   \log_2\frac{\mathcal L(\widehat p;n_y)}{\mathcal L(p_0;n_y)}
    =\sum_y n_y\log_2\frac{\widehat p_y}{p_{0y}}
    =nD_{\rm KL}(\widehat p\Vert p_0).

The multinomial coefficient cancels; zero-count terms contribute zero. This connects descriptive information to statistical likelihood, but is not a P-value: the alternative was fitted to those counts, the background is usually estimated, and donor/expansion dependence can violate independent sampling. For a binary population label :math:`U` with target probability :math:`\pi`, set :math:`M=\pi P_e+(1-\pi)P_0`. Then mutual information is

.. math::

   I(U;Y)=\pi D_{\rm KL}(P_e(Y)\Vert M(Y))
          +(1-\pi)D_{\rm KL}(P_0(Y)\Vert M(Y)).

It measures label information under that explicit mixture and differs from :math:`D_{\rm KL}(P_e\Vert P_0)`. Neither quantity supplies an automatic performance-optimal weighting of overlapping features.

Which feature carries the information?
--------------------------------------

The probability chain rule gives an exact decomposition, without feature independence:

.. math::

   \begin{aligned}
    i_e(L,V,J,Z)={}&\log_2\frac{P_e(L)}{P_0(L)}
    +\log_2\frac{P_e(V\mid L)}{P_0(V\mid L)}\\
    &+\log_2\frac{P_e(J\mid V,L)}{P_0(J\mid V,L)}
    +\log_2\frac{P_e(Z\mid J,V,L)}{P_0(Z\mid J,V,L)}.
   \end{aligned}

Taking the target expectation decomposes :math:`I_e` into length information and target-weighted conditional divergences for V, J and the sequence remainder. Each expected term is nonnegative, although individual query contributions can be negative. The allocation depends on the conditioning order; it must be reported, not interpreted as a unique causal partition. A V-gene contribution is constant across queries with the same V and conditioning variables: it cannot resolve their junction-specific recognition.

More generally, for a deterministic feature map :math:`Y=\phi(x)`,

.. math::
    D_{\rm KL}(P_e(X)\Vert P_0(X))
    =D_{\rm KL}(P_e(Y)\Vert P_0(Y))
    +\mathbb E_{Y\sim P_e}D_{\rm KL}(P_e(X\mid Y)\Vert P_0(X\mid Y)).

This follows by factoring each probability into its feature marginal and conditional remainder. Thus a projection retains some of the discriminating information and loses the nonnegative conditional remainder. It motivates comparing V/J, length, motif and neighbourhood projections against the same control. It does not justify using their marginal divergences as additive score weights: correlated features can encode the same information. Finite-sample plug-in divergences are descriptive estimates, not unbiased information estimates or proofs that a rule will generalise.

For a genuine pair :math:`(a,b)`, define :math:`C_j(a,b)=P_j(a,b)/[P_j(a)P_j(b)]` for :math:`j\in\{e,0\}`. Then

.. math:: Q_e(a,b)=Q_e(a)Q_e(b)\frac{C_e(a,b)}{C_0(a,b)}.

This separates marginal chain evidence from excess pairing dependence. Independent random TRA/TRB pairing sets :math:`C_0=1` in the specified null; it does not assert independence in natural paired repertoires. When one chain is absent, use its observed partner’s marginal evidence; an unobserved alpha-chain germline signal cannot be assigned to a beta-only query.

Masked four-mers as an information-retrieval score
--------------------------------------------------

The implemented ``motif-evidence`` diagnostic uses the MF–IKF recipe; its equations below specify the calculation without requiring that precursor code. Map amino acids to the eleven classes AILV, RHK, CM, ST, DE, NQ, G, F, P, W and Y. Slide a four-residue window along the full junction and replace each of its four positions in turn by a wildcard. Let :math:`\mathcal T(x)` be the resulting set of distinct tokens; repeated occurrences of one token in a junction count once. V alleles are pooled to gene identity; the reduction acts on amino acids, not on V-gene names.

For token :math:`t`, let :math:`x_{lve}(t)` count distinct reference junctions with length :math:`l`, V gene :math:`v` and full-pMHC label :math:`e` containing :math:`t`. Write :math:`x_{ve}=\sum_lx_{lve}`, :math:`x_e=\sum_{l,v}x_{lve}` and :math:`x(t)=\sum_e x_e(t)` over the declared labelled reference universe of :math:`X` keys. A key includes its target label, so a receptor with several labels contributes to several labelled keys. The fixed support set is

.. math:: \mathcal S_{ve}(t)=\{l:x_{lve}(t)\ge2,\ x_{ve}(t)\ge3,\ x_e(t)\ge5\}.

For ordinary annotation the score is

.. math::
    S_e^{\rm motif}(q)=\frac{1}{L_q}
    \sum_{t\in\mathcal T(q):\,\mathcal S_{V_qe}(t)\ne\varnothing}
    \underbrace{\frac{x_{V_qe}(t)}{x(t)}}_{\mathrm{MF}}
    \underbrace{\ln\frac{X}{x(t)}}_{\mathrm{IKF}}.

The implementation takes the maximum supported length-stratum contribution per token before summing; the displayed unfiltered weight is identical across those strata. MF is the target/V share of token-bearing labelled keys, not the usual within-target term frequency :math:`x_{ve}/N_{ve}`. IKF is the surprisal of token presence under the labelled reference-key distribution, in nats; dividing by :math:`\ln2` expresses it in bits. This is a TF–IDF analogue, not a fitted likelihood ratio against AIRR control. Rare tokens receive more weight only when supported and concentrated in the target/V combination. Overlapping windows are dependent; their sum is a retrieval score, not a product of independent probabilities. Its support cutoffs are inherited implementation choices, not biophysical constants or a newly fitted parameter grid.

In benchmark mode only exact-hit contributions are filtered from the numerator: replace :math:`x_{V_qe}(t)` there by the remaining hit count. The support set, :math:`x(t)`, :math:`X`, priors, graph and all fitted weights remain fixed. No-hit token contributions vanish; no reference model is refitted for an individual query. The motif score can expose a useful sequence axis, but selecting it for an epitope still requires evidence that this axis differs from the declared control.

Control-relative information in a motif projection
--------------------------------------------------

For a fixed masked token :math:`t`, define the binary feature :math:`Z_t(x)=\mathbf1\{t\in\mathcal T(x)\}`. Within a V-gene/length stratum :math:`(v,l)`, let :math:`n,N` be reference keys with the motif and all reference keys, and :math:`m,M` the corresponding control counts. Each population uses distinct full-junction/V-gene keys; repeated windows, V alleles and alternative J annotations do not multiply a key. With :math:`p=n/N` and :math:`q=m/M`, the presence contribution and binary information are

.. math::

   \begin{aligned}
    i_t^+(v,l)&=\log_2\frac{p}{q},\\
    D_t(v,l)&=p\log_2\frac{p}{q}+(1-p)\log_2\frac{1-p}{1-q}.
   \end{aligned}

The first can be negative for a depleted motif; the second is the Kullback–Leibler divergence of one Bernoulli projection, in bits per conditional junction/V key. Zero target-probability terms contribute zero by continuity. An empty control stratum, unobserved control motif with positive target mass, or unobserved control absence with positive target absence leaves the corresponding empirical comparison unavailable; the counts remain explicit. The all-length scope conditions on V only, while the exact-length scope also controls the opportunity to contain more windows in a longer junction. J usage and motif position are not additional conditioning variables in this diagnostic.

The shipping ``reference-structure --control`` command writes these reference-observed token projections to ``.motifs.tsv``. It reuses the existing tokenizer in bounded blocks and reduces once-per-key presence counts, with one control model per species/locus and one fixed model per target. Tables are written one target at a time without retaining every target’s motif rows. Overlapping tokens are dependent: their binary divergences are not an additive estimate of joint sequence information. These descriptive projections help identify informative families and motifs; the retained MF–IKF retrieval score, final rule selection and complete-score calibration remain distinct calculations.

A depleted germline category can still contain a strongly selected junction
solution. Let :math:`N_{\rm all},M_{\rm all}` count all reference/control keys in
the same declared populations. The joint presence evidence factors exactly as

.. math::

   \log_2\frac{n/N_{\rm all}}{m/M_{\rm all}}
   =\underbrace{\log_2\frac{N/N_{\rm all}}{M/M_{\rm all}}}_{\text{V/length stratum}}
   +\underbrace{\log_2\frac{n/N}{m/M}}_{\text{conditional motif}}.

Only matched population units and the same stratum permit this addition. A
pooled marginal V factor cannot replace the V/length term in an exact-length
comparison. A positive conditional term can outweigh stratum depletion;
neither term alone is a receptor veto. This is one joint feature projection,
not a sum over overlapping motifs or a new score-selection rule. Unobserved
control support remains unavailable instead of receiving an invented floor.

Germline-loop information with allele ambiguity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``reference-structure --control`` also writes ``.germline-map.tsv``, ``.germline.tsv`` and ``.germline-information.tsv``. These project the existing V-gene counts onto individual aligned CDR1, CDR2 and CDR2.5 columns, using full amino acids and the eleven IMGT classes separately. Column numbers are one-based columns in the bundled alignment, not structural residue numbers. A gene has a confirmed residue only where every bundled allele agrees. Missing genes, missing loop templates, stop/unknown symbols and allele-dependent columns remain explicit; an omitted allele never selects ``*01``. A dot denotes aligned absence within a populated loop, while an entirely empty template is unavailable.

The sampling unit remains a distinct full-junction/V-gene key. Keys from two genes retain their separate mass after those genes map to the same residue, even if their junctions agree. Let :math:`w_e(g)` and :math:`w_0(g)` count such keys for gene :math:`g`, and let :math:`N=\sum_gw_e(g)` and :math:`M=\sum_gw_0(g)`. For one loop column and alphabet, let :math:`\mathcal A_g` be the set of residues across the modeled alleles; missing models or template symbols admit the entire alphabet. For residue :math:`a`, define

.. math::

   n_e(a)=\sum_gw_e(g)\mathbf1\{\mathcal A_g=\{a\}\},\qquad
   u_0(a)=\sum_gw_0(g)\mathbf1\{a\in\mathcal A_g\}.

For any allele assignments within those allowed sets, the actual residue fractions obey :math:`p_e(a)\ge n_e(a)/N` and :math:`p_0(a)\le u_0(a)/M`. Thus, when the displayed counts are positive,

.. math::

   \log_2\frac{p_e(a)}{p_0(a)}
   \ge \log_2\frac{n_e(a)/N}{u_0(a)/M}.

``residue_bits_lower`` is this lower bound: every compatible or unresolved control gene can contribute to its denominator. It is an empirical bound over allele assignments under the declared template model, not a sampling confidence interval, population guarantee or P-value. Zero compatible control counts leave it unavailable instead of assigning infinite evidence. The exact-length scope uses the corresponding length-stratum counts.

``consensus_bits`` and ``consensus_information_bits`` instead describe the deterministic gene-consensus annotation, with unresolved reference mass left unavailable. Where that expectation is defined, the projection cannot carry more information than the original V-gene distribution on the same key population. This annotation information differs from the unknown allele-resolved residue information; the lower-bound field handles compatible background ambiguity explicitly. Correlated residues, loops and reduced alphabets are separate views and must not be summed as independent evidence. None of these diagnostics selects a score automatically.

Geometry as a projection of selection
-------------------------------------

For a fixed dissimilarity :math:`d` and radius :math:`r`, define :math:`B_d(q,r)=\{x:d(q,x)\le r\}`. Local enrichment is another likelihood ratio, now for the binary neighbourhood event:

.. math:: Q_{e,d,r}(q)=\frac{P_e(B_d(q,r))}{P_0(B_d(q,r))}.

For a nonnegative kernel :math:`K`, the corresponding ratio is

.. math::

   Q_{e,d,K}(q)=
    \frac{\mathbb E_{X\sim P_e}K(d(q,X))}
         {\mathbb E_{X\sim P_0}K(d(q,X))},

when the denominator is positive. These are population definitions of candidate geometry evidence, not identities with every shipped reducer. Distances, radii, kernels and sampling units must match on both sides. A smooth substitution score can retain within-neighbourhood variation; a discrete motif or edit graph can expose separated recognition families. Reduced alphabets can collapse distinct junctions to zero distance, and an alignment dissimilarity need not define a Euclidean manifold.

For :math:`n` reference keys, the empirical pair-mass curve and its scale-dependent slope are

.. math::

   C_e(r)=\frac{2}{n(n-1)}\sum_{i<j}\mathbf1\{d(x_i,x_j)\le r\},
    \qquad D_e(r)=\frac{d\log C_e(r)}{d\log r}.

The slope describes local growth only where the curve has an estimable nonzero range; discrete plateaux, small samples and disconnected families do not yield a single intrinsic dimension. Compare target and control curves under the same geometry and units before attributing density to antigen selection. Nearest-neighbour distances :math:`u_i=\min_{j\ne i}d(x_i,x_j)` give a complementary sorted curve. A Kneedle diagnostic locates a change in its rescaled shape; its knee is a geometric scale proposal, not a significance threshold or a confidence probability. A constant curve has no identified knee. The existing nearest-distance diagnostic and the pair-mass slope are distinct summaries; the latter definition does not claim that it has already been computed.

At a declared radius, graph components with masses :math:`p_k=n_k/n` have entropy :math:`H_C=-\sum_kp_k\log_2p_k`. This describes fragmentation, not antigen information without a matched control comparison. The cross-length nearest-neighbour fraction at query length :math:`l` is

.. math::

   g_l=\frac{\#\{i:L_i=l,\ \min_{j:L_j\ne l}d(x_i,x_j)
                         <\min_{j\ne i:L_j=l}d(x_i,x_j)\}}
                {\#\{i:L_i=l,\ \text{both comparisons available}\}}.

Strict ties are not cross-length wins; an empty denominator is unavailable. Unlike :math:`1-\sum_lP_e(L=l)^2`, this statistic measures an actual geometric relation, not length diversity alone. Neither statistic is, by itself, a validated gap-score weight.

Human-readable rules and comparable tails
-----------------------------------------

A germline-anchored class has informative V/loop usage relative to control, a motif-constrained class has informative conditional junction motifs, and a distributed-neighbourhood class has informative smooth local enrichment. Multiple recognition families can coexist within a target; these labels describe measured axes rather than force each peptide into one exclusive class. MHC class, allele, peptide chemistry and available structural contacts motivate an interpretation of those axes; peptide composition alone does not determine receptor recognition or its score. Species, chain coverage, length and family support remain explicit conditions. No informative supported projection means use the retained general score and report the available evidence, rather than invent specificity.

Let :math:`A` denote the complete procedure, including data-dependent rule choice, and :math:`T_e(q,R)=T_{A,e}(q,R)` its higher-is-better score. The two directional tails are

.. math::

   \begin{aligned}
    p_{R,e}&=\frac{1+\sum_{b=1}^{B_R}\mathbf1\{T_e(q,R^{*(b)})\ge T_e(q,R)\}}{B_R+1},\\
    p_{Q,e}&=\frac{1+\sum_{b=1}^{B_Q}\mathbf1\{T_e(q^{*(b)},R)\ge T_e(q,R)\}}{B_Q+1}.
   \end{aligned}

The reference direction uses composition-matched mock references; the query direction uses control junctions matched to query Pgen strata. A stratum conditions on generation probability, not on exact sequence identity or all biological selection effects. Each mock reference gets its own fixed model and the same rule-selection algorithm. The real reference model remains fixed across queries; the query direction applies the same procedure to every random query. Ties enter the upper tail, all null no-hits remain in the denominator, and observed no-hits receive one. Under the respective conditional exchangeability assumptions, these are conservative Monte Carlo rank P-values. Matching Pgen alone does not prove those assumptions.

If both directional P-values are valid for nulls :math:`H_R,H_Q`, then :math:`p_{{\rm both},e}=\max(p_{R,e},p_{Q,e})` is valid for :math:`H_R\cup H_Q`. Indeed, under :math:`H_R`, :math:`\Pr(p_{{\rm both},e}\le\alpha)\le\Pr(p_{R,e}\le\alpha)\le\alpha`, and the same holds under :math:`H_Q`. No independence assumption is needed. This gives a common tail scale for different epitope scoring rules, conditional on their own correctly reproduced nulls; it does not make the scores themselves universal likelihoods. Selecting the smallest tail across labels or the strongest paired/marginal axis requires calibration of that complete selection event, not reuse of a component tail. The same principle applies when reporting any confident specificity in a repertoire: the reported query/label search scope must match the calibrated event.

With :math:`B=199`, a plus-one tail cannot be smaller than :math:`1/200=0.005`. Many distinct strong raw scores therefore share the same empirical tail. Preserve raw within-epitope ranking alongside cross-epitope tails; finite calibration resolution must not silently erase ranking information. The present formulation supplies exact identities and conditional validity statements, not a guarantee of future ranking accuracy, biological binding or state-of-the-art performance.

Implementation entry points
---------------------------

Use :doc:`cli` for ``motif-evidence`` and ``reference-structure --control``. The latter produces pointwise selection bits in ``.genes.tsv`` and conditional/marginal information in ``.information.tsv``.
See :doc:`scoring` for the current default reducers and :doc:`control-calibration`
for the two control directions. The information-directed selector described here
is a mathematical development objective, not an additional CLI mode.
