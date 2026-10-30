Command-line interface
======================

The console script provides reference management, annotation and benchmark commands. Use
``vdjmatch match -h`` for the installed version's complete options.

vdjmatch update
---------------

Fetch and retain a latest or pinned VDJdb release input, and prepare its adjacent
native index and epitope/pMHC profiles for CLI and API use.

.. code-block:: bash

   vdjmatch update --asset default --cache reference-inputs
   vdjmatch update --pin YOUR_RELEASE_TAG --asset legacy --cache reference-inputs

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Option
     - Meaning
   * - ``--asset``
     - ``default``, ``primary``, ``legacy``, ``slim``, ``full`` or ``airr``;
       select a role or legacy projection supplied by the release.
   * - ``--cache``
     - Directory retaining downloaded immutable reference inputs.
   * - ``--pin``
     - Specific release tag; omit to resolve latest.
   * - ``--force``
     - Fetch the requested release again.

vdjmatch index
--------------

Prepare a fetched or custom reference using the annotation ingestion rules.
No query repertoire, Pgen model or external control is required.

.. code-block:: bash

   vdjmatch index --threads 4
   vdjmatch index --vdjdb reference.tsv --threads 4
   vdjmatch index --vdjdb reference.zip --species HomoSapiens --threads 4

The printed sibling directory ``REFERENCE.vdjmatch/IDENTITY/`` contains saved
native indexes, normalized observations, ``profiles.tsv``, ``profiles.parquet``,
``joint-profiles.tsv``
and a checksummed manifest. Unchanged references reopen these assets without
rebuilding or reclassifying. Source bytes, filters, implementation and native
format determine identity; incomplete or corrupt assets rebuild from the source.
The reference's parent directory must be writable. ``Annotator.from_path`` and
``match`` also prepare these assets automatically. ``annotate`` scores its
reference directly without preparing unused profile sidecars. ``index`` and
``update`` accept a positive ``--threads`` budget (default 4).

Use ``reference-profiles`` for a fresh scientific classification rebuild:

.. code-block:: bash

   vdjmatch reference-profiles --vdjdb reference.zip --threads 4 --output-prefix fresh

It reloads raw observations and rebuilds graphs and label-shuffle diagnostics
on every invocation, without reopening saved profile/index bundles. Outputs are
``fresh.profiles.tsv``, ``fresh.joint-profiles.tsv`` and a manifest with raw,
package-source and actual native hashes checked before and after calculation.
``--shuffles`` defaults to 999 (minimum 200). The shuffle count and recipe
are recorded; earlier 199-shuffle outputs remain historical comparisons.
``--species human`` or ``mouse`` selects one population; omission retains all
represented species. ``--reference-population`` defaults to curated benchmark
validity; ``all-valid`` retains all valid junction records. A new output prefix
keeps completed tables recoverable. This command does not score receptors.

Profiles retain separate species and TRA/TRB chains, at pooled-peptide and
observed full-pMHC levels. Eligibility requires at least 30 valid source records;
distinct-junction count is shown alongside. Classification tests two information
axes: V-usage KL and enriched core 3-mer presence conditional on V. The core trims
three residues from each end of the full junction. Each motif is counted once
per junction/V key; overlapping motifs are not added as independent evidence.

The reference itself supplies the noise comparison: 999 fixed label shuffles
(seed 0). Whole key-membership vectors shuffle globally for V usage and within
V for the motif test, preserving each target's V composition in the latter.
The strongest motif is selected again inside every shuffle. Each axis requires
a Monte Carlo tail probability at most 0.005 and positive information above the
null mean. The two tests give a per-profile error bound of 0.01 under their
respective exchangeability nulls. This is not a database-wide guarantee or
binding-score calibration. See :doc:`algorithms` for definitions.

The decision leaves are ``v_and_motif``, ``v_only``, ``motif_only`` and
``no_supported_axis``. Missing V calls are imputed once per original distinct
junction/V key from resolved-gene frequencies within the same species/locus,
using a fixed seed0 stream independent of epitope labels. Original calls, key
identities and unresolved-key counts/fractions are retained; imputation never
merges original keys. A wholly unresolved population uses one neutral missing
category, providing no invented V-usage signal. Missing V does not veto classes.
Insufficient source records, fewer than two distinct keys or missing
other-peptide background remain explicit.
The selected motif's observed counts, raw/excess bits, null mean/maximum and
Monte Carlo tail are included alongside each decision. V concentration and
neighbour density remain descriptive, without a binary concentration/enrichment
class flag. A small supported effect need not provide
useful binding prediction, and these classes do not change the scoring default.

``length_database_kl_bits`` and ``length_database_js_bits`` describe each
eligible target's exact length law against the complete same-species/locus
database junction/V-key population. They add no support test or class switch.
The separate control-based ``reference-structure`` length fields use
junction/V/J keys and identify their control denominator explicitly.
``joint-profiles.tsv`` retains both marginal classes, each chain's valid
record count and eligibility, and ``any_chain_motif_supported``. A supported
chain suffices; false plus unknown remains unknown. The table marks whether
both chains are observed and eligible, without asserting physical pairing.

The separate ``full20_profile_class`` and ``full20_motif_*`` fields replace
the core motif feature set with the union of full20 3-mers and 4-mers with
one X at each possible position. The three-residue terminal trim is unchanged.
Maximum selection over the entire union repeats inside every one of the
999 V-conditional label shuffles. Its per-axis cutoff is also 0.005. This
alternative preserves the original core3 class and does not combine the two
motif tests into a larger family with the original error bound. Joint profiles
also export ``any_chain_full20_motif_supported``.

Neighbour enrichment is ``(w/Pw)/(b/Pb)``. ``w`` counts unordered distinct
within-target Levenshtein-one pairs; ``b`` counts target-to-background pairs.
``Pw`` and ``Pb`` count nonidentical pair opportunities with length difference
at most one. Background is the distinct-junction union annotated to any other
peptide in the same species/chain, including shared-label junctions. Exact pairs
contribute to neither count nor opportunity; exact sharing is reported separately.
Source duplicates do not increase graph or V-key evidence.

The table also reports the fraction of distinct target junctions with a
within-target neighbour. Missing background, no neighbour evidence or inadequate
coverage remains explicit. These classes describe reference structure; they do
not fit assay labels, select a scoring default or supply binding probabilities.
Pooled modal MHC annotations are distinct from actual full-pMHC profiles;
unlinked chains do not create a paired-receptor class.

vdjmatch annotate
-----------------

The targeted annotation default is the stem4 kernel-ladder scorer:

Aligned CDR1/CDR2/CDR2.5 V-loop costs are precomputed once per annotation
invocation and reused across target populations. Different V genes retain their
loop-similarity contribution; the same-V shell is an additional kernel component.
The legacy ``match`` command retains edit ``--scope`` and optional ``--match-v``.

.. code-block:: bash

   vdjmatch annotate sample.tsv --threads 4 -o annotation
   vdjmatch annotate airr_directory/ --threads 4 -o cohort
   vdjmatch annotate --metadata samples.tsv --vdjdb reference.zip --threads 4 -o cohort

The reference is fetched automatically unless ``--vdjdb`` supplies an offline
table/archive. Metadata has ``path`` and optional ``sample_id``; relative paths
resolve beside the metadata TSV. Directory discovery accepts TSV, CSV, text,
gzip and Parquet repertoire tables with unique sample names.

A custom reference can contain ``junction_aa,locus,epitope`` with ``pair_id``
linking paired chain rows, or wide ``junction_aa_alpha/beta,v_alpha/beta,j_alpha/beta``
and ``epitope`` columns. Reference pairing and source metadata survive normalization.

Input uses anchor-inclusive ``junction_aa`` with AIRR ``v_call`` or simple ``v``
and optional ``j`` fields. Wide pairs use ``junction_aa_alpha/beta`` and
``v_alpha/beta``, with optional ``j_alpha/beta``. Cell/pair identity takes priority
over clone linkage; clone and clonotype IDs remain distinct. Missing chains are
retained and ambiguous multiple chains require resolution before scoring.
``--locus auto`` infers populations; ``TRA``, ``TRB`` and ``paired`` select explicitly.
The command retains exact evidence in ordinary annotation. ``--exclude-exact``
filters exact hits for benchmarking; model populations, priors and graph weights stay fixed. Combine it
with the default ``--reference-population curated-benchmark`` for the manuscript
primary reference view. Use ``--reference-population all-valid`` to reproduce the
accepted all-valid baseline; the stem4 formula is unchanged.

One-control complete-score mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With the default ``--background vdjdb-other``, ``--control RAW_TABLE`` opts into
combined-score calibration. The control supplies full junction/V/J events;
productive and amino-acid eligibility are reported separately. Use
``--mock-replicates N`` (default200), ``--seed N`` (default0), repeatable
``--evidence-type TYPE`` and ``--min-score N`` (default0). Reference filters apply
before the null populations are created. Query and control must be different
inputs. Unsupported nondefault geometry/selection options raise an error.

The outputs are ``.queries.tsv``, ``.database.tsv``, ``.epitopes.tsv`` and
``.calibration.tsv``. No JSON or fitted calibration input is needed. Every raw
chain row survives with a generated calibration identity and disposition;
genuine incomplete pairs route eligible chains to single-chain calibration.
Database tails include all-label and available-axis maxima. Fixed full-pMHC
rows retain components, weights, directional tails and separate null match rates.
See :doc:`control-calibration` for the conditional laws and candidate prior policy.
``--background real``/``generative`` retain their existing control behavior.

Default stem4 score
~~~~~~~~~~~~~~~~~~~

The geometry combines threefold 3/3-trimmed junction distance, aligned V loops,
and the rounded mean of four hairpin-level costs on untrimmed junctions. All
terms use unclipped BLOSUM62 Gram penalties. Each level uses the maximum of its
two opposite-residue costs; the stem is not multiplied by the matrix maximum.

For target population size ``N``, shell counts ``dn(t)``, same-V shell counts
``dv(t)`` and other-peptide cumulative counts ``m(t)``, the score is
``sum_t exp(-t/h) * (0.25*dn(t) + 0.75*dv(t)) / (1+m(t)) / N``.
Here ``h = 22.5 * Pmax`` and ``Pmax`` is the maximum substitution penalty.
The fixed ladder contains 90 inner intervals to ``h`` and 62 outer intervals
through ``100 * Pmax``. Empty target/control populations remain unavailable.
Kernel ``n_reference`` and ``n_control`` columns count through this full ladder
reach; the manifest records their cutoff.
Explicit peptide/MHC selectors restrict targets while retaining the full
other-peptide reference background. By default every full peptide–MHC target
is scored; query assay labels are ignored.

Default output includes query dispositions, selection TSVs and a manifest.
Every eligible target is scored; restrictions of each peptide are reduced before
retaining the winner, strongest competing peptide, score ties and exact-label
ambiguity. Tables retained in memory grow with query count rather than the
query-by-target product. ``--details`` additionally writes candidate and full
all-label enrichment TSVs for continuous accuracy or detailed inspection.
The advanced ``neighbours`` kernel route enables compact output with
``--selection-only``; its default detailed exports remain unchanged.
``evidence --evidence-only`` exports all-label enrichment, selection and query
dispositions (plus linked and full-reference marginal evidence for pairs), while
omitting candidate distance matrices. It requires uniform ``vdjdb-other`` kernel
scoring with non-null peptide targets and excludes other diagnostic reducers.
Scores and fixed models are
identical to detailed output; this option adds no calibration. Add
``--gzip-evidence`` to write enrichment and paired-evidence as ``.tsv.gz`` with
native gzip level1; query dispositions and selection remain plain TSVs.
Selection pools restrictions of one peptide, retains score ties and conflicting
exact-label ambiguity, and reports winner/(winner+strongest other peptide) as
``kernel_support_share``. Kernel scores and shares are relative evidence. Use
``calibrate-calls`` with matching whole-selection null scores before making a
false-call recommendation; component calibration does not transfer.

Paired runs emit linked alpha+beta evidence and separate ``.TRA``/``.TRB``
outputs using all compatible single-chain references, including unlinked ones.
With ``--details``, ``.paired-evidence.tsv`` keeps those three scores together. Selection reports
three component labels/statuses and ``has_support``; missing linked references
cannot veto a supported marginal. There is no inferred joint probability.

The original PSSM-density scorer remains available as ``historical-density``. Existing
advanced alternatives remain in ``historical-density`` and ``neighbours``;
see their installed ``--help`` for options. Enrichment exports all/V-weighted/
same-V-only ladder scores with bounded and rare-control denominators. These
alternatives retain their own identities and are not calibrated by the default.
``--kernel-scorer`` selects ``v-bounded`` (default), ``all-bounded``,
``vonly-bounded``, ``v-rare``, ``all-rare`` or ``vonly-rare`` for selection.
Historical V-similarity ratios are also precomputed per invocation, preserving
query/reference direction and the original unknown-gene behavior.
Mixed compact single/paired cohorts share marginal chain searches; the V-loop
tables are prepared once for the cohort. ``--v-regions cdr12`` retains the
``vcdr12`` alternative, omitting only CDR2.5. Its recorded geometry also uses
``--stem-depth 0 --distance-matrix clipped --junction-ends trim3``.
Default geometry retains CDR1, CDR2 and CDR2.5. Calibration belongs to the entire
selected geometry and scorer.
The matrix extras are ``--distance-matrix pam250``, ``mj`` and ``vdjam``;
their recorded geometries use ``--stem-depth 0 --junction-ends trim3``.
They retain unclipped Gram penalties and use the same matrix units for V loops,
junctions, gaps and the ladder.
``--position-weighting significance-trimmed`` selects the recorded ``sig``
geometry with ``--distance-matrix clipped --stem-depth 0 --junction-ends trim3``;
``sigb62`` changes the matrix to ``blosum62``. These require compact output.
The full-junction profile is restricted to the retained window, normalized to
mean one and rounded to integer weights; the junction component is floored
before V-loop and stem addition, independently for each paired chain.

``--geometry-extra`` selects the fixed retained ``flank``, ``framev``,
``framemid``, ``framej``, ``freegap``, ``jw``, ``stem4j``, ``jlen``, ``vjlen``,
``norm``, ``jnorm``, ``sw`` or ``nw`` recipe. Each sets its matrix, trim3 ends
and stem term and requires compact uniform full-alphabet kernel output with
all V regions. J recipes compare literal J families after removing allele
suffixes; length recipes use full-junction length. Local/global alignment uses
the native affine implementation. Frame recipes retain their historical
callback ordering. The historical ``norm``/``jnorm`` denominator includes the
pre-exclusion same-V target count; exact-excluded runs therefore reproduce that
diagnostic rather than a completely punctured model. Calibration cannot transfer
from the targeted default to these alternatives.

vdjmatch match
--------------

.. code-block:: bash

   vdjmatch match [options] SAMPLE [SAMPLE ...]

The command preserves query identities and reports all peptide–MHC candidates.
It requests background calibration by default; ``--no-evalue`` selects uncalibrated
ranking. Reference and background loading errors propagate.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Option
     - Meaning
   * - ``-o, --output-prefix``
     - Output path prefix.
   * - ``--vdjdb``
     - Local reference ZIP, extracted directory or supported table; omit to fetch latest.
   * - ``--asset``
     - Reference role/projection: ``default``, ``primary``, ``legacy``, ``slim``,
       ``full`` or ``airr``.
   * - ``--cache``
     - Location for immutable downloaded reference inputs.
   * - ``--pin``
     - Specific release tag.
   * - ``--species``
     - Reference species filter; default ``HomoSapiens``.
   * - ``--scope``
     - Maximum substitutions, insertions, deletions and total edits;
       default ``1,0,0,1``.
   * - ``--matrix``
     - Bundled ``vdjam`` scoring or ``none`` for unit costs.
   * - ``--min-score``
     - Minimum database confidence score, separate from search/ranking scores.
   * - ``--input-format``
     - Source convention selector ``auto``, ``airr``, ``legacy`` or ``custom``.
   * - ``--sequence-convention``
     - ``junction`` declares anchor-inclusive sequences in ambiguous custom tables;
       it cannot convert a bare AIRR CDR3 into a junction.
   * - ``--paired``
     - Read linked TRA/TRB rows and require same-reference-complex paired evidence.
   * - ``--link``
     - Source linkage column for paired rows, such as ``cell_id``.
   * - ``--match-v`` / ``--match-j``
     - Require matching V/J genes. Supply gene-labelled raw ``--control`` tables
       for calibration, or use ``--no-evalue``. Sequence-only controls cannot
       calibrate these predicates. Paired matching rejects
       both flags until chain-specific V/J predicates are supported.
   * - ``--no-evalue``
     - Use explicitly uncalibrated ranking without requesting a background.
   * - ``--no-align``
     - Omit per-hit alignment and CIGAR computation for single-chain and paired output.
   * - ``--threads``
     - Native sequence-search threads; default ``0`` resolves the job allocation
       before execution. ``1`` selects serial operation. This does not reconfigure an already initialized Polars thread pool.

Output
~~~~~~

For each sample, ``match`` writes TSV tables for detailed hits, all candidates,
per-query calls, descriptive summaries and ingestion diagnostics, plus a JSON
manifest. The manifest identifies reference inputs, software, filters, source
convention and resource settings. Sample identities determine output names;
colliding names are rejected rather than overwriting another sample.

The hit table retains independent reference observations and their metadata.
Paired hits are detailed same-complex alpha/beta observation matches, with chain
metadata and optional alignment fields prefixed ``alpha_`` and ``beta_``. Paired
candidates are reduced separately from those observations using distinct junction
pairs; the hits table is not a copy of the candidate table. Both come from the same
reference searches.
Candidates retain competing peptide–MHC restrictions, ranking scores and any
requested background evidence. Calls retain invalid and unmatched rows with
explicit statuses. Ingestion diagnostics state input, missing, invalid, retained
and dropped row counts. Summary counts describe matches; they are not sample-level
statistical enrichment tests.

See :doc:`explanation` for counting units and :doc:`how-to` for paired/custom input
examples. ``score``, ``vdjdb_score``, NED, ``E`` and ``p_enrichment`` represent
separate quantities; none is a posterior specificity probability.

``vdjmatch precursor``
----------------------

T-cell **precursor frequency** and unseen-junction diversity for a set of TCRs — how much
repertoire mass can see an epitope. Needs the optional extra::

   pip install 'vdjmatch[precursor]'

.. code-block:: console

   $ vdjmatch precursor --vdjdb -o precursor.txt              # every VDJdb epitope, both chains
   $ vdjmatch precursor --vdjdb --min-junctions 10 -r 2       # well-sampled epitopes, radius 2
   $ vdjmatch precursor tcrs.tsv --group-by epitope --locus TRA
   $ vdjmatch precursor --vdjdb --q 9.41 --n-eff 1e8          # calibrated F, plus P(>=k precursors)
   $ vdjmatch precursor --vdjdb --species MusMusculus --organism mouse   # mouse (auto-picks arda)

.. warning::

   Sequences must be **junctions** (Cys104…Phe/Trp118 inclusive), not IMGT CDR3s — VDJdb's column
   is named ``cdr3`` but holds junctions. An anchor-stripped CDR3 scores exactly ``0.0`` with no
   error, so it is dropped and counted in ``n_dropped`` rather than silently scored.

.. note::

   **Where the method is described.** This page documents the tool. The estimator, its benchmarks
   and the paper live in separate repositories:
   `repseq/2026-precursor-freq <https://github.com/repseq/2026-precursor-freq>`_ (benchmarks, result
   tables, the 16-study literature compendium of measured naive precursor frequencies) and
   `repseq/2026-precursor-freq-ms <https://github.com/repseq/2026-precursor-freq-ms>`_ (the
   manuscript).

.. note::

   **Which model set.** ``--source`` selects a bundled model family, and the choice matters more
   for a mass over a *set* than for scoring one sequence, because a junction whose ``Pgen`` is
   exactly zero contributes nothing and raises nothing.

   ``olga`` (default)
      A bit-faithful import of OLGA's published models — ``Pgen`` matches OLGA's own to machine
      precision. Human, seven loci. Faithful includes faithful to OLGA's deletion-bin grid, on which
      an allele shorter than the grid carries probability on trims it cannot reach; **5.7% of human
      TRA junctions in VDJdb score exactly zero** as a result, and for some epitopes it reaches 13–15%.
   ``isalgo``
      Current gene-level fits for seven human loci and mouse TRA/TRB/IGH/IGK/IGL.
      ``learned`` and ``arda`` are compatibility aliases to these same weights.
      Mouse TRG/TRD are unavailable and fail explicitly.

   The legacy ``olga`` default aliases ``olga-legacy`` for historical reproduction.
   Historical learned/ARDA fits had 270/52,363 and 19/52,363 zero-Pgen human TRA
   junctions respectively, and two zero-Pgen human TRB junctions each. The old
   mouse ARDA TRA/TRB fits had no zero-Pgen junctions in that accepted cohort.
   These recorded results describe the original weights, not the current alias
   targets; historical reproduction requires the pinned original parameters.
   The ``n_zero_pgen`` output reports the observed loss per group and warns
   above 1% overall.

.. note::

   ``--species`` selects the **records**, ``--organism`` selects the **model**, and they have to
   agree. Scoring mouse junctions against the human recombination model does not error and does not
   produce zeros — it returns a plausible number for every group, a median 0.157× of the right
   answer with an 11× spread across epitopes, so the ranking is wrong too. That mismatch is
   therefore refused. Mouse models use ``isalgo`` (also selected by the ``arda`` alias), which ``--organism mouse``
   selects for you.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Option
     - Meaning
   * - ``--vdjdb``
     - score VDJdb itself instead of an input file (with ``--table`` / ``--pin`` / ``--species`` / ``--mhc-class``)
   * - ``--group-by``
     - column to group by; one output row per group (e.g. ``epitope``)
   * - ``--chain-col``
     - locus column; each chain is scored with its own recombination model
   * - ``--capture-col``
     - capture-unit column (e.g. ``reference_id``) enabling the unseen-species estimate
   * - ``-r`` / ``--radius``
     - neighbourhood radius in substitutions (default ``1`` — part of the estimator, not a knob)
   * - ``--alpha``
     - cognacy retention per edit (default ``0.1``, Mayer & Callan 2023)
   * - ``--q``
     - selection constant; ``1`` = uncalibrated raw mass. ALICE's published TRB value is ``9.41``
   * - ``--n-cells`` / ``--compartment``
     - pool size and the fraction the restriction addresses, for the expected precursor count
   * - ``--n-eff``
     - independent rearrangements; switches on ``P(>=k precursors)`` and the seen/unseen counts
   * - ``--selection``
     - depth factor for the occupancy counts: a float, or ``auto`` for the measured per-chain
       values (TRB 4.62, TRA 1.07). These differ and are not pooled
   * - ``--min-junctions``
     - skip groups with fewer distinct junctions

Output
~~~~~~

One tab-separated table, one row per group. Key columns:

- ``n_dropped``, ``n_zero_pgen`` — junctions rejected by the anchor check, and junctions that
  passed it and still scored ``Pgen`` exactly zero. The second is a property of the model set, not
  of your input (see the note above); both contribute nothing to any mass, so a large count means
  the totals in that row are short.
- ``pgen_log10_span``, ``top1_share``, ``top10_share`` — how far the cognate set's generation
  probabilities spread, and how concentrated the sum is. A mean is the wrong summary here.
- ``naive_sum``, ``union``, ``overlap`` — the per-sequence ball masses added up, the exact union,
  and the share of the naive sum that double-counting would have invented.
- ``retained``, ``F``, ``cells`` — the shell-weighted estimate, the calibrated frequency and the
  expected precursor count at ``--n-cells``.
- ``lambda``, ``p_ge_1``, ``p_ge_10`` — Poisson precursor probabilities (needs ``--n-eff``).
- ``S``, ``n_seen``, ``n_unseen``, ``seen_fraction`` — the occupancy model (needs ``--n-eff``):
  the effective cognate-set size, how many of its clonotypes a repertoire of that depth shows, and
  how many it does not. ``n_seen`` **saturates**, so a cognate set concentrated on a few high-``Pgen``
  junctions is exhausted at shallow depth while a broad one keeps accumulating.
- ``n_ball``, ``n_unseen_ball``, ``unseen_ball_mass`` — the raw neighbourhood census, unweighted by
  cognacy and independent of depth.
- ``n_unseen_ht``, ``unseen_ht_mass``, ``rarity_ratio``, ``richness_reliable`` — the
  Horvitz–Thompson extrapolation. Read the **mass**; the count diverges in the tail and
  ``richness_reliable`` says when that is happening.

Reference selectors
-------------------

Repeat ``--epitope``, ``--mhc-a``, ``--mhc-b``, ``--reference-id``,
``--exclude-reference-ids`` or ``--evidence-type`` to filter the reference.
Values within a selector are alternatives; selectors combine conjunctively.
Study exclusions are explicit and recorded in the manifest.

For AIRR reference exports without species metadata, use ``--species any``.
Calibration then requires ``--control-species human`` or ``mouse``; alternatively
select ``--no-evalue`` for uncalibrated ranking.

For a fresh raw-input run, use ``--fresh-control``. Human TRB is bundled; other
loci require ``--control TRA=control.tsv`` or ``--control TRB=control.tsv``.
Supplied tables explicitly declare junctions from that locus and control species.
Controls are rebuilt once and their raw hashes/counts enter the manifest; supplied
mappings must cover every active locus. No persisted search index is reused.

``--exclude-exact`` punctures zero-edit target and control hits before candidate
reduction in single-chain mode. It preserves reference/control population sizes.
In paired mode it removes only the exact/exact junction pair; an exact alpha
with a neighbouring beta, or the reverse, remains evidence. Target counts use
unique accepted junction pairs, while detailed output preserves accepted source
observations. Reference pair population sizes remain unchanged.

Paired exclusion reports empirical Cartesian ``n_control_joint = na*nb - za*zb``
and ``E_raw`` with ``p_poisson_raw``. Here ``na,nb`` count accepted chain controls
and ``za,zb`` their exact query memberships. Controls must contain unique
junctions. The reported ``E`` uses a symmetric upper bound from four one-sided
Clopper–Pearson marginal category bounds, applied at all counts with total failure
budget ``delta=1e-6``. ``p_enrichment`` is
``min(1, delta + BinomialSF(n_hits-1, n_reference, p_upper))``; this differs from
the default paired Poisson calibration. The calibration label is
``paired_independent_joint_punctured_binomial_bound``. ``rule_of_three`` is false
for this method; ``finite_control_delta`` identifies its bound budget.

This is model-conditional: the query and search predicate are fixed, marginal
control categories and target pairs are IID, background chains are independent,
and controls are independent of targets. Deduplicating raw controls does not
establish IID sampling. Neither this bound nor its tail guarantees unconditional
calibration for biological repertoires, estimates natural chain co-occurrence,
or reports posterior specificity confidence. The manifest records these
assumptions and the fixed budget. Default paired matching without exclusion
retains its existing independent-chain Poisson calculation.

``--verbose`` reports coarse native search stages and per-sample timing/resource
information. It keeps one native query batch per locus; it does not split the batch
to render a progress bar.

Global first-hit calibration
----------------------------

``first-hit`` exposes the existing adaptive first-hit calculation against the
whole selected reference. It reports enrichment relative to a raw control;
this is distinct from peptide–MHC specificity confidence.

.. code-block:: bash

   vdjmatch first-hit sample.tsv --vdjdb vdjdb.txt --locus TRB \
     --scope 5,2,2,5 --min-refs 2 --exclude-exact --threads 4 \
     --output-prefix results/first_hit

The command writes ``.evidence.tsv`` with one row per original query, including
invalid and no-hit rows, and ``.manifest.json`` with input hashes, software
versions, selection counts and resources. Target/control counts use unique
junctions; query identity and abundance remain separate. ``--min-refs`` selects
reference clonotype/epitope combinations supported by distinct studies.

Human TRB uses the raw control bundled with seqtree. Other species/loci require
``--control`` with a matching raw junction repertoire. Every invocation builds
fresh target and control indexes and sends one native query batch to each.
AIRR inputs require ``junction_aa``. Default threads are one.

Sequence batch search
---------------------

``search`` exposes native amino-acid matching for reproducible graph and scaling
calculations. It preserves original query/reference row positions, including
duplicate sequences, and builds a fresh index for one native query batch.

.. code-block:: bash

   vdjmatch search queries.tsv --reference reference.tsv \
     --sequence-col junction_aa --scope 1,1,1,1 --threads 4 \
     --same-key-col v_call --output-prefix results/search

Outputs are ``.pairs.tsv``, ``.counts.tsv`` and ``.manifest.json``. Pair rows report
``query_row``, ``reference_row``, score and available edit counts. ``--counts-only``
omits pair rows; counts include reference duplicates. ``--same-key-col`` requires
equal non-null keys in both tables and filters hits before output construction
and positional rescoring. It does not reduce the native candidate search memory.
``--exclude-exact`` removes sequence identities. Generic amino-acid columns are
accepted; this command does not infer receptor locus or CDR3/junction conventions.

``--engine seqtm`` enforces individual substitution/insertion/deletion caps.
``--engine seqtrie`` uses only the total unit-cost radius and reports null edit
decomposition. Matrix scoring requires seqtm: choose ``none``, ``vdjam`` or
``blosum62``. Scores follow the selected matrix's integer cost units.
``--position-significance`` requires an explicit matrix and zero indel caps;
it performs one unit-cost candidate search, then vectorized rescoring using the
existing end-anchored positional model and integer weights
``max(1, round(100 * significance_weight))``. These weighted scores have different
units from unweighted matrix scores. The manifest records parameters, model/input
hashes, counts, software versions and separate build/search/rescoring times.

Graded neighbourhood annotation
-------------------------------

``--search-mode ball`` reports fixed edit balls through ``--max-radius``
(default 5, permitted values 1–5). ``--max-radius 3`` limits substitutions,
insertions and deletions together to three per chain and reports radii 1–3.
Unit-cost balls share one native batch at the declared maximum. Weighted scoring
uses one complete query batch per radius: its best alignment can use more edits
than the minimum-edit alignment, so a broad weighted hit cannot simply be
filtered into smaller balls. Each ball permits at most two insertions/deletions
and caps total edits at its radius. Matrix penalties rank accepted neighbours.
No global neighbour-count cap truncates the evidence. The calls table has one
row per retained query and radius; it does not select the radius with the smallest
P-value. Ball output can be large and does not include fixed-mode CIGAR strings.

.. code-block:: bash

   vdjmatch match queries.tsv --vdjdb reference.zip --search-mode ball \
     --max-radius 3 --threads 4 --output-prefix results/balls
   vdjmatch match cells.tsv --vdjdb reference.zip --paired --scope 3,1,1,3 \
     --threads 4 --output-prefix results/paired

Paired matching retains its linked-only fixed-mode default. The explicit
``--paired-marginals`` opt-in adds fixed scopes and graded balls for unit,
hard-V unit, weighted and hard-V weighted families. Each chain receives the
same declared budget; paired ball radius is ``max(edits_alpha, edits_beta)``.
The linked axis retains the existing unique-junction-pair exponential mass;
TRA/TRB axes retain existing NED ranking over all references, including unlinked
observations. Label assignment retains every positive linked maximum peptide,
falling back to the union of supported chain maxima only without linked support.
Missing/invalid partners retain their valid marginal and explicit validity flags.
Paired marginal selection requires explicit human or mouse species, even without
calibration; ``--species any`` is rejected before input loading or search.

This opt-in exports linked ``candidates``, ``TRA.candidates``, ``TRB.candidates``
and ``calls`` tables with separate continuous axes and radius labels. It retains
no detailed hit table. Linked E-values and final selection E/P fields remain null;
``uncalibrated_paired_scope_selection`` identifies the final event. Marginal
fixed-predicate calibration is reported only on its own axis and never transferred
to linked or selected evidence. Final independent-pair empirical calibration
remains required; these outputs do not establish a calibrated paired detector.
The native build must provide ``Index.search_paired_batch_consume`` with
``skip_empty``. Compatible finite-cap controls use existing histograms; custom
control adapters and other geometries can retain materialized control hits.

.. code-block:: bash

   vdjmatch match cells.tsv --vdjdb reference.zip --paired --paired-marginals \
     --search-mode ball --max-radius 3 --matrix none --match-v --no-evalue --no-hit-table \
     --output-prefix results/paired-v-unit

Both fixed and graded single-chain modes accept ``--match-v`` and ``--match-j``
with an explicitly gene-labelled raw ``--control`` table or ``--no-evalue``.
Gene matching pools allele suffixes and preserves IMGT slash-containing names;
ambiguous genes are not reduced to the first call. Missing/ambiguous genes use
label-independent observed modes (bundled model modes when no calls exist),
with original calls and imputation flags retained for every query, including no-hits.
Supplied control species/locus and known gene locus select the population before
imputation. Distinct full junctions within each gene predicate define fixed control
sizes; abundance is ignored, and exact-hit filtering does not decrement model sizes.
One native histogram batch per gene-conditioned population supplies control counts.
Reference sizes follow the identical gene predicates. Sequence-only controls remain
rejected for this calibration. Paired hard-gene predicates require the explicit
``--paired-marginals`` opt-in; their linked and selected-event calibration remains
unavailable.
These fixed-ball E-values describe expected neighbour counts, not final selected-label
false-call probabilities; final-call calibration is a separate procedure.

The experimental fixed-K annotation mode and ``--top-k`` option have been
removed from the CLI and sample runner. The lower-level
``Annotator.ranked_candidates`` research API remains available for reproducing
historical diagnostics; it is not the production annotation workflow.

Experimental background-mass components
---------------------------------------

``match --search-mode unified --matrix none --unified-distance edit`` evaluates
an experimental five-edit ball (at most two insertions and two deletions).
``--unified-distance gapblock`` instead evaluates an anchored BLOSUM62
single-gap penalty ball at every query length. There is no length switch or
neighbour cap. These are component diagnostics, not a validated combined scorer.
The default fixed-mode scorer is unchanged.

Both components rank a pMHC by ``M/N * sum(1/(b(d)+1))`` over its distinct matched
junctions: ``M`` is the number of unique control junctions, ``N`` the pMHC
reference junction count, and ``b(d)`` the number of controls at distance at most
that neighbour's distance. The existing ``ned_score`` column carries this score;
``estimator=background-mass-v1`` distinguishes it from NED. Fixed-ball Poisson
statistics are separate model-conditional evidence, not posterior specificity.
The gap-block cutoff is five BLOSUM62 scale units, gap opening two units and
extension one unit, with allowed starts 3, 4, -4 and -3. Target and control use
the same predicate and exact-hit exclusion. Native dense distance buffers are
bounded to 64 MiB per query batch and reduced immediately.

These experimental components require unique sequence controls. Paired input,
V/J restrictions, matrix overrides and disabling controls are rejected. They do
not implement germline-loop calibration or establish improvements over existing
results. Use their labelled output tables to evaluate both geometries on the
same observations before developing a combined score.

Historical PSSM density baseline
--------------------------------

``vdjmatch historical-density SAMPLE --vdjdb REFERENCE --locus TRB
--epitope NLVPMVATV --mhc-a 'HLA-A*02:01' --threads 4 --output-prefix OUT``
recomputes the historical single-chain density component from raw controls.
Use ``--control TABLE`` for a supplied control (required for nonbundled loci).
It uses positional BLOSUM62 penalties within five substitutions, soft V-loop
weights, temperature 400 and the historical control-density floor of 0.01.
The original unit-distance significance calculation is reported separately.
Exact junction matches are excluded on both sides. Query identities remain
separate even when their sequences coincide. For reproduction only, the first
source-order V representative is used within each pMHC/junction group.

``--legacy-components --reference-population all-valid --germline-background RAW
--control RAW`` exports the original held-out calculation: NED with the raw
control CDF, separate V/J/length log odds and the union-level first-hit gate.
Omit peptide selectors to preserve the full reference union; an explicit MHC
selector is supported. This mode requires single-chain TRA or TRB input and
keeps fixed model populations. Empirical P-values, four-channel Sidak fusion
and cohort ROC summaries consume these flat outputs separately.

Benchmark exact exclusion is a hit-list filter only. Never remove a query
junction from the reference model, population counts, germline priors, motifs
or graph weights. Build each real or mock reference model once. Earlier
population-puncture experiments are withdrawn; their recorded results are
historical diagnostics and do not define benchmark eligibility.

``--mhc-match compatible --pool-reference`` pools intersecting family and allele
restrictions before counting unique reference junctions. Source restriction
labels remain intact in the reference. ``--targets-from-sample`` instead selects
each row's declared peptide, both MHC chains and class; it preserves observation
identity and reports unavailable targets explicitly.

The original BLOSUM kernel already represents likelihood evidence. Summed BLOSUM62
entries are log odds in half-bit units. With the actual quantised positional weights,
let A(q,r) be their weighted similarity sum. The native penalty is
P = 100 [A(q,q) + A(r,r) - 2 A(q,r)], so the existing exp(-P/400) equals
exp([A(q,r) - (A(q,q)+A(r,r))/2]/2). This restates the existing tempered,
self-normalised likelihood ratio without changing scores. Its reference aggregation
and empirical background correction stay intact. The radius-one enrichment P-value
remains a count test; optional final-call filtering uses ``calibrate-calls`` and
``apply-call-calibration`` with nulls from the matching complete scorer.

``historical-density --gzip-scores`` writes ``scores.tsv.gz`` using native gzip
level1; decompressed contents match the default plain ``scores.tsv``.

``--likelihood-score`` also exports ``likelihood_score`` while preserving every
historical score/status field. Single-chain output is log(density + 1e-6), plus
the same original germline term for sparse references. Paired output adds alpha
and beta log densities and the original combined germline term. It uses the same
reference/control searches, V weights and prior counts; it introduces no extra
sequence search. The paired log score is independent of the other queries in a
cohort. It is opt-in original-model evidence; its task regressions and missing
linked-pair evidence still prevent promotion.

``--all-targets --likelihood-score`` instead emits one original-model evidence
row per query and full pMHC in the selected species reference. Input epitope,
MHC and binder labels are ignored. Single TRA/TRB and paired inputs are supported;
paired scoring requires the same alpha control and raw germline background as
above. Missing or invalid queries and empty references retain explicit dispositions.
It uses the existing original components, including independent same-target
marginal references when linked pairs are unavailable. This is exhaustive evidence
for benchmark/calibration workflows, not an accepted winning-label call or a new
annotation default. The original kernel still excludes exact chain junctions.

``neighbours --all-targets --unbounded-nearest --nearest-rank
--background vdjdb-other --neighbour-weighting local-rank`` exports the matched
conditional-rank component for every full pMHC, with the usual radius, matrix,
V-loop and exact-exclusion options. Paired distances still refer to the same
linked reference pair. Nearest evidence is joined by query plus full pMHC, so
one target cannot supply another target's nearest rank.

Both commands accept repeated ``--withhold-epitope PEPTIDE`` arguments in
all-target or sample-declared target mode. The latter also supports null
normalization of compatible coarse-MHC hypothesis unions.
A withholding fold hides that peptide from candidate hits across every MHC
restriction. The full database remains the model: other target priors, complement
background populations, raw controls and reference geometry parameters stay fixed.
Thus each remaining candidate retains its original score. This differs from
physically deleting source records and refitting the model. Invoke human and mouse separately.
The existing source/reference/native recipes record all options and withheld
source-chain record counts; component P-values do not calibrate final selection.

For the prior-activation audit, ``--likelihood-score`` also exports
``likelihood_density_only`` and ``likelihood_full_prior``: the same log density
without the germline term, and with the complete original germline term at every
reference size. These diagnostic formulas have no 500-key switch and use identical
hits, controls and model counts. Paired cohort output also exports both columns
with ``_alpha`` and ``_beta`` suffixes for separate chain evaluation. They do not
replace the historical/default score or establish an optimal activation rule.

``--germline-background RAW_TABLE`` adds the original sparse-reference unpaired
formula: log(density + 1e-6) plus the V/J/length log likelihood ratio when the
reference has fewer than 500 unique junctions. Background multiplicity is retained
for these priors; control search counts unique productive junctions. These are
different populations. Alpha priors omit junction length. The manifest names
omitted components; a density-only run is not the complete historical scorer.

``--locus paired --targets-from-sample --germline-background RAW_TABLE
--alpha-control RAW_ALPHA_TABLE`` reproduces the original per-task cohort rank
sum of alpha density, beta density and the combined germline prior. ``--control``
selects the beta search background. For exact historical reproduction, the alpha
V/J prior uses the supplied beta prior table. This command uses independent
chain reference pools and excludes exact junctions per chain. It reports separate
chain P-values; it does not invent a combined P-value. Missing partners and
invalid inputs retain their original pair/chain IDs and explicit statuses.

``--gapped-extension`` is an experimental single-chain ranking candidate for pooled
or sample-declared targets. It preserves the original density and adds reference
edges outside the equal-length five-substitution ball, at every supported length.
Each additional edge contributes w(V) exp(-total distance / tau) divided by the
empirical CDR3 background mass, with the original 0.01 denominator floor. Here
w(V) is 1 for the same allele-stripped V gene and 0.25 for different genes. This
is a categorical gene prior: different genes retain the 0.25 factor even when
their native loop templates are identical. Total distance already includes
continuous V CDR1/CDR2/CDR2.5 differences; no additional loop-similarity multiplier
is applied. Background counts use CDR3 distance alone.
``--gap-radius`` defaults to 90 and ``--gap-temperature`` to 12. The output retains
``historical_density`` and, when priors are supplied, ``historical_score``.
Flat diagnostics report the pre-prior, pre-V-weight same-gene and cross-gene
kernel sums and accepted-edge counts. ``gapped_density`` equals the same-gene
sum plus 0.25 times the cross-gene sum. ``gapped_floor_density`` is the weighted
contribution of edges whose empirical denominator is strictly below 0.01. Best
accepted-edge total, weighted CDR3 and V-loop distances retain native distance units
units, with the first source-order reference chosen on total-distance ties.
Empty or unavailable extension outputs have zero sums/counts and null best-edge
distances; rows that never enter scoring keep null diagnostics and their status.
``p_enrichment`` remains the original component test; the combined ranking has
no calibrated P-value. This option does not change production annotation defaults.

``--gapped-extension --gap-geometry historical-pssm`` instead extends the original
full-junction positional BLOSUM Gram kernel. Integer positional weights use the
longer junction as a symmetric frame; equal-length penalties reproduce the
original kernel exactly. The sole gap start is Cys-relative column6: six
matched prefix residues, with the conserved Cys at index0. On junctions shorter
than six residues the native rule uses ``min(6, shorter length)``.
Equal-length pairs have no gap and ignore the gap prior. For length difference
``d > 0``, the gap charge is ``2800 + (d - 1) * 1400``; at ``d = 0`` it is zero.
The ball cutoff
is7000, and the exponential kernel scale remains400. The estimator is
``pssm-plus-historical-pssm-apex6-extension-v2``.
These fixed choices apply across the length range; there is no length switch or
hit cap. This mode retains original raw V-gene weights, including their historical
alias handling. Controls count the same weighted full-junction geometry, with
full exact identities excluded. Highest-contribution reference junction/V,
penalty, control count, contribution and gap length are exported. Block placement
is explicitly unavailable from the native scalar output. The original density,
score and P-value remain separate; this experimental addition has no combined
P-value. It requires positional gap-block support in the installed seqtree binary
and rejects older implementations before loading controls. Clipped-geometry
radius/temperature options do not apply.

``--gapped-extension --gap-geometry matched-pssm`` uses the same positional
geometry for every accepted nonexact reference neighbour, including original
substitution-ball edges within the weighted cutoff of 7000. It computes each contribution as
``V_weight * exp(-penalty / 400) / max(N * control_count(penalty) / M, .01)``.
``N`` is the number of unique reference junction representatives, ``M`` the
number of unique full control junctions, and ``control_count(penalty)`` counts
nonexact controls under that identical weighted predicate. The same legacy
expected-count floor is retained; this is not a finite-sample P-value.
The score is this matched-kernel sum, rather than original density plus an
extension. ``historical_density`` remains the unchanged original comparator;
``gapped_density`` reports the whole matched sum in this mode. The estimator is
``matched-positional-kernel-apex6-v1``. Original germline-prior fusion, if requested,
remains explicit. Neither this mode nor an additive extension assigns a new
confidence or changes the original component significance gate.

``--pssm-kernel-scale`` sets the positional decay scale (default 400); it must
be finite and positive. It changes kernel weighting, not the 7000 retrieval
cutoff, control counts, original comparator or component significance.
See :doc:`scoring` for the equations, statistical interpretation and limits.


Native junction/V-loop distances
--------------------------------

Experimental background enrichment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For single-chain declared-target comparisons, add ``--targets-from-sample
--background real --control RAW_TABLE`` to ``neighbours``. The same
control-table interface supports ``--background generative``; this flag declares
its origin and does not generate sequences. ``--background vdjdb-other`` builds
controls from other-peptide annotations in the selected species/locus reference,
without a separate control file. Other MHC restrictions are retained in that
background. A key annotated to both the target and another peptide remains in
the other-peptide population through its other annotation.

The control table must provide full ``junction_aa`` and V calls (legacy
junction columns are accepted only outside AIRR). Present species/locus columns
are filtered; without them the caller declares those properties. Reference and
real-control counts use distinct full-junction/resolved-V-allele keys, ignoring
abundance. Generated controls retain every eligible draw, including repeated
keys, to preserve their empirical generative mass. Ambiguous
V calls remain unavailable. ``--exclude-exact`` removes every full-junction
identity from counts and population sizes on both sides.

The additional ``OUT.enrichment.tsv`` exports ``n_reference``,
``reference_population``, ``n_control``, ``control_population``,
``reference_weight``, ``control_weight``, ``expected_count`` and
``enrichment=reference_weight/expected_count``, where
``expected_count=reference_population*(control_weight+0.5)/(control_population+1)``.
With the default ``--neighbour-weighting ball``, weights equal ball counts.
With ``--neighbour-weighting linear`` and a positive radius, each receptor
contributes ``max(0,1-distance/radius)``. Counts and weights are reduced in
one native traversal; boundary matches count but have zero weight.
The manifest states the kernel and smoothing contract; ``zero_control_weight``
distinguishes zero weighted support from ``zero_control_hits``.
Empty populations produce null enrichment with ``empty_population`` status;
zero control support is reported separately. This is a regularised ranking
statistic, with no P-value or confidence calibration.

Use ``--junction-ends trim3`` for symmetric end exclusion or
``--position-weighting significance`` for the complete empirical flank/core
profile. They cannot be combined. The profile is quantised to 0.01, with a
minimum 0.01; it is not a gene-specific NDN boundary model. The existing candidate
minimum distances/counts retain their unweighted comparator geometry. The
manifest distinguishes that output from the new enrichment geometry. It also
records control loading, count-batch wall time (including Python adapter
preparation), and total enrichment-stage wall time; these are separate from
the complete command wall time. Paired ``--background vdjdb-other`` is supported
with uniform positional weights. It uses the summed alpha/beta distance to each
linked pair, distinct full-pair/V-allele keys and exact exclusion only when both
full junctions match. The paired enrichment IDs map through ``OUT.queries.tsv``
to the original ``pair_id``. Real/generated paired backgrounds remain rejected;
independent chain products are not a joint null.


``--stem-depth`` adds the joint V--J hairpin stem term to the junction distance.
The CDR3 loop folds as a hairpin, so junction position ``i`` sits opposite position
``L-1-i``; a stem level is that pair, and a level counts as conserved only when
**both** of its residues are conserved, so the two residue costs combine as a
maximum rather than a sum. The term is evaluated on untrimmed junctions, which is
why it adds information the trimmed alignment cannot see, and it uses the same
substitution table selected by ``--distance-matrix`` as the junction alignment
and V loops. The clipped geometry uses bounded substitution costs; targeted annotation uses unclipped
BLOSUM62 Gram penalties.

The default ``--stem-depth 0`` omits the terminal stem term.
``--stem-depth 4`` is the configuration selected by the four-task benchmark
objective; ``--stem-weight`` defaults to ``1.0``, which spans the same range as one
aligned residue. Depth is a declared residue count, not a per-clonotype boundary.

For ``--locus paired`` the stem is applied to **both** chains and the two terms
are summed. The distance carries one stem per chain, using the same sum as
the other paired distance terms. Per-chain benchmark views from paired sources
do not establish true linked-pair accuracy or joint calibration.

======================  =========  ===============================================
Parameter               Default    Meaning
======================  =========  ===============================================
``--stem-depth``        ``0``      stem levels; 0 omits the stem, 4 is the retained depth
``--stem-weight``       ``1.0``    weight in shared substitution-matrix cost units
======================  =========  ===============================================

Both are rejected at the command line if negative.

The stem applies to dense distances and native single/linked-pair ball counts and
linear mass. Native counts evaluate factored junction-end codes alongside the
small V-loop table, avoiding a query-by-reference offset matrix. Full junction
identity still defines exact exclusion. Python and native capability are checked
before counting.

``--junction-ends trim3`` is an experimental comparison that removes three
residues at both junction ends. The default ``--junction-ends trim32``
preserves the clipped BLOSUM62 geometry three-N/two-C trim. Both keep full-junction exact
exclusion and V-loop distances; paired distances use the same option on both
linked chains. The manifest records the trim and distinguishes the experimental
geometry from the native default. No confidence calibration is added.

``vdjmatch neighbours SAMPLE --vdjdb REFERENCE --locus TRB --threads 4
--radius 90 --exclude-exact --output-prefix OUT`` exports minimum distance and
neighbour counts per pMHC, plus a separate query-availability table. It reproduces
the default clipped BLOSUM62 geometry distance: three times the trimmed, restricted-gap CDR3
distance, plus aligned V-gene CDR1, CDR2 and CDR2.5 distances. The bundled loop
models and their license/provenance are under ``resources/vloops``. Use
``--species human|mouse`` to choose the organism. Missing allele
suffixes explicitly resolve to ``*01``. Model-proven co-locus aliases resolve
uniquely without changing raw calls, including underscore spellings such as
``TRAV14_DV4*02`` for ``TRAV14/DV4*02``. The allele is preserved and must
exist in the selected organism/locus model; ambiguous calls remain unavailable. The
query table retains raw calls separately from resolved alleles, including chains
with missing junctions. Unknown alleles and junctions shorter
than eight residues remain unavailable in the query table. This command reports
distances, not specificity probabilities.

The default radius is ``22.5`` times the maximum effective substitution cost
of the shared clipped BLOSUM62 matrix (currently 4, giving radius 90).
``--radius-factor`` changes that factor; ``--radius`` supplies an absolute
historical distance instead. The CDR3 gap cost is three matrix units
(currently 12); CDR3 substitutions retain weight 3 and V-loop substitutions
weight 1. Distance bins use the positive matrix-cost greatest common divisor.
Common integer matrix rescaling therefore rescales gaps, V-loop distances,
bins and the factor-derived radius together. Fractional cutoffs round down
to native units (1/200); linear weights use that same effective cutoff.
The manifest records the requested and effective radius and matrix scale.
The Python shortlist audit also derives its default radius from the selected
matrix. Stem4 reference-call radii, kernel widths and manifest units use the
prepared geometry's scale rather than independently recorded numeric constants.

``--locus paired`` accepts linked AIRR rows with ``clone_id`` and declared
TRA/TRB roles. It sums alpha and beta distances against the same linked database
receptor; ``--radius`` applies to that sum. Exact exclusion removes a reference
only when both full junctions match. Original independent marginal minima must
be computed from separate full single-chain reference pools; they are a distinct
baseline. Paired output retains missing partners and original gene calls.

Unresolved V calls use the lexically first complete co-locus loop template
from the fixed species model. The template does not depend on queries, mocks
or target labels. Original calls remain in the output with
``v_geometry_imputed`` flags (one flag per chain for paired inputs).
Distinct unresolved source keys remain distinct when their geometry coincides
with a supplied template allele. Invalid junctions and missing paired partners
retain their existing unavailable dispositions.

``--targets-from-sample`` selects compatible declared restrictions and counts
unique junction/V search keys within their union. Query status distinguishes
missing references, unavailable sequence inputs, and searches without
neighbours. Explicit ``--epitope``, ``--mhc-a`` and ``--mhc-b`` selectors support
``--mhc-match compatible``; exact selection remains the default.

The existing seqtree gap-block batch kernel computes CDR3 distances; vectorized
lookups add germline-loop distances. Single-chain distance outputs are bounded
to 64 MiB per batch; paired matrix work uses a 32 MiB budget. These allocation
bounds do not establish whole-process RSS. The clipped BLOSUM62 matrix has nonnegative bounded substitution penalties; ordinary
seqtree BLOSUM Gram penalties are a different metric.

Exhaustive clipped BLOSUM62 geometry nearest reference
------------------------------------------------------

``neighbours --targets-from-sample --unbounded-nearest`` reports the
minimum native clipped BLOSUM62 geometry distance over every usable, nonexcluded reference in each
query's declared compatible pMHC union. A minimum may exceed ``--radius``;
``n_neighbours`` continues to count only references within that radius. With
``--exclude-exact``, exclusion compares the complete junction. In paired mode,
only equality of both complete junctions excludes a reference pair, and distance
is the alpha-plus-beta sum to the same linked reference pair.

The command emits one candidate per available query/declared target when any
nonexcluded reference remains. An exact-only excluded pool produces no candidate
and keeps ``no_neighbours`` disposition. Unknown genes, invalid junctions,
unselected species and missing partners retain their query rows and existing
statuses. Without this flag, the original radius-ball candidate output is
unchanged. The new flag requires ``--targets-from-sample``.

This mode reduces each existing bounded native distance-matrix batch to query
minima and radius counts; it does not materialize all query/reference hit pairs.
The manifest records the selection mode and all CLI parameters. The resulting
nearest distances are uncalibrated distances, suitable for separately defined
comparator ranking; the command does not assign a specificity probability.

Shared original calculations across cohorts
-------------------------------------------

``historical-density --locus paired --targets-from-sample --cohort-column cohort``
shares the existing per-target chain searches across observations while computing
paired rank fusion separately for each explicit cohort. The raw column must be
nonempty, agree across linked rows, and avoid canonical input and score names.
The option also emits original single-chain scores and separate chain statuses;
isolated chain rows retain their identities and stay outside paired ranks.
Without the option, existing output columns and rank cohorts are unchanged.
The supplied prior remains the original beta background for both chains.
No combined paired P-value or new confidence interpretation is introduced.

For commands with an explicit thread budget, the CLI caps the Polars thread pool
at that budget and honors a smaller parent ``POLARS_MAX_THREADS`` allocation.
The native sequence engine still receives the requested ``--threads`` value.
Malformed parent thread allocations fail before loading inputs.

vdjmatch apply-call-calibration
-------------------------------

The `mouse TRA operating-point recommendation <https://huggingface.co/datasets/isalgo/airr_benchmark/blob/72fb387cf6098f4135732c4e7693ea951698cb2a/vdjmatch/releases/1.0.0/calibration/README.md>`_
uses ``kernel_score >= 0.17`` for the recorded historical all-valid,
exact-retained stem4 event. It calibrates a nominal 1% null-query any-call rate;
this is not a per-panel release gate or false discovery rate among called receptors.
The four same-source observed subsamples contribute 189 calls among 19,045
panel-query observations (0.9924%); they are calibration inputs. The separate
generated panel has 872 calls among 100,000 queries (0.8720%). Correct accepted
reference calls remain 6 among 6,791 queries, unchanged from the previous cutoff.
Complete per-panel counts and the recall limitation accompany the recipe.
Application requires its matching score-producer manifest; this cutoff does not
transfer to curated references or another database build.

Apply a previously computed cutoff to independent query-level scores:

.. code-block:: bash

   vdjmatch apply-call-calibration queries.tsv --calibration calibration.calibration.json --output calls.tsv

The input has one unique ``query_id`` and the recipe's numeric score column per query. The command
does not read query outcome labels or fit a threshold. It preserves all input
columns and adds ``call_status`` (``accepted`` or ``no_supported_match``) and
``call_score`` (null for rejected queries). A missing score is rejected; a tie at
the frozen cutoff is accepted. An ``abstain_all`` calibration rejects every query.
The manifest records both input hashes and exact accepted/rejected/no-score counts.
Input label columns, if present, remain metadata or candidate evidence; they are
not accepted predictions when ``call_status`` is ``no_supported_match``.

Add ``--null-cases calibration-input.tsv`` to report per-query ``p_null`` and
``E_null`` from the frozen null scores. This must be the original combined
``calibrate-calls`` input: its SHA-256 must match the saved recipe. Only rows
with ``outcome=null`` define tails. For each source, let N be all null queries
(including no-hits), and k the number whose selected scores are at least as
extreme as the query score (including ties). Then P=(k+1)/(N+1) and E=k/N.
The output uses the largest value across sources and exports every source's k,
N, P and E to ``calls.tsv.null-tails.tsv``. An unsupported query has P=1 and E=0.
Input scores and the frozen cutoff decision remain unchanged.

``calibrate-calls --score-column likelihood_score`` freezes that column in the
recipe. Application and null tails then require the same column, retaining any
historical ``score`` column as evidence. The default remains ``score``.

P is a conservative rank P-value under exchangeability of query and null scores
from the same complete selection pipeline. E is the observed expected false
calls per query at that fixed score threshold; it is not a confidence bound.
These tails do not establish biological binding probability or independence of
peptide-withheld observations. Query outcome labels do not enter this calculation.

Use a calibration produced by the same entire scorer, reference, geometry and
winning-label selection pipeline. A recipe containing ``pipeline_binding``
requires ``--input-manifest producer.manifest.json``. The command validates its
single-chain species, reference population/hash, exact-hit policy, geometry,
selection and model/software identities before writing outputs. Different or
missing identities fail explicitly. Generic unbound empirical recipes remain
supported; their caller owns pipeline compatibility. Application supplies the
rejection step and preserves the annotation API's existing decisions. Matched
open-set calibration still requires fresh null scores from that exact pipeline.
Retained exact-reference evidence is not proof of
biological binding. Rejection means insufficient support for a VDJdb assignment,
not proof that a receptor cannot bind anything.

Open-set performance depends on repertoire composition. For illustration, among
10,000 queries with 99% unknown specificity, 1% false calls among the 9,900 unknown
queries give 99 false calls. Even perfect recall of the 100 covered queries then
gives only 100/199 accepted-call precision. The empirical null-call rate and the
biological binder probability are different quantities.

vdjmatch calibrate-calls
------------------------

Reduce final query-level decisions into tie-aware empirical operating points:

.. code-block:: bash

   vdjmatch calibrate-calls cases.tsv --target-fpr 0.01 --output-prefix calibration

Add ``--input-manifest producer.manifest.json`` to freeze the complete supported
single-chain kernel pipeline identity in the recipe. Apply that recipe with
``apply-call-calibration --input-manifest query-producer.manifest.json``; the
producer formats may differ, but their canonical scoring identities must match.
This keeps exact-retained ordinary calibration separate from exact-excluded
benchmark calibration. Paired and mixed-population bindings require a separate
declared decision and null contract.

The flat input has one unique ``query_id`` per query, a numeric ``score``, an
``outcome`` (``correct``, ``incorrect`` or ``null``), and ``null_source`` for
null queries. Correct/incorrect describe the winning label on labelled VDJdb
queries. Resolve candidate competition and ambiguity before making this input;
per-candidate rows would estimate a different error rate. Use null scores for
abstention, retaining those queries in denominators. Invalid queries must have
an explicit disposition in the upstream source record.

Higher scores are stronger by default; ``--lower-is-better`` supports distances
or P-values. Ties are indivisible. The command selects the most permissive
cutoff for which each null source independently meets the requested empirical
any-call rate. If none is usable, it reports ``abstain_all`` with a null cutoff.
It writes ``.operating_points.tsv`` and ``.calibration.json`` with separate
null counts/rates, labelled-query correct-call recall and call precision.
Generated nulls do not enter the labelled-query precision denominator.
Labelled queries are optional. A null-only input fits the same cutoff and
source rates, with recall and precision left unavailable.

This command computes an empirical recommendation candidate. The input hash
is recorded; callers must also pin the producing scorer, geometry, reference,
model/seed, identity exclusion and query populations before recommending it.
Calibration data must match the actual winning-label sample annotation path.
A target of 1% measured null calls is not a guarantee of 1% biological error.
No reference-specific recommended numeric cutoff is shipped yet.

Recommended P-value/probability cutoffs are optional sample-query guidance.
Principal benchmark evaluations use all eligible values: ROC/AP use continuous
scores, and accuracy/F1/confusion state their decision rule. Optional metrics at
a recommendation cutoff remain separate and retain below-threshold/no-hit rows.
The command does not change annotation defaults or enforce cutoffs on benchmarks.

Single-chain ``match --exclude-exact`` filters exact target/control hits while
retaining the full target and control population sizes in the evidence reducer.
The historical ``fixed_ball_punctured`` output label is retained for compatibility;
it denotes hit filtering, not removal from the model.

Exact-self exclusion is a benchmark option. Ordinary sample annotation retains
exact VDJdb hits; recommendations for nonexact extrapolated matches should not
screen out that direct evidence. Exact sequence identity is certain. A shared
single-chain sequence with multiple specificity labels still requires an
ambiguous disposition rather than a unique antigen call from identity alone.

vdjmatch reference-calls
------------------------

This benchmark command queries distinct full-junction/V/J keys from peptides with
at least 30 valid records for the selected species and chain (before self-exclusion), retains every query disposition, and excludes full-junction self hits.
Ordinary ``match`` continues to retain exact hits. ``--leave-epitope-out`` instead
queries each distinct key/truth-peptide appearance and removes that entire peptide
across all MHC restrictions before candidate reduction. One native query batch
covers all eligible folds; the result agrees with physically removing each peptide.
The searchable reference remains complete, as in ordinary sample annotation.
``--min-target-records`` changes benchmark target eligibility only; target support
and retained/excluded record counts are exported.

.. code-block:: bash

   vdjmatch reference-calls --vdjdb reference.zip --locus TRB \
     --threads 12 --output-prefix labelled
   vdjmatch reference-calls --vdjdb reference.zip --locus TRB \
     --leave-epitope-out --threads 12 --output-prefix withheld

The default ``--scorer original`` retains the existing scoring and curated
reference-valid policy. ``--reference-population all-valid`` explicitly includes
all canonical valid junctions, including source records excluded by benchmark
curation; ``curated-benchmark`` remains this command's default.

``--scorer stem4`` exports the final selected ``kernel_score`` for single-chain
calibration. It uses raw BLOSUM62 costs, trim3, depth 4 stem geometry, radius 585,
the recorded 153-cutoff ladder and V-bounded evidence. Each peptide first takes
its maximum over MHC restrictions, then all peptides compete; ties retain an
explicit ambiguous disposition. The entire selected-policy reference includes
sparse peptides. Peptide-withheld folds remove full peptide memberships before
counts, population sizes and shared-key exclusivity, in one native batch.
Distinct anchored-AA-junction/V/J identities are retained, with truth peptide
also part of withheld query identity. Exact exclusion removes the full junction
across V/J identities. Invalid genes and no-support queries remain in outputs.

``--exact-policy exclude`` is the benchmark default. For a separately calibrated
ordinary annotation pipeline, ``--exact-policy retain`` keeps exact evidence
while preserving peptide withholding and all query identities. The manifest
records the policy; a cutoff from excluded queries cannot be applied to retained
queries. This option produces calibration inputs and does not supply a cutoff.

.. code-block:: bash

   vdjmatch reference-calls --vdjdb reference.zip --locus TRA --scorer stem4 \
     --reference-population all-valid --leave-epitope-out --threads 12 \
     --output-prefix stem4-withheld

Stem4 uses the other-peptide reference background. Explicit ``--control`` and
nondefault ``--scope`` are rejected; ``human``/``mouse`` species aliases are accepted.
Only ``--calibration-score score`` is compatible; component E/P values cannot
calibrate this final scalar. Its outputs are calls, calibration cases, truth,
source-record target support and a geometry/model/source-hashed manifest.
Calibration must bind this scorer, selection and reference policy.
Unresolved V calls use the same fixed co-locus template as ``neighbours``;
calls retain original V identities and geometry-imputation flags.

Paired reference nulls require both ``--locus paired`` and
``--leave-epitope-out`` with ``--scorer stem4``. Queries retain genuine TRA/TRB
linkage, both V/J identities and their source complex IDs in the truth table.
Eligibility requires at least 30 source records for either chain and peptide
before exclusion; paired or unique-key counts do not replace source counts.
The truth peptide is removed from all linked and marginal hit populations,
counts, sizes and shared-label memberships. Unlinked compatible chain references
remain available to the two marginal axes.

.. code-block:: bash

   vdjmatch reference-calls --vdjdb reference.zip --locus paired --scorer stem4 \
     --leave-epitope-out --reference-population all-valid --exact-policy retain \
     --threads 12 --output-prefix paired-withheld

The paired statistic is ``kernel_any_axis_score``, the maximum of the final
selected linked, alpha and beta scores. This defines an any-axis call event and
introduces no combined label. Paired ``annotate`` outputs preserve the linked
``kernel_score`` and existing labels, and add ``kernel_alpha_score``,
``kernel_beta_score`` and ``kernel_any_axis_score`` with marginal ties.
Fit/apply a matching recipe with ``--score-column kernel_any_axis_score``, or
use a flat ``score`` column projected from that statistic with its matching
``--input-manifest``. A paired recipe rejects the linked ``kernel_score``.
Reference population, exact policy, geometry, sources and model must match;
the command provides null inputs and does not supply a precomputed cutoff.

``paired-score`` exports a separate continuous linked-first candidate from a
native paired selection, with no new search or fit:

.. code-block:: bash

   vdjmatch paired-score native.selection.tsv --input-manifest native.manifest.json \
     --output-prefix linked-first

Positive linked evidence supplies ``score``; otherwise the existing any-axis
score retains marginal support. No supported axis gives an explicit no-score
row. All original identities, component values and metadata remain unchanged.
The reducer validates the native any-axis maximum and records input/output
hashes. Its ``.scores.tsv`` is continuous ranking evidence, unlike the decision
indicators from ``paired-labels``. It changes no default or label assignment;
the existing any-axis call calibration does not apply to this candidate.

The original scorer requires a matching raw ``--control`` repertoire for nonbundled
species/loci. Its output
includes continuous candidates, calls, full truth labels, calibration input and a
manifest. Wrong labels and no-hit/ambiguous dispositions remain in denominators.
Calibration input uses the existing winning-label ranking ``score`` by default.
``--calibration-score E`` exports the expected number of background hits in the
candidate reference population; ``--calibration-score p_enrichment`` retains
the enrichment-tail diagnostic. Neither diagnostic measures winning-label error.
The enrichment P-value depends on the observed
number of neighbours and is not the probability of a false winning-label call;
it is not interchangeable with an E-value or a biological posterior.
``--calibration-score score`` exports the existing winning-label ranking score.
``calibrate-calls`` reports its empirical expected false calls per query as the
largest observed null-source any-call tail at each cutoff. For a sample of Q
queries evaluated at a common cutoff, multiply this value by Q to obtain the
empirical expected number of false calls. This expectation requires no independence
assumption, but transport to a new population remains unvalidated. It is distinct
from the Poisson neighbour-count E-value and is not a confidence bound.
Neither this command nor ``calibrate-calls`` supplies a biological posterior or
filters principal accuracy/ranking evaluations with a recommended cutoff.

Fixed V/J versus NDN geometry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The experimental ``neighbours --position-weighting ndn`` arm requires
``vdjmatch[ndn]`` and a seqtree native build with per-input positional weights.
Use it with ``--background`` and default ``--junction-ends trim32``. It weights
verified ARDA V/J junction flanks by0.1 and NDN by1, retaining the clipped
BLOSUM62 substitution matrix, V-loop contribution and positive length-gap costs.
Aligned substitutions use the larger of the query and reference position weights.

Annotation uses one shared query/reference/control catalog. Conflicting or
ambiguous J calls, missing calls, changed junctions or unusable boundaries retain
uniform full weight. No benchmark sequence is repaired. The CLI emits
``.boundaries.tsv`` with every eligible junction/V key and its disposition; the
manifest records geometry, counts and source identities. Generated draw
multiplicity and observed junction/V population units remain unchanged.

This option changes background enrichment/local-rank geometry only. The nearest
distance output retains its clipped BLOSUM62 geometry comparator geometry. It provides an
experimental ranking contrast, with no default annotation or confidence change.

Fixed V-loop information contrast
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``neighbours --v-loop-weight 0`` is an experimental CDR3-only
contrast against the default ``--v-loop-weight 1``. Both retain identical
query/reference/control V-call eligibility, full identity keys and optional
exact puncturing. Only CDR1/CDR2/CDR2.5 distance offsets are removed;
CDR3 substitution, gap, trim and positional geometry remain the same.
The option applies to nearest distances and background count/mass/rank
batches, including linked alpha/beta pairs against the same reference pair.
No arbitrary fitted weight or task-specific blend is exposed.

This contrast isolates V-loop information within the eligible population;
it does not establish a recall-preserving CDR3 shortlist, recover queries
with missing V calls, or transfer calibrated recommendations between
geometries. The manifest records the experimental geometry and weight.
Native rival manuscript parameters remain unchanged.

Reference-selected distance radii
---------------------------------

``reference-radius`` derives a distance radius from nonexact nearest neighbours
within each reference full-pMHC group::

   vdjmatch reference-radius --vdjdb REFERENCE --locus paired --threads 4 \
      --output-prefix reference-radii

The strongest confirmed normalized Kneedle bend wins. Outputs ``.radii.tsv``,
``.nearest.tsv`` and a provenance manifest. Distinct full
junction/resolved-V keys define the search population; paired mode uses genuine
linked receptors and the summed distance to the same reference pair. Source
observation counts and unavailable keys remain explicit.
Unresolved V calls receive the fixed co-locus template used by ``neighbours``.
The nearest-distance table retains ``original_v0`` (and ``original_v1`` for
pairs) and geometry-imputation flags; imputation does not merge original keys.

Use ``--matrix blosum62`` for full BLOSUM62 costs; the default remains
``clipped``. The selected matrix supplies junction and V-loop costs, gaps,
coverage units and the knee-radius ratio. The manifest records this geometry;
radii from different matrices are not interchangeable. ``pam250``, ``mj`` and
``vdjam`` are also supported by the existing distance backend.

``neighbours --reference-radius`` applies this selection to each
sample-declared compatible reference union. It requires uniform geometry,
V loops, ``--nearest-rank --unbounded-nearest --targets-from-sample``,
``--neighbour-weighting local-rank --background vdjdb-other``. Its additional
``.radii.tsv`` records the raw knee, selected integer radius and disposition;
``selected_radius`` accompanies each enrichment row. The exhaustive nearest
search remains uncensored.

Both commands use ``--knee-floor-factor 0.5`` times the median positive knee
within the same species, locus and geometry. If there are no positive knees,
the floor uses the median positive nonexact nearest distance. Missing knees
remain missing in ``reference-radius``; the scoring command explicitly falls
back to its positive ``--radius``. No positive scale means no inferred radius.
Rounding upward to an integer preserves the selected distance boundary.
These are experimental ranking radii, without transferred probability calibration.

Cohort-relative paired evidence ordering
----------------------------------------

The opt-in ``paired-evidence-order`` command reduces existing fixed-radius
clipped BLOSUM62 geometry evidence; it performs no sequence search and changes no default scorer::

   vdjmatch paired-evidence-order evidence.tsv --output ordered.tsv --threads 1

Supply one row per linked query, cohort, species and full pMHC. Required columns:

* ``query_id``, ``cohort``, ``species``, ``epitope``, ``mhc_a``, ``mhc_b``,
  ``mhc_class`` identify the observation and ranking population. Unspecified MHC fields may be null; query, cohort, species and epitope must
  be non-null and nonempty. Duplicate identities within a population are rejected.
* ``paired_score`` is the fixed-radius paired rank Bayes-factor statistic at
  distance 90; ``fallback_score`` is the original paired score.
  ``alpha_score`` and ``beta_score`` are the original single-chain scores.
  Higher values mean stronger evidence. These fields must be finite numbers or
  null; null scores tie below finite values. Convert unavailable infinite scores
  to null before invoking the command.
* ``n_reference_pair``, ``n_reference_alpha``, ``n_reference_beta`` are
  nonnegative UInt32 integer neighbour counts from the same fixed-radius geometry.
  Pair support counts common linked reference pairs, not independent chain hits.
  Null counts mean unavailable support and are treated as zero for routing.
* ``paired_available`` is a non-null Boolean indicating that both query chains
  are present and valid. A missing chain makes it false; an observed chain with
  no reference hits remains available even if its score is null.

Block 1 requires common-pair support and orders by ``paired_score``. Block 2
requires no common-pair support and at least one supported chain; it sums the
average ranks of alpha and beta scores computed separately within block 2.
Block 3 orders the remaining queries by ``fallback_score``. Unavailable pairs
always enter block 3. Equal statistics and equal rank sums remain tied.

Input fields and row order are preserved, with counts normalized to nullable
UInt32. Output adds ``evidence_block``
(UInt8, 1/2/3), ``evidence_status`` (``common_reference_pair``,
``single_chain_only``, ``no_reference_neighbour`` or
``paired_evidence_unavailable``), ``within_block_rank`` and ``evidence_score``.
Within each cohort/species/full-pMHC block, average ranks run from 1 to the
block size n. The score is ``(3 - evidence_block) + within_block_rank / (n + 1)``;
higher is better, and blocks have strictly disjoint score ranges.

Provide the complete cohort in one invocation: changing its membership can
change ranks and the block-2 order. Scores are not comparable across batches or
pMHC populations. Labels are never consumed, and no P-value, E-value, posterior
or biological specificity probability transfers to this empirical ordering.

Curated independent source support
----------------------------------

``reference-support`` exports distinct nonempty source-reference counts for exact
full junction/V/J, locus, species and full pMHC. Curated excluded records never
count. Repeated observations from one source count once. Distinct reference IDs
do not establish independent experimental validation or guaranteed binding.

.. code-block:: bash

   vdjmatch reference-support --vdjdb REFERENCE --min-refs 2 --output support.tsv

``--gene-level`` explicitly strips V/J allele suffixes for gene-resolution
comparisons. The default retains source allele identities. Missing receptor
fields do not establish replicated support. The output JSON records reference
identity, source/valid record counts, retained keys, table hash and elapsed time.
This is an evaluation subset tool; it never filters the searchable reference.
V-conditioned motif evidence
----------------------------

``vdjmatch motif-evidence SAMPLE --vdjdb REFERENCE --locus TRA|TRB
--species human|mouse --output-prefix RESULT`` computes experimental masked-motif
MF×IKF evidence for the sample's declared full pMHC. It uses the union of compatible
reference restrictions for target motif counts and all valid same-species/locus
references for global motif frequency. Full junctions, all four wildcard positions,
IMGT amino-acid classes and the original 2/3/5 motif support gates are retained.
``--exclude-exact`` is benchmark-only and filters exact junction hit contributions.
Motif support gates and global frequencies remain fixed. Repeated query identities remain separate.
Zero supported motifs, invalid queries and absent references have explicit statuses.
This command supplies a candidate component, not a switched scorer, joint paired
likelihood, P-value or specificity probability. Pairing remains in the output IDs;
run both loci and join by the original linkage to compare genuine paired queries.
``--threads`` sets a positive table-thread budget (default ``1``); BLAS stays
serial. Output tables and the manifest must not overlap inputs or lie inside
the reference directory.

Global evidence combination
---------------------------

Both commands accept a positive ``--threads`` table budget (default ``1``),
respect a smaller parent table cap, and keep numerical BLAS work serial.

``fit-evidence TABLE --model MODEL`` fits one experimental logistic formula from
``likelihood_score``, ``conditional_rank_score``, ``reference_coverage``,
``knee_shape``, ``knee_radius_ratio``, ``original_support`` and ``paired`` columns. Labels default to
``label``; ``--weight-column`` supplies explicit population weights. Shared
nonnegative evidence coefficients preserve monotonicity at fixed structure.
Missing evidence receives an indicator; unsupported training rows require an
explicit disposition before fitting. Model JSON records coefficients, scaling,
training groups and input hash.

The two scores are the preserved original log score and log conditional rank.
For reference coverage ``c``, normalized bend ``k``, radius ratio ``e`` and
original chain support ``u``, the shared sparsity gate is
``s = 1-u*c*(1-(1-k)*e)``. Single-chain support is zero or one; paired support
is the fraction of chains with positive original density. The GLM fits four
nonnegative dense/sparse score coefficients, plus mode and missing indicators.
Standalone masked-motif motifs remain a separate comparator.

``combine-evidence TABLE --model MODEL --output SCORES`` preserves every row and
exports a linear ranking margin with an unsupported status when all components
are absent. ``--require-heldout`` rejects overlap with training epitope groups
for benchmark evaluation. This experimental margin has no transferred component
P-value and is not a calibrated probability of binding.

``combine-evidence TABLE --fixed --output SCORES`` applies the unfitted formula
``likelihood_score+s*log(conditional_rank_score)``. If only one component
is available it uses that component; if neither is available it retains an
unsupported row. Fitted applications also export this formula as ``fixed_score``
for matched comparisons. Both require calibration of the complete selected score.

``combine-evidence TABLE --null-evidence NULL_TABLE --output SCORES`` puts
the original and conditional-rank components on matched empirical null-tail
scales before choosing the stronger component. Both tables require ``query_id``,
``species``, ``locus``, all four pMHC fields, ``likelihood_score`` and
``conditional_rank_score``. Null rows additionally require ``outcome=null`` and
``null_source``. Query IDs must be unique within each species/locus/full-pMHC
population. Query outcome labels are not read.

For each component and null source, its tail is ``(k+1)/(N+1)``, where ``k``
counts null scores at least as large as the query score, including exact ties,
and ``N`` includes every null query, including no-hits. The most conservative
source tail is retained. Output ``normalized_score`` (also ``combined_score``)
is ``-log(min(p_component_likelihood,p_component_rank))``. Missing components
have tail one; two missing components remain unsupported. This avoids blending
raw component scales and uses no fitted class coefficients or epitope exceptions.
Empirical steps can introduce ties; accuracy improvement is not assumed.

Add ``--structure-conditioned`` to use the reference-only gate
``s=1-reference_coverage*(1-(1-knee_shape)*knee_radius_ratio)``, with
``normalized_score = max(-log(p_L), -log(p_R)+log(s))``. Query-local
``original_support`` is not required and does not switch routes.
Here ``p_L`` and ``p_R`` are the component null tails above. This retains the
original normalized evidence and discounts rank evidence for dense structure;
``s=0`` disables the alternative route. A missing original score uses the
available rank tail without discount. This option also requires the four
structure/support columns used by the fixed formula. It uses no query labels,
fit or epitope-name switch; it is experimental, with accuracy and leave-query-out
stability requiring validation.

This maximum is a ranking statistic, not a combined P-value. Freeze the matching
scorer/reference/null recipe, apply the same complete selection to null queries,
reduce all candidate labels with the same tie/ambiguity rule, then calibrate that
selected statistic with ``calibrate-calls``. Both input hashes and exact grouping
keys are recorded. Earlier prediction tables are not inputs to a fresh benchmark.

Three-channel normalization (experimental)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``combine-evidence TABLE --null-evidence NULL_TABLE --three-channel --output SCORES``
normalizes ``original_score``, ``gapblock_score`` and ``stem4_score`` with the
same conservative null-tail rule. Both tables use the population keys above;
all three components must have the same reference, exclusion and null-query
identities. The evidence table additionally supplies ``original_weight``,
``gapblock_weight`` and ``stem4_weight``. Each row's weights must be finite,
nonnegative and sum to one.

The output retains inputs and weights, adds ``p_component_original``,
``p_component_gapblock`` and ``p_component_stem4``, each channel's
``log_evidence_*`` in nats and its weighted ``contribution_*``, and reports
``combined_score = sum(weight * -log(component_tail))``. Higher scores rank
stronger evidence. A missing component has tail one and contributes zero;
its weight is retained. If every positively weighted component is missing,
``combined_status`` is ``unsupported`` and the score is unavailable.
When ``original_supported`` is declared, a finite zero-density prior still
contributes to continuous ``combined_score``. The separate
``combined_call_score`` requires supported original evidence or a positive
kernel with positive weight; ``combined_status`` describes this call support.
Use ``combined_score`` for continuous ranking diagnostics, retaining all rows.
A finite rank does not authorize a recognition call.

This command applies explicit weights; automatic weights from graph, length and
V usage remain under investigation. It fits no coefficients and does not read
assay labels. The component tails describe their declared null populations;
the sum requires its own complete selected-label and paired-axis calibration
before reporting a match P-value.

``reference-radius --targets-from-sample SAMPLE`` profiles the compatible pMHC
reference unions declared by a sample. Its output also includes neighbour
coverage and normalized convex bend. ``--coverage-radius-factor`` defaults to
22.5 times the shared substitution-matrix scale; the bend uses unit-square
normalization. These describe the reference without changing search radii.
The radius ratio is regularized epsilon divided by epsilon plus the
matrix-relative coverage radius; no knee is encoded as zero with its status retained.

Legacy matrix and retrieval reproduction
----------------------------------------

``scoring-audit --task matrix --pmhc RAW --locus TRB --output-prefix OUT``
rebuilds the accepted positional substitution tables and their off-diagonal
Pearson correlations from an explicit raw VDJdb export. ``--task loo`` rebuilds
the legacy epitope-withheld matrix training and candidate retrieval tables.
Both accept ``--threads`` and record source/model hashes and options.
TRA and TRB are separate invocations.

These commands preserve the legacy unique-junction eligibility and, for LOO,
composition caps, first-label assignment and tied-score balanced trapezoidal
PR area. This metric is distinct from average precision and the current
full-VDJdb source-record eligibility rule. Supply ``--structural-inc`` to
include the structural matrix; without that model its row is omitted.

vdjmatch fixed-scope
----------------------------

Compare naive unit-edit search against declared peptide–MHC targets:

.. code-block:: sh

   vdjmatch fixed-scope assays.tsv --vdjdb reference.zip --scope 1,0,0,1 \
     --exclude-exact --reference-population curated-benchmark --threads 4 \
     --output-prefix fixed

The sample contains junction/V/J/locus fields and the assayed ``epitope``,
``mhc_a``, ``mhc_b`` and ``mhc_class`` columns. Compatible reference restrictions
are pooled for each declared target. ``--paired`` uses supplied ``pair_id``
links and retains incomplete pairs with their available chain evidence.

The four scope entries bound substitutions, insertions, deletions and total
operations. The default admits one substitution and no indels. These are
operation budgets; a rule based only on unrestricted Levenshtein distance and
net length difference can admit a different set of matches.

``.scores.tsv`` preserves every case. Nearest negative unit-edit distance and
neighbor counts are separate scores. Marginal counts distinguish unique full
junctions, joint junction/V/J clonotypes and source observations; linked counts
distinguish full junction pairs, paired clonotypes and complexes. Linked and
TRA/TRB marginal rankings remain separate. ``any_hit`` is their support union;
missing linked references cannot veto marginal support. No-hit counts are zero
and nearest distances are unavailable. No V/J similarity, abundance, background
or fitted probability weights affect these naive scores.

PSSM-density remains available through ``historical-density`` with its declared
model assets. It is a separate scoring option from fixed-scope distance/counts,
NED-v1, kernel scores and reference-weighted combinations. Comparisons must
preserve the same observations, target definitions and exact-exclusion policy.

vdjmatch trained-score
----------------------------

Score declared peptide--MHC targets with the bundled frozen alignment model:

This command requires native ``seqtree.Index.alignment_max_batch`` and histogram
aggregate-indel support. Older native builds fail explicitly before scoring.

.. code-block:: sh

   vdjmatch trained-score assays.tsv --vdjdb reference.zip --scope 4,2,2,6 \
     --max-indels 2 --exclude-exact --reference-population curated-benchmark \
     --threads 4 --output-prefix trained

Inputs use the same source-preserving boundary as ``fixed-scope``. The four
scope entries retain their substitution/insertion/deletion/total meaning;
``--max-indels`` additionally bounds the sum of insertions and deletions.
The defaults therefore admit four substitutions and two aggregate indels.
This candidate predicate is separate from the subsequent alignment optimization.
Matching uses full junctions including conserved anchors; bare AIRR ``cdr3_aa``
does not supply those anchors.

For each candidate, global Needleman--Wunsch maximizes raw VDJAM similarity
with linear gap charge :math:`\lambda=\max_{a,b}|M(a,b)|`. Let :math:`U` be
the selected alignment's residue-pair similarity sum without gap charges,
:math:`g` its gap-residue count, and :math:`H(q),H(r)` the full-junction diagonal
self sums. The CDR3 term is :math:`C_3=U-\max(H(q),H(r))`. The frozen model is

.. math::

   \eta=b+\beta_V S_V+\beta_1 S_1+\beta_2 S_2+\beta_J S_J
          +\beta_3 C_3+\gamma g,
   \qquad P_{\mathrm{fit}}=-\operatorname{expm1}(-\exp\eta).

The bundled V/CDR1/CDR2/J terms :math:`S_V,S_1,S_2,S_J` are already
self-normalized. The fitted coefficient :math:`\gamma` is applied after
alignment selection; it does not replace :math:`\lambda` in the NW search.
All coefficients and similarities retain their bundled precision. No model
file or fitting step is supplied by the user.

The target ranking is the best hit's :math:`P_{\mathrm{fit}} I`, with
:math:`I=-\log_{10}((k+1)/N)`. Within each species/chain, :math:`N` counts
distinct full-junction/epitope roles and :math:`k` counts distinct junctions
carrying another epitope within the same candidate scope. Cross-reactive
junctions retain every epitope role. Benchmark exclusion deletes the complete
query junction from every role and adjusts both counts before scoring.

Unknown or ambiguous V/J calls remain explicit unscored rows; the command
does not invent an allele or use an unspecified missing-gene fallback.
Rows with missing or other-species metadata remain ``unselected_species``;
``--species`` selects the population scored in this invocation.
``--paired`` retains source pairing and both trained marginal rows, including
incomplete pairs. A trained joint probability is unavailable: no paired model
is fabricated from the single-chain coefficients.

Outputs are ``.scores.tsv`` and flat ``.provenance.tsv``, retaining source
identities, no-hits, model/native hashes and the exact budgets. These are
ranking/model scores, not calibrated match P-values or binding posteriors.
Calibration from another scorer cannot be transferred to them.

vdjmatch reference-structure
----------------------------

Inspect the reference composition that can guide chainwise scoring weights:

.. code-block:: sh

   vdjmatch reference-structure --vdjdb reference.zip --species human --locus TRB \
     --epitope NLVPMVATV --threads 4 --output-prefix structure

The command writes ``structure.structures.tsv`` with one row per full peptide–MHC
population, ``structure.degrees.tsv`` with its degree histogram, and
``structure.lengths.tsv`` with coverage at each exact junction length, and
``structure.genes.tsv`` with V/J usage overall and by exact length.
``--targets-from-sample AIRR_TABLE`` uses the compatible restriction unions declared
by that table. ``--reference-population`` selects ``curated-benchmark`` (default)
or ``all-valid``.

Add ``--control RAW_REPERTOIRE`` to compare V/J composition against the supplied
rearrangement control. The existing control reader and event normalizer preserve
joint junction/V/J recombination-event weights and report excluded raw rows.
Read abundance does not define event counts. The additional
``structure.control-genes.tsv`` reports event counts and unique junction/gene-key
counts, overall and within each exact junction length, with alleles pooled.
``structure.genes.tsv`` retains every reference row and denominator, and adds
the matching control counts/fractions and two descriptive ratios:
``record_to_control_event_ratio`` divides the reference record fraction by the
control event fraction; ``key_to_control_event_ratio`` uses the reference
junction/gene-key fraction as numerator. The sampling units remain explicit.
Unresolved reference genes, genes not observed in the control, and absent control
lengths have separate status values. Ratios with zero or unavailable control
frequency are null; no pseudocount or cross-length fallback is applied.
These ratios describe composition relative to the supplied control and do not
establish antigen selection, structural contacts, binding confidence or a scorer
assignment. The manifest records the raw-control identity, retained/excluded rows,
event counts and normalizer source identity.

The same control enables ``structure.spectratype.tsv`` and the exact histogram
``structure.spectratype-lengths.tsv``. These use distinct full-junction/allele-pooled
V-gene/J-gene keys with equal weights in both populations, retaining unresolved
genes as explicit null categories. Their empirical CDF maximum distance is a KS
effect size without a p-value; Wasserstein-1 distance is measured in residues.
Signed mean and ordinary sample-median shifts subtract the control value from
the reference value. The control's inclusive central interval uses inverse-CDF
2.5% and 97.5% endpoints (smallest observed length whose cumulative fraction
reaches the probability). Reference tail mass counts only lengths strictly
outside those endpoints. Discrete length gaps and void bins need no smoothing,
or Gaussian assumption. ``length_kl_bits`` is KL(reference || control) in
bits per key, infinite if any reference length lacks control support.
``length_js_bits`` is equal-mixture Jensen–Shannon divergence, between zero
and one bit, including disjoint supports. These describe length selection;
neither field alone establishes above-noise class support.

``structure.degrees.tsv`` also exports ``inclusive_degree=degree+1``: self
counts once, so orphans have inclusive degree one. The accompanying
``structure.clonotype-degrees.tsv`` retains distinct junction/V-gene/J-gene
keys with unambiguous alleles pooled. Unresolved V/J categories stay in the
denominator and their counts are reported. Different genes sharing an exact
junction remain separate vertices. The legacy ``.degrees.tsv`` and
``hamming_*`` summaries retain unique-junction vertices.
``structure.graph-acquisition.tsv`` gives exact expected edge counts and
connected/orphan sampled clonotypes under uniform sampling without replacement
from this observed Hamming<=1 graph. It uses up to 25 geometric depths through
the full distinct-clonotype count. This is interpolation of observed graph
structure, without extrapolation to unseen families. It does not treat
overlapping neighbourhoods as disjoint taxa for Chao coverage.

The structure summary additionally reports similarity-aware ``graph_hill0``,
``graph_hill1`` and ``graph_hill2``. With uniform mass over N gene-level clonotypes and
binary self-inclusive Hamming similarity, ordinariness is ``z=(degree+1)/N``.
The three summaries are ``mean(1/z)``, ``exp(-mean(log(z)))`` and
``1/mean(z)``. Shannon-equivalent entropy uses natural logs; ``N/Hill`` is
reported as enrichment. An orphan graph has Hill=N/enrichment=1; a complete
graph has Hill=1/enrichment=N. No additional search or dense matrix is needed.
Separate ``degree_chao1`` and ``degree_chao_coverage`` apply existing vdjtools
estimators to the synthetic count vector ``degree+1``. Their exported
``degree_mass`` is the sum of these counts; degree singletons are orphans
and degree doubletons have one neighbour. ``degree_chao_status`` labels this
degree-mass proxy. Ordinary occurrence-based Chao coverage remains available
through ``reference-coverage``.

Add ``--publicity`` with ``--control`` to measure epitope record publicity by
one additional native seqtree batch. Full-junction Hamming distance at most
one includes exact matches and excludes indels. Control junctions are unique;
reference source-record multiplicity is retained. ``publicity_record_fraction``
is the fraction of valid reference records with at least one control match.
``publicity_mean_control_fraction`` divides the matching record/control-junction
pair count by valid reference records times distinct control junctions.
Both integer numerators and denominators are exported. These pooled-control
statistics do not estimate donor prevalence or determine germline origin.

Use ``--neighbor-draws 200`` for same-N background comparisons. N counts
distinct allele-pooled full-junction/V-gene/J-gene clonotypes. Uniform draws
without replacement use all valid same-species/locus VDJdb clonotypes and,
with ``--control``, observed control clonotypes, ignoring event/read abundance.
Each independent replicate supplies nested prefixes for every requested N;
draws across N are correlated. One native Hamming batch searches the union of
sampled junctions per background. Different clonotypes sharing an exact
junction remain neighbours; self is excluded.
``structure.neighbor-draws.tsv`` retains every replicate. The separate
``structure.neighbor-background.tsv`` gives observed values, draw means,
medians and 10th/90th percentiles for pair density, connected fraction and
mean degree, plus observed/mean ratios. Zero/zero is unavailable and a
positive value over zero is infinite. Empty or undersized backgrounds remain
explicit. These are random-background distributions, not bootstrap confidence
intervals or a knee-selection procedure. No Pgen matching is applied here.

The control also enables ``structure.motifs.tsv``. Select its projection with
``--motif-mode exact3`` (full-amino-acid 3-mers), ``wild4`` (full-amino-acid
4-mers with one wildcard at each of the four positions), ``imgt3`` (3-mers
in the eleven IMGT classes), or ``imgt_wild4`` (the existing default).
Each row labels ``motif_mode`` and ``motif_alphabet``; the manifest records
window size and mask positions. All views use full anchor-inclusive junctions
and once-per-junction/V-gene presence, with the same V/length denominators and
missing-control statuses. Their information values describe separate, dependent
projections and must not be added together. The bounded counter reuses the
existing tokenizer and rejects a count buffer above 256 MiB before allocation.
The separate ``idf`` and ``tfidf`` fields use control junction/V keys as
documents: ``idf=1+ln((M+1)/(df+1))`` and ``tfidf=(n/N)*idf``. Here ``M``
is the number of distinct control junction/V keys, ``df`` their global motif
document frequency, ``N`` the target's total junction/V-key count and ``n``
the keys containing the motif in the row's V/length stratum. Repeated windows
count once per document. No vector norm is applied. The same global IDF is
used for both length scopes; TF-IDF and binary selection bits remain separate.

Use ``--motif-mode full20`` for just full20 3-mers and single-X 4-mers at all
four mask positions. Add ``--epitope-idf`` to keep a second IDF based on the
complete same-species/locus reference, independent of ``--epitope`` filters.
Its documents are epitopes with at least ``--min-records`` valid source records.
Within each V/length stratum, ``epitope_documents`` counts eligible epitopes
and ``epitope_document_frequency`` counts those containing the motif at least
once. ``epitope_idf=1+ln((epitope_documents+1)/(epitope_document_frequency+1))``;
``conditional_epitope_tfidf=reference_fraction*epitope_idf`` uses conditional
motif presence. ``control_idf`` preserves the original control IDF, while
``conditional_control_idf_bits=log2(control_keys/control_with_motif)`` measures
conditional control rarity (infinite if absent; null if no control stratum).
These are distinct descriptive axes. Control commonness alone does not prove
germline origin, and epitope IDF alone does not establish held-out accuracy.

Add ``--mock-replicates 200 --seed 0`` to export ``control_tfidf_mean``,
``control_tfidf_lo`` and ``control_tfidf_hi`` in the same motif table. These
are the mean and 2.5/97.5 percentiles from the existing Pgen-matched control
sampler, with IDF fixed. Source-record counts are matched within one-bit Pgen
strata before both real and mock populations collapse to junction/V keys.
Joint junction/V/J event slots are sampled without replacement; V/J and length
are not matching constraints. An absent mock token/stratum contributes zero.
These are matched-control envelopes, not confidence intervals for selection
bits or final specificity probabilities. Replicates default to zero (off),
and otherwise must be at least two. By default one fresh native vdjtools exact
junction-Pgen batch covers the selected reference/control union. Per-target
seeds are stable under target ordering.

``reference-structure``, ``reference-coverage`` and ``reference-diagnostics``
also accept ``--pgen-features PREFIX.features.parquet`` with ``--control`` and
``--mock-replicates``. This explicitly reuses the requested durable Hamming-1
Pgen features written by ``pgen --control CONTROL --mismatches 1``. The companion
``PREFIX.manifest.json`` travels with the table; its filename and SHA256 bind
the feature data. It must identify species/locus, ``mismatches=1`` and
uncollapsed model files. The consumer verifies the complete consumed model-file
set against its bundled model and requires sorted unique canonical junctions,
finite probabilities in [0,1], and coverage of every valid reference/control
junction. It performs no Pgen recomputation when features are supplied. The
matching bins then use closed Hamming-1 probability mass. Reference bins absent from the controls,
including zero mass, share the lowest occupied control bin with enough event
slots for all mapped source rows plus its existing occupants;
the stored probabilities remain unchanged. The default fresh
route retains exact-junction mass. Diagnostic manifests record the feature and
producer-manifest hashes and model identity. These are reusable model features;
control draws and diagnostic calculations remain fresh within each invocation.

Compute fixed control reductions once with ``prepare-control --control RAW
--species human --locus TRB --threads 12 --output-dir NEW_DIRECTORY``. It stores
the normalized event table, gene counts and four motif-count models with raw-input
and file hashes. Supply ``--control-model NEW_DIRECTORY`` alongside the original
``--control RAW`` in later structure, coverage or combined diagnostics to reuse
these explicit prepared inputs. The consumer validates their population and
hashes; matched draws, reference features, graphs and bootstrap outputs are fresh.

With the same ``--control`` option, ``structure.genes.tsv`` also reports
``gene_bits``, ``length_bits`` and ``joint_bits``. These use distinct
junction/gene keys on **both** sides, rather than the event fractions used by
the two ratios above. Positive pointwise bits indicate enrichment and negative
bits indicate depletion. The companion ``structure.information.tsv`` contains
one row per full-pMHC/species/locus/segment and length scope:

.. list-table:: Information fields (bits per junction/gene key)
   :header-rows: 1

   * - Field
     - ``all_lengths``
     - ``exact_length``
   * - ``length_information_bits``
     - Zero
     - :math:`D_{\rm KL}(P_e(L)\Vert P_0(L))`
   * - ``gene_information_bits``
     - :math:`D_{\rm KL}(P_e(G)\Vert P_0(G))`
     - :math:`\mathbb E_{P_e(L)}D_{\rm KL}(P_e(G\mid L)\Vert P_0(G\mid L))`
   * - ``joint_information_bits``
     - Marginal gene information
     - Sum of length and conditional gene information

Here :math:`G` is V or J, :math:`L` is anchor-inclusive junction length,
:math:`P_e` is the target distribution and :math:`P_0` the empirical control.
V and J are separate projections. Counts include all original mass. Missing/nonfinite pointwise gene or length
contributions receive neutral zero-bit imputation, retaining per-row
``gene_bits_imputed``/``length_bits_imputed`` flags and original
``unresolved_comparison_keys``/``unavailable_length_keys`` counts. The summary
status is ``available_neutral_imputation`` when needed; other contributions
remain usable. The displayed KL identities apply to fully observed estimates;
imputed summaries are neutral-filled descriptive contributions. Missing MHC fields remain explicit
context values. Reference/key/control counts and the sampling units are retained.
These descriptive information estimates do not select a score or supply a
P-value. See :doc:`selection-information` for the likelihood and chain-rule
interpretation, and :doc:`control-calibration` for final-score tails.

The same control-aware run writes ``structure.motifs.tsv`` with reduced
single-wildcard four-mer presence within V genes, both pooled over lengths and
at each exact junction length. Each row retains motif-present counts and full
key denominators on both sides, ``presence_bits`` and the binary
``motif_information_bits``. The fixed eleven-class alphabet, all four masks and
full-junction windows match ``motif-evidence``. Repeated windows count once per
key; V alleles and alternative J annotations do not multiply a key. Only
reference-observed tokens are listed. Empty control strata and unobserved motif
presence/absence are explicit unavailable states; no infinite evidence or
pseudocount replaces them. The manifest counts short reference keys without a
four-residue window. J usage and motif position are not conditioned here.
The tokenizer runs in bounded blocks, and each target table is written directly
to the final output. Overlapping motif information is not added as independent
evidence or used to change the current scorer. See :doc:`selection-information`
for the Bernoulli likelihood-ratio interpretation.

The germline outputs ``structure.germline-map.tsv``, ``structure.germline.tsv``
and ``structure.germline-information.tsv`` project those V-key counts onto
individual CDR1/CDR2/CDR2.5 alignment columns, with full residues and eleven-class
alphabets kept separate. All bundled alleles must agree before a residue is
confirmed; missing or allele-dependent states retain their mass. Positions are
one-based model alignment columns. ``consensus_bits`` and
``consensus_information_bits`` describe the gene-consensus annotation.
``residue_bits_lower`` instead allows every compatible or unresolved control gene
in the background count; it bounds empirical residue enrichment over the
allowed allele assignments. Zero compatible control counts remain unavailable.
The original junction/V-key sampling unit is preserved after projection, so
merging genes with shared residues does not deduplicate their keys. These are
separate correlated views, without an inferred allele mixture, P-value or fitted
score choice. See :doc:`selection-information` for the bound and its assumptions.

The graph has one vertex per unique anchor-inclusive junction and an edge for
one amino-acid substitution, with no indels. One native batch searches the entire
selected reference population. The complete reference graph and its composition
remain fixed for all queries; no query-specific graph is rebuilt or exported.

Outputs report exact source-record, junction, clonotype-key, component, isolated-
junction and edge counts; largest-component and clustered fractions; component
entropy; size-weighted component size; junction-length mean/SD/entropy; and V-gene
entropy/effective diversity. Entropies use natural logarithms. Lengths include the
conserved junction anchors. V usage counts unique junction/resolved-V-allele keys
with allele-stripped gene categories. Missing or unmodelled records are counted
explicitly. An empty population has unavailable composition fractions and entropies.

``structure.degrees.tsv`` reports the exact Hamming-one degree histogram as flat
``degree,junction_count`` rows for each full reference population. These nodes
include every valid full junction, even when its V gene is unresolved. The summary
table adds ``hamming_`` counts, mean/variance/maximum degree, median, 10th and 90th
percentiles, and isolated fraction for that population. Quantiles use linear
interpolation. These degree distributions describe the fixed full reference; the existing modeled connectivity columns retain their stated
V-resolution and length requirements.

``structure.lengths.tsv`` reuses those same fixed graph nodes and edges. Each row
reports one exact ``junction_length`` in amino acids, ``valid_reference_records``,
``modeled_reference_records``, unique ``junctions``, ``graph_edges``,
``isolated_junctions``, ``mean_degree`` and ``isolated_fraction``. Duplicate source
records count separately from graph nodes; valid short junctions and unresolved
V genes remain in the graph. Modeled records additionally require a resolved V
allele and junction length of at least eight. The inherited ``source_observations``
and invalid/unavailable counts describe the whole pMHC population, not the length
stratum. An empty population has one row with null length, zero counts and null
degree/fraction. Hamming edges never connect unequal lengths. Thus a rare length
can have neighbours, and a missing neighbourhood does not establish nonrecognition.

``structure.genes.tsv`` counts V and J separately among all valid reference
records, including short junctions and unresolved genes. Alleles are pooled;
ambiguous, missing or wrong-locus calls remain the null gene category with
``gene_status=unresolved``. ``length_scope`` distinguishes ``all_lengths`` from
``exact_length``; a null length denotes only the full population.
``record_fraction`` divides ``reference_records`` by ``total_reference_records``.
``key_fraction`` divides distinct ``junction_gene_keys`` by
``total_junction_gene_keys`` within each segment and length scope, including the
unresolved category. Counts are not read abundance. These compositions are not
background-adjusted germline enrichment or evidence of direct structural contacts.

These are reference morphology features. The command fits no weights, makes no
recognition claim and supplies no P-value. Its manifest records population policy,
raw reference identity, feature definitions, source hashes and output hashes.

vdjmatch database-evidence
--------------------------

Report literal database presence and two control directions for a declared
full-junction statistic. Supply one raw control; Pgen, mock references and
same-Pgen query draws are computed on demand:

.. code-block:: sh

   vdjmatch database-evidence query.airr.tsv --vdjdb reference.zip \
     --locus TRB --species human --control control.tsv.gz \
     --scope 1,0,0,1 --statistic count --threads 4 --output-prefix presence

The command writes ``presence.database.tsv`` with every original query row,
``database_exact_junction``, ``database_exact_gene_clonotype``,
``database_nearest_edits``, ``database_score``, ``p_reference``, ``p_query``,
``p_both`` and ``database_status``.
``database_exact_gene_clonotype`` pools alleles and requires single supplied V/J
calls. It is unavailable for missing/ambiguous query gene calls. A false flag
means no complete matching gene key was found; an exact junction can still
occur in records with incomplete reference gene annotations. The statistic uses
full-junction unit edits and ignores V/J. ``--scope`` separately enforces
substitution, insertion, deletion and total budgets; otherwise ``--max-edits``
selects symmetric Levenshtein distance. ``--statistic nearest`` (default) uses
negative nearest distance; ``count`` uses distinct matching junctions.
A no-hit query receives both tails and their maximum equal to one;
invalid/mismatched query rows retain unavailable diagnostics.

``p_reference`` compares the observed query against mock references preserving
the selected reference's record counts and one-bit log2-Pgen composition.
``p_query`` compares the observed score against same-Pgen control queries on
the real reference. ``p_both`` is their maximum. Each direction uses
``--mock-replicates`` fresh draws (default 200), seeded by ``--seed`` (default 0).
Exact overlaps stay eligible in random draws. Read abundance is ignored;
control weights count available rearrangement events. Species and chain are
selected separately. Bundled native Pgen models are selected automatically.

``null_query_match_rate`` estimates the probability that a same-Pgen control
clonotype has any match within scope, and ``null_reference_match_rate`` estimates
the corresponding rate against mock references. These ``k/B`` estimates differ
from the ``(k+1)/(B+1)`` score-tail P-values; both counts and denominators are
exported. They are not posterior probabilities of a false specificity.

``--epitope``, ``--evidence-type`` and ``--min-score`` filter the loaded reference
before both calibration directions. Evidence types use the reference's own
declared vocabulary. A custom reference table uses the same loader and its
own selected source records and full-pMHC groups.

``presence.epitopes.tsv`` exports supported fixed-pMHC candidates; an omitted
valid query/target has no hit and both tails equal to one.
``presence.calibration.tsv`` records per-target Pgen strata, source/key counts,
control event counts, scope, statistic, seed, model/native identities and timings.
The control route writes these three flat TSVs, with no JSON calibration input
or output. One native population batch handles query and control draws; mock
references are streamed and reduced. Insufficient control strata raise an error.
Sparse join expansions are checked; native hit allocation remains a resource
limit for broad scopes. Native workers use ``--threads`` and table/BLAS helpers
stay serial. See :doc:`control-calibration` for the derivation and scope of the
conditional-null interpretation.

The advanced legacy ``--null-sample`` route retains its original unconditional
query-null diagnostic and ``database_p``/null-table/manifest outputs. It does not
provide the two new conditional directions. Neither route calibrates the
V-loop/combined VDJmatch score, a winning antigen label or genuine paired
recognition; those require their complete scoring and selection adapters.

vdjmatch reference-coverage
---------------------------

Install ``vdjmatch[diagnostics]`` for the seeded Leiden adapter. Compute database
occurrence coverage at exact-clonotype and Leiden-community resolution:

.. code-block:: sh

   vdjmatch reference-coverage --vdjdb reference.zip --species human --locus TRB \
     --epitope DEMO --geometry levenshtein --scope 1,1,1,1 \
     --nboot 200 --seed 0 --threads 4 --output-prefix coverage-edit
   vdjmatch reference-coverage --vdjdb reference.zip --species human --locus TRB \
     --epitope DEMO --geometry sw --sw-radius 40 --cdr-loops cdr12plus \
     --gap-open 11 --gap-extend 1 --nboot 200 --seed 0 --threads 4 \
     --output-prefix coverage-sw

``DEMO`` and radius 40 illustrate syntax; the command does not select a
scientific radius or a sufficiency cutoff. Supply the intended graph scale.
``--membership-budget-bytes 0`` disables the optional estimated raw/mock
membership abort; its default is 2 GiB. This budget does not change graph
membership, mock draws or reported statistics.
The edit graph uses the declared substitution/insertion/deletion/total budgets;
insertion and deletion budgets must be symmetric for an undirected graph.
SW uses full local BLOSUM62 alignment independently for the full junction and
the selected germline CDR1/CDR2/CDR2.5 loops, summing
``SW(a,a)+SW(b,b)-2*SW(a,b)`` across loops. ``--cdr-loops cdr12`` omits CDR2.5;
``junction`` explicitly requests junction-only SW. There is no trimming, stem4
or bounded-kernel substitution. Germline templates reuse the existing allele
resolver and flagged complete co-locus V-template imputation. Original V calls
and clonotype identities remain retained alongside the imputation flags. Gaps cost
``gap_open+(length-1)*gap_extend`` in raw BLOSUM62 score units.

Vertices are distinct full-junction/V-call/J-call keys, with supplied alleles
retained. Every valid source occurrence remains an abundance count. Leiden
uses unweighted undirected graphs, RBConfiguration modularity, the explicit
``--leiden-resolution`` (default 1), a stable seed and iterations until stable.
Isolates remain singleton communities. Organism, chain and full pMHC stay
separate. ``--targets-from-sample`` uses the same compatible-union reader as
``reference-structure``. ``--reference-population`` has the same curation choices.

The five final TSVs are ``.populations.tsv`` (source/retained/excluded counts),
``.communities.tsv`` (vertex keys, occurrence counts and memberships),
``.edges.tsv`` (vertex IDs and the inclusive ``distance_ceiling``, not exact
edge distances), ``.coverage.tsv`` (Chao1 richness, coverage, unseen mass,
unseen-type fraction and entropy) and ``.rarefaction.tsv`` (richness and coverage
against occurrence effort, with uncertainty). The latter also supplies richness
against coverage without an additional calculation. See :doc:`selection-information`
for the estimators. The manifest records inputs, model/source/native identities,
geometry, partition parameters, seeds and versions. Missing V/J calls remain
explicit in the exact key; full SW with germline loops requires resolvable V.

Add ``--control RAW_CONTROL --mock-replicates 200`` for Pgen-matched reference-size
control populations, using the same event-slot law as the motif diagnostic.
Each mock gets its own graph partition; each population's occurrence bootstrap
then holds its complete partition fixed. One native graph batch per invocation
serves all real/mock populations, with unique search keys shared during that
operation. Full SW emits bounded sparse native hits within each fixed population, without
a dense control matrix. The CLI spools edge chunks to temporary population files,
loads one sorted population for Leiden, and removes scratch files on exit. The
final edge table retains every edge. The default ``--pair-budget 2000000`` limits candidate
pairs before search: within-population unordered pairs for sparse SW, union
unordered pairs for indexed edit search, and query-by-union pairs for the knee
count traversal. Raise it explicitly
for a declared larger workload. Memberships, graph
expansion and Chao bootstrap buffers have separate allocation guards. A small
synthetic check does not establish full-cohort time/memory acceptance.


vdjmatch reference-density
---------------------------

Export kth-other-clonotype SW distances and pooled epsilon knees without
refitting communities:

.. code-block:: sh

   vdjmatch reference-density --vdjdb reference.zip --species human --locus TRB \
     --epitope DEMO --cdr-loops cdr12plus --threads 4 --output-prefix density

The default neighbour orders are 1, 2, 5 and 10; repeat ``--neighbor-order K``
to choose others. All orders share one native cumulative-count batch. Distinct
full junction/V/J clonotypes are neighbours only within the same exact pMHC;
self is counted once and excluded by the condition ``count > K``. Zero-distance
ties remain neighbours. Unresolved geometry, insufficient neighbours and missing
confirmed knees have explicit null/status outputs. Knees pool the selected
species/chain distances across its exact pMHC groups. This diagnostic leaves the
community graph's existing k=1 radius rule unchanged.


vdjmatch reference-pairing
---------------------------

Join genuine linked TRA/TRB observations to fixed real community memberships:

.. code-block:: sh

   vdjmatch reference-pairing --vdjdb reference.zip --species human \
     --epitope DEMO --geometry sw \
     --alpha-communities alpha-sw.communities.tsv \
     --beta-communities beta-sw.communities.tsv --output-prefix pairing

The command verifies each community table and companion manifest against the
reference, species, chain, geometry and population policy. Joins retain exact
pMHC and junction/V/J keys. Outputs separate complex links, distinct paired
clonotypes, community pairs and all marginal chain clonotypes, including unlinked
and unavailable memberships. Pairing fractions report source-record and
clonotype denominators separately; missing linked evidence leaves marginal
support intact. Compatible-MHC union tables are rejected.


vdjmatch motif-position
-----------------------

Export a descriptive positional amino-acid contrast against an existing prepared
control within one explicitly chosen species/chain/epitope/V-gene/J-gene/junction-length
stratum. MHC contexts are pooled; reference source records retain their metadata.
This prepared rearrangement control is empirical: its frequencies include the
biology and ascertainment of the observed repertoire, rather than representing
a pure V(D)J generation law.

.. code-block:: sh

   vdjmatch motif-position --vdjdb reference.zip --species human --locus TRB \
     --epitope GLCTLVAML --v-gene TRBV6-1 --j-gene TRBJ2-1 --junction-length 17 \
     --control raw-events.tsv.gz --control-model prepared-control --motif GIV \
     --threads 1 --output-prefix motif-position

Add ``--generative-background olga`` to append a separately labelled exact
generation contrast. The immutable bundled human OLGA model is loaded without
allele collapsing; every allele pair belonging to the selected genes contributes
its joint Pgen. One native masked transfer-matrix batch evaluates the all-wildcard
length denominator, every position/residue numerator, and each optional exact
motif start. ``generative_fraction`` is numerator Pgen divided by the common
V/J/length Pgen; ``generative_difference`` is reference fraction minus that
conditional probability. The model sums generation-event probability mass;
it does not assign equal weight to all possible junction strings. No simulation,
refitting or model regeneration occurs. The bundled OLGA source currently
supports human; other organisms fail explicitly.

The optional motif's anywhere union has no generative estimate (null,
``unsupported_anywhere_union``), because overlapping starts cannot be added as
disjoint events. Each start has its own exact conditional estimate. A zero
generation denominator produces ``unsupported_generative_stratum`` with null
conditional probabilities and differences. Model files, allele lists, native
vdjtools binary, DP wrapper and loader identities are pinned in the manifest.

The existing prepared-model loader validates the original raw-control identity,
manifest, population and event table. Nothing is generated or fitted. Both
populations count unique full-junction/V-gene/J-gene keys with equal weight;
read abundance and prepared event weights are ignored. Alleles are pooled by
the shared gene normalizer (uppercase unambiguous locus-matching calls, allele
suffix removed); ambiguous calls cannot enter the selected stratum. The default
``curated-benchmark`` population retains existing curation exclusions.

``.positions.tsv`` contains all 20 canonical residues at every 1-based full-junction
position, including conserved anchors, with observed and control counts,
denominators, proportions and their signed difference (observed minus control).
No pseudocounts or log ratios are used. A zero control denominator produces
``unsupported_control_stratum`` and null control fractions/differences; an empty
reference stratum fails explicitly. Optional ``.motifs.tsv`` reports an exact
motif at every possible start and at position zero for presence anywhere, counting
each key once. ``.manifest.json`` records input/output hashes, arguments, counting
rules and software/native identities. These are background frequency contrasts,
not fitted ALICE selection coefficients, calibrated significance or binding claims.


vdjmatch reference-diagnostics
------------------------------

The classification unit is species/chain/epitope. This command pools all MHC
contexts by default, labels each epitope with the most frequent observed MHC-A/B
pair across both chains, and requires at least 30 valid source records per chain
(``--min-records``). Ties use lexical order; MHC class is the mode within the
chosen pair. Original MHC fields remain in source provenance. The small epitope
tail is reserved for C-index and rarefaction continuum, not singleton profiles.
Flat ``.eligibility.tsv`` and ``.epitopes.tsv`` record counts and labels.

The standalone structure, coverage, density and pairing commands accept
``--epitope-profiles`` for the same pooling/eligibility rule. Their legacy exact
MHC behavior remains available for earlier calculations.

Run all four motif views and both occurrence/community geometries in one fresh
invocation. The raw control, native Pgen and immutable event-pool preparation
are shared in memory only for that invocation; no computed intermediates are saved or reused.

.. code-block:: sh

   vdjmatch reference-diagnostics --vdjdb reference.zip --species human --locus TRB \
     --epitope DEMO --control raw-events.tsv.gz --motif-mode all --geometry both \
     --sw-radius-knee --cdr-loops cdr12plus --pair-budget 8000000000 \
     --mock-replicates 200 --nboot 200 --seed 0 --threads 12 --output-prefix diagnostics

``--sw-radius-knee`` selects one radius from the real reference only. Each
modeled full-junction/V/J vertex contributes its nearest other vertex's full
SW distance within its full pMHC stratum. The distances are pooled over the
selected species/chain reference, sorted, and passed to the existing convex
Kneedle helper (S=1, strongest confirmed interior knee, no smoothing). Singleton
strata and unavailable V loops have no nearest distance. If no knee is confirmed,
the command fails rather than inventing a radius. The selected radius remains
fixed for all matched controls. It is a geometric diagnostic, not a specificity
or adequacy threshold. Numeric radii from other geometries cannot be transferred to raw
SW units.

The reference nearest distances and full strata are exported as ``.sw.knee.tsv``.
Other coverage outputs receive ``.levenshtein`` or ``.sw`` prefixes; structure
outputs retain the main prefix. All four motif views are explicitly labelled.
The sparse-hit and compressed-bootstrap routes require native seqtree and
vdjtools builds exposing those capabilities; version strings alone do not
establish that support.

vdjmatch reference-rule-scores
------------------------------

Rebuild all 15 retained fixed-rule scores from raw VDJdb for one species and
chain (or genuine paired receptors with ``--locus paired``), pooling MHC by
epitope. Eligibility uses valid source records before
exact-hit filtering; sparse labels remain search candidates. Models remain
fixed. Unresolved V calls receive the fixed co-locus template used by
``neighbours``, with original V/J query identities and imputation flags retained.
Unavailable junctions and absent linked target populations remain explicit.

.. code-block:: sh

   vdjmatch reference-rule-scores --vdjdb reference.zip --species human \
     --locus TRB --threads 4 --output-prefix fresh-rules

``fresh-rules.scores.tsv`` contains every winning label, tie and no-hit for
``score-confusion``. ``fresh-rules.continuous/`` contains flat Parquet parts
for every query, eligible epitope and scorer, including unavailable cases.
Queries, complete known truth memberships, target eligibility and a source/raw/
native/output hash manifest accompany these tables. Existing output paths are
rejected. Add ``--components --germline-background PRIOR`` for the six retained
single-chain components (21 scorers total). Supply the accepted matching density
``--control`` except for the bundled human TRB control; the prior is always
explicit. These raw model inputs are hashed before and after scoring.
Add ``--kernel-only`` to export just the six retained current stem4 kernel
alternatives from one native traversal, using the fixed VDJdb-other background.
This skips the fixed rules and historical components. The geometry is unchanged:
BLOSUM62, trim3 junction ends, stem depth4 and the standard kernel ladder.
These are continuous relative evidence scores; they have no transferred
P-value or E-value calibration. All six exports exactly match the detailed
``neighbours`` CLI under the same MHC-pooled population.

For a baseline-preserving candidate, add ``--v-model-policy resolved-first``
with ``--kernel-only``. Queries whose original V calls resolve use the fixed
original resolved-allele reference population; queries with unresolved V calls
use the fixed extended imputed population. Paired queries require both original
V calls to resolve for the first model. Original query identities, every eligible
case and explicit unsupported target models remain in the exports. Synthetic
fallback alleles therefore cannot alter established resolved-query rankings.
This deliberately forgoes unresolved-reference support for resolved queries.
Each population uses one native batch; neither model is refitted per query.
The default remains ``extended``. The manifest records the policy; existing
extended-model final-call calibration does not apply to ``resolved-first``.

``--position-profile PROFILE.tsv`` with ``--kernel-only`` scales substitution
costs per junction column instead of weighting every column equally. The file is
a TSV of ``locus``, ``side``, ``offset`` and ``multiplier``. ``side`` is ``cys``
for offsets counted from the conserved Cys104 anchor and ``fw`` for offsets
counted back from Phe/Trp118; both sides are required for a declared locus and
each runs contiguously from offset 0. The last declared offset is the pooled
core, reused by every deeper position, as in the end-anchored
``position_significance.tsv``.

The two sides combine by **minimum**, the convention of the end-anchored
``significance_weights`` whose ``omega = 1 - max(p_same_V, p_same_J)`` says the
same thing about costs rather than proximities: either anchor's germline
constraint is enough to hold a column down. The two counts therefore converge
once the junction is longer than the two declared depths, when every interior
column carries the minimum of the two pooled cores. Anchoring at both ends rather than tabulating
``(length, position)`` is what makes the file estimable and complete: valid human
TRB carries 105,677 source records over 29 junction lengths, 12 of which hold
fewer than 30 records each and 6 exactly one, so a per-length table would both
demand numbers no cohort supports and silently miss any length it did not
enumerate. Pooled anchor offsets cover every length, including unseen ones.

Multiplier ``1`` leaves a column unchanged, so an all-ones profile reproduces the
unweighted scores exactly. Only CDR3 substitutions are scaled: gaps, V-loop
offsets and stem offsets are not, and the search radius and cutoff ladder stay
derived from the substitution-matrix maximum. Normalising the file to mean 1
over scored positions preserves average nominal cost scaling. It does not
preserve neighbour identities: substitutions at upweighted positions become
more costly, while downweighted positions become cheaper. Only an all-ones
profile guarantees equivalence with the flat geometry. The manifest records the profile hash, its anchor depths
and its multiplier range. This is an opt-in experimental geometry; the default
and its calibration are unchanged. Explicit matrix and profile files are both
hash-bound and checked for changes during scoring.

.. code-block:: sh

   vdjmatch reference-rule-scores --vdjdb reference.zip --species human \
     --locus TRB --kernel-only --matrix blosum62 \
     --position-profile contact-mask.tsv \
     --reference-population curated-benchmark --min-records 30 \
     --threads 8 --output-prefix arm-contact

``--class-policy MODEL/reference-policy.tsv`` exports only the experimental
``reference_class_policy`` scorer. It is exclusive with ``--components`` and
``--kernel-only``. The sibling ``manifest.json`` must bind the policy table hash,
its reference-only provenance, the existing ``FAST_ROUTES`` and the same raw
reference identity. Both model files are hashed before and after fresh scoring.
The annotation default is unchanged; no query labels fit or tune these routes.

Existing fixed-rule and kernel traversals supply V and core-exact3 scores plus
``kernel_v_bounded`` as the configured default. The unchanged
``apply_reference_policy`` selects V for ``v_only``, exact3/V for ``motif_only``,
exact3 for ``v_and_motif`` and the configured default otherwise. Missing tail
profiles, unavailable components and paired routes retain the existing finite
default behavior. Kernel zero/no-hit support stays unchanged. No cross-axis
fusion or component P/E transfer is introduced; final calibration remains open.

Every valid tail model participates in winner selection. Eligibility restricts
query sampling and continuous exports, so an unprofiled tail can win or tie
against primary labels. Continuous parts include ``reference_class``,
``selected_score_column`` and ``policy_status``; ``.policy-audit.tsv`` counts
these states over all query/target coordinates, including tails and unavailable
models. Only routed component columns are retained during traversal; the other
component scores are not exported as new results.

.. code-block:: sh

   vdjmatch reference-rule-scores --vdjdb reference.zip --species human \
     --locus TRB --class-policy MODEL/reference-policy.tsv \
     --threads 4 --output-prefix fresh-class-policy

In paired mode, eligibility uses the maximum valid chain source-record count;
``paired_source_records`` reports genuine complete complex observations
separately. Eligible labels without linked references remain explicit null
coordinates. Pair query identities retain both junctions and supplied V/J calls.
The fixed-rule and current-kernel exports score the complete-pair component;
native marginal evidence and label fusion remain separate. Historical
``--components`` supports single chains only. These exports do not establish
final-call calibration or the complete manuscript resource gate.

vdjmatch score-confusion
------------------------

Reduce complete scorer candidate outputs to one-vs-other annotation-retrieval
confusion tables. This command does not score receptors, refit models or choose
thresholds. Benchmark score producers must filter matching exact hits only and
retain fixed reference models. Ordinary exact-retained annotation is a separate
population.

.. code-block:: sh

   vdjmatch score-confusion candidates.tsv --truth truth.tsv \
     --profiles profiles.tsv --threads 1 --output-prefix retrieval

Candidate columns are ``species,locus,query_id,epitope,scorer,score,supported``.
Scores increase with evidence; ``supported`` is an explicit Boolean. Every
scorer must represent every truth query in every supplied context, including
no-hit rows (null epitope/score, false support). Truth columns are
``species,locus,query_id,epitope`` with every known membership retained.
Multiple truth labels are allowed. All supported finite maximum-score labels
win, including exact ties. Unsearched or unsupported queries cannot disappear.

``retrieval.confusion.tsv`` gives TP/FN/FP/TN, accuracy, balanced accuracy and
MCC per scorer and pooled epitope with at least 30 valid source records.
Undefined denominators produce null metrics. All query, positive/negative,
no-hit and tied-query counts are explicit; ``retrieval.calls.tsv`` retains
winning labels and no-hit dispositions. When full index profiles are supplied,
the metric table includes alpha, beta and joint marginal classes for plotting.
An other-class retrieval label is not proof of experimental nonbinding.

Paired target eligibility requires at least one chain with 30 valid source
records, before exact-hit exclusion; two smaller chain counts cannot be added
to qualify a target. ``valid_source_records`` then reports both chains' total
source observations and ``eligibility_chain_max_records`` names the qualifying
quantity. Missing marginal chains remain explicit in joint class metadata.

``retrieval.class-groups.tsv`` supplies the heatmap summaries: alpha classes
for TRA scores, beta classes for TRB scores, and joint classes for paired scores.
Organism species and core3/full20 definitions remain separate. Mean and median
accuracy, balanced accuracy and MCC give every epitope equal weight; minimum,
maximum and the number of defined epitopes accompany each. Exact summed
confusion counts and query-target appearance counts remain available. A null
or nonfinite metric contributes no defined value and is not replaced by zero.
The individual-epitope confusion table is always retained.

Class means also include an Efron nonparametric bootstrap, controlled by
``--nboot`` (default 1000, minimum 200) and ``--seed`` (default 0). Each
replicate samples N epitope units with replacement from the N observed epitopes
in that class group. The same multinomial multiplicities apply across scorers
and metrics. Scores, reference models and class memberships stay fixed: the
estimand is the conditional class-group macro mean, not uncertainty in model
training or class selection. Graph edges, overlapping motif tokens and repeated
source records are not treated as independent bootstrap units.
``bootstrap_se_*`` uses sample standard deviation (ddof=1), ``bootstrap_bias_*``
is bootstrap mean minus observed mean, and ``bootstrap_lo_*``/``bootstrap_hi_*``
are 2.5/97.5 percentiles. ``bootstrap_valid_*`` and ``bootstrap_status_*`` retain
undefined-resample information. Intervals are unavailable if fewer than two
epitopes have a defined metric or any replicate is undefined; no zero fill or
silent rejection of resamples is used. This empirical bootstrap is separate
from label-permutation nulls, control draws without replacement and the Chao
augmented-assemblage law used for unseen-taxon coverage.
Method: `Efron (1979), Bootstrap Methods: Another Look at the Jackknife
<https://doi.org/10.1214/aos/1176344552>`_.

Linked-first paired labels
~~~~~~~~~~~~~~~~~~~~~~~~~~

``paired-labels paired-evidence.tsv.gz --selection selection.tsv
--input-manifest native.manifest.json --output-prefix labels --threads 4``
reduces existing all-target stem4 VDJdb-other kernel outputs without rerunning
search, controls or models. Positive finite linked evidence retains every
linked maximum peptide. Only absent linked support permits the union of
alpha/beta maximum peptides; disagreement remains ambiguous. This label policy
does not compare raw magnitudes across axes or change the any-call scalar.
Native maxima, tie counts and lexical winners must match the evidence exactly.
The reader projects maximum rows before streaming collection.

``.calls.tsv`` repeats unchanged native query metadata for every selected
peptide, or emits one null-label row for no support. ``.scores.tsv`` is accepted
by ``score-confusion``: its score is a selected-label indicator (1, otherwise
null), **not continuous evidence for AP/ROC, a probability or calibrated
confidence**. The manifest pins inputs and sources; existing any-call
calibration does not establish confidence for this label policy.


Strict repertoire enrichment
----------------------------

``repertoire-enrichment`` evaluates each repertoire independently against fresh
Pgen-matched real-control draws. The metadata TSV must have unique ``sample_id``
and ``path`` columns; paths resolve relative to the metadata file. Input counts
are read counts. Source receptor records and their counts remain distinct even
when their amino-acid junctions coincide. Invalid amino acids and other loci are
counted in the ingestion table; noncanonical anchor residues and zero model Pgen
do not themselves discard records. Same-gene V alleles are pooled. Missing or
unresolved V uses the label-independent observed V-gene mode, falling back to the
bundled model mode; missing locus uses the requested locus. Original source calls
remain in the raw input, and imputed-record counts accompany ingestion. Control
event identities retain distinct original missing-gene keys while preserving the
established allele-pooled event law for complete inputs.

.. code-block:: console

   vdjmatch repertoire-enrichment samples.tsv --vdjdb reference.zip \
     --control control.tsv.gz --species human --locus TRB --threads 4 \
     --mock-replicates 200 --seed 0 --output-prefix fresh

The reference index is built once. One native batch searches each real repertoire
or each freshly drawn mock, with full-junction Hamming distance at most one and
matching V genes (alleles pooled). Only one repertoire is loaded at a time; the
whole control bank is never searched. Sample seeds are the requested seed plus
the zero-based metadata row number and are recorded in the manifest. Ordinary
exact reference evidence remains present. A receptor contributes once per full
pMHC regardless of repeated reference observations. Both matched receptor-record
counts and their summed read counts are compared with the mock distributions.

Every mock has the source repertoire's record count, the same rounded-log2-Pgen
histogram and the same read-count multiset within every bin. Bins follow the
floored rule ``round(log2(Pgen + 2**-60))``, with ties to even, applied
identically to real and control probabilities. The existing
control-event slot sampler draws without replacement, then randomly assigns the
source count roles. An absent or undersized bin fails explicitly; it cannot be
remapped or fixed by dropping source observations. Zero Pgen falls in bin -60.
Pgen bins use inclusive Hamming≤1 full-junction probability by default.
``--pgen-mismatches 0`` selects exact-sequence probability; the supplied
feature bank and fresh calculations must use the same definition.
``--pgen-features`` accepts the validated durable control feature bank with the selected probability definition; only
additional query masses are computed, in memory, without writing a feature cache.

Outputs are ``.summary.tsv``, ``.mocks.parquet``, ``.matching-audit.tsv``,
``.ingestion.tsv`` and ``.manifest.json``. Each endpoint's empirical upper-tail
P-value is ``(1 + number of mocks at least as large as observed)/(B + 1)``, with
``B`` mocks. At 200 mocks its minimum is 1/201. Lower/upper fields are the central
95% mock interval, not a confidence interval over biological donors. Enrichment
ratios are undefined when the mock mean is zero. No detection threshold,
multiple-testing correction or comparison between samples is imposed.

The manifest pins raw inputs, source/native bytes, model inputs, prepared-feature
identity, model-preprocessing source hashes, seed, stages and output hashes. Use an external summed-process-tree RAM
monitor for a large research cohort; this command does not reuse predictions,
mock samples or fitted nulls from another run. Row-preserving tcR input uses the
owning vdjtools ``read_mitcr(collapse=False)`` capability when available.

For durable probability preparation, ``pgen INPUT --features-only`` writes the
unique columnar ``.features.parquet`` bank and manifest while omitting the large
row-wise TSV. Existing probability columns are ignored. ``--mismatches 0`` and
``--mismatches 1`` produce separate, explicitly defined probability banks using
vdjtools; the normal row-wise output remains the default.

Streaming complete kernel evidence
----------------------------------

``evidence INPUT --vdjdb REFERENCE --locus paired --all-targets
--evidence-only --stream-evidence --gzip-evidence --output-prefix fresh``
uses the existing native batch consumer. Every detailed evidence field and
full pMHC row is retained. Selection, ties, competing scores, exact-hit policy
and radius-based query dispositions match the retained-table workflow. Rows
follow native execution order; compare or join by query/full-pMHC keys.

Paired streaming emits the linked ``.enrichment.tsv.gz`` and independent
``.TRA.enrichment.tsv.gz`` / ``.TRB.enrichment.tsv.gz`` tables, plus unchanged
compact selection and query tables. It omits the eagerly merged
``.paired-evidence`` table; all three complete axes remain available.
Marginal-only targets remain present. This changes storage, not scoring,
model populations or zero-background semantics. The retained detailed path
fails at8GB of accumulated table buffers and directs large jobs to streaming.

``examples/profile_cli.py`` retains summed RSS measurements and enforces its
memory cap against the larger of RSS and macOS physical footprint per process.
The latter includes compressed/swapped allocations; RSS alone can miss an
overflow. Resource reports expose both totals.

Missing control gene metadata
-----------------------------

V-aware control-event preparation imputes a missing V/J call with the selected
population's modal supplied call (lexical tie-break), or the bundled vdjtools
model mode when no call is supplied. Original normalized gene annotations remain
part of event identity before weights are aggregated. Affected outputs carry
``v_imputed``/``j_imputed`` indicators. Analytic control ingestion retains these
records for imputation. Input files and source event weights remain unchanged.
Explicit conflicting locus annotations still fail validation.

Fresh reference matching
------------------------

``match --fresh-reference`` builds its native index directly from the raw
reference in memory, bypassing prepared bundles without reading or writing them.
It supports the existing fixed and graded-ball modes, including exact-hit
exclusion, and records the fresh-reference policy in each output manifest.

With ``match --search-mode ball --no-hit-table``, reference observations are
expanded and reduced in bounded blocks of complete queries after the same one
native full-population search. Candidate values, exact-hit filtering and no-hit
calls are retained. This avoids retaining the full wide reference-hit table;
native hit memberships are still held until their blocks are consumed, so actual
whole-process-tree memory must still be measured for a large cohort.

For epitope-only reference diagnostics, ``match --fresh-reference
--epitope-profiles`` pools MHC contexts and annotates each species/epitope with
its modal observed MHC pair across both chains. All valid tail labels remain
in the fixed model (internal minimum one record); query eligibility is selected
separately. Original MHC fields remain in source metadata and the manifest
records the pooling rule and inventory. Ordinary matching keeps full pMHC.

Hard-gene reference diagnostics also impute missing/ambiguous reference calls
once for the fixed locus population, retaining original calls, flags, source
identities and junctions. All queries share that same annotation view; neither
exact-hit exclusion nor query membership changes reference model counts. The
ordinary source model remains unchanged.
