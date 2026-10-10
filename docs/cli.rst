Command-line interface
======================

The console script provides ``update``, ``match`` and ``precursor``. Use
``vdjmatch match -h`` for the installed version's complete options.

vdjmatch update
---------------

Fetch and retain a latest or pinned VDJdb release input for CLI and API use.

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
     - Require matching V/J calls; sequence-only control calibration cannot support
       these predicates, so combine with ``--no-evalue``. Paired matching rejects
       both flags until chain-specific V/J predicates are supported.
   * - ``--no-evalue``
     - Use explicitly uncalibrated ranking without requesting a background.
   * - ``--no-align``
     - Omit per-hit alignment and CIGAR computation for single-chain and paired output.
   * - ``--threads``
     - Native sequence-search threads; default ``1``. ``0`` requests native automatic
       selection. This does not reconfigure an already initialized Polars thread pool.

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

   **Which model set.** ``--source`` picks between three bundled sets, and the choice matters more
   for a mass over a *set* than for scoring one sequence, because a junction whose ``Pgen`` is
   exactly zero contributes nothing and raises nothing.

   ``olga`` (default)
      A bit-faithful import of OLGA's published models — ``Pgen`` matches OLGA's own to machine
      precision. Human, seven loci. Faithful includes faithful to OLGA's deletion-bin grid, on which
      an allele shorter than the grid carries probability on trims it cannot reach; **5.7% of human
      TRA junctions in VDJdb score exactly zero** as a result, and for some epitopes it reaches 13–15%.
   ``learned``
      Refit from real 5'RACE reads on arda germline, so it does not inherit that grid. Human, seven
      loci. 0.5% on TRA.
   ``arda``
      The same refit on the arda IMGT allele namespace, and the only set with a non-human organism:
      human for seven loci plus mouse TRA/TRB. 0.0% on both human loci and both mouse loci.

   Use ``olga`` when exact agreement with the reference implementation is the point. Prefer
   ``learned`` or ``arda`` for precursor work, and use ``arda`` on both species when comparing them
   so the model family is held fixed. The ``n_zero_pgen`` output column reports the loss per group,
   and the CLI warns when it exceeds 1% overall.

.. note::

   ``--species`` selects the **records**, ``--organism`` selects the **model**, and they have to
   agree. Scoring mouse junctions against the human recombination model does not error and does not
   produce zeros — it returns a plausible number for every group, a median 0.157× of the right
   answer with an 11× spread across epitopes, so the ranking is wrong too. That mismatch is
   therefore refused. Mouse models live only in ``--source arda``, which ``--organism mouse``
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

``--search-mode ball`` reports five fixed edit balls (radii 1–5). Unit-cost
balls share one radius-5 native batch per locus and control. Weighted scoring
uses one complete query batch per radius: its best alignment can use more edits
than the minimum-edit alignment, so a broad weighted hit cannot simply be
filtered into smaller balls. Each ball permits at most two insertions/deletions
and caps total edits at its radius. Matrix penalties rank accepted neighbours.
No global neighbour-count cap truncates the evidence. The calls table has one
row per retained query and radius; it does not select the radius with the smallest
P-value. Ball output can be large and does not include fixed-mode CIGAR strings.

.. code-block:: bash

   vdjmatch match queries.tsv --vdjdb reference.zip --search-mode ball \
     --threads 4 --output-prefix results/balls
   vdjmatch match cells.tsv --vdjdb reference.zip --paired --scope 3,1,1,3 \
     --threads 4 --output-prefix results/paired

Paired receptors currently use fixed mode with an explicit edit scope. Graded
paired balls require a separately validated joint calibration before use.
Both fixed and graded single-chain modes accept ``--match-v`` and ``--match-j``
with ``--no-evalue``. Gene matching strips allele suffixes and preserves IMGT
slash-containing gene names. Reference sizes follow the same gene predicates.
Sequence-only controls cannot calibrate these restrictions; calibrated requests
fail explicitly. Repertoire-level statistical inference remains future work.

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

``--mhc-match compatible --pool-reference`` pools intersecting family and allele
restrictions before counting unique reference junctions. Source restriction
labels remain intact in the reference. ``--targets-from-sample`` instead selects
each row's declared peptide, both MHC chains and class; it preserves observation
identity and reports unavailable targets explicitly.

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
Each additional edge contributes w(V) exp(-total TCRdist / tau) divided by the
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
accepted-edge total, weighted CDR3 and V-loop distances retain native TCRdist
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
and rejects older implementations before loading controls. TCRdist-specific
radius/temperature options do not apply.

Native TCRdist3-compatible distances
------------------------------------

``vdjmatch tcrdist-neighbours SAMPLE --vdjdb REFERENCE --locus TRB --threads 4
--radius 90 --exclude-exact --output-prefix OUT`` exports minimum distance and
neighbour counts per pMHC, plus a separate query-availability table. It reproduces
the default TCRdist3 distance: three times the trimmed, restricted-gap CDR3
distance, plus aligned V-gene CDR1, CDR2 and CDR2.5 distances. The bundled loop
models and their license/provenance are under ``resources/tcrdist``. Use
``--species human|mouse`` to choose the organism. Missing allele
suffixes explicitly resolve to ``*01``. Model-proven co-locus aliases resolve
uniquely without changing raw calls; ambiguous calls remain unavailable. The
query table retains raw calls separately from resolved alleles, including chains
with missing junctions. Unknown alleles and junctions shorter
than eight residues remain unavailable in the query table. This command reports
distances, not specificity probabilities.

``--locus paired`` accepts linked AIRR rows with ``clone_id`` and declared
TRA/TRB roles. It sums alpha and beta distances against the same linked database
receptor; ``--radius`` applies to that sum. Exact exclusion removes a reference
only when both full junctions match. Original independent marginal minima must
be computed from separate full single-chain reference pools; they are a distinct
baseline. Paired output retains missing partners and unavailable gene calls.

``--targets-from-sample`` selects compatible declared restrictions and counts
unique junction/V search keys within their union. Query status distinguishes
missing references, unavailable sequence/gene inputs, and searches without
neighbours. Explicit ``--epitope``, ``--mhc-a`` and ``--mhc-b`` selectors support
``--mhc-match compatible``; exact selection remains the default.

The existing seqtree gap-block batch kernel computes CDR3 distances; vectorized
lookups add germline-loop distances. Single-chain distance outputs are bounded
to 64 MiB per batch; paired matrix work uses a 32 MiB budget. These allocation
bounds do not establish whole-process RSS. A clipped BLOSUM-derived matrix reproduces TCRdist penalties; ordinary
seqtree BLOSUM Gram penalties are a different metric.

Exhaustive TCRdist nearest reference
------------------------------------

``tcrdist-neighbours --targets-from-sample --unbounded-nearest`` reports the
minimum native TCRdist distance over every usable, nonexcluded reference in each
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
