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
       these predicates, so combine with ``--no-evalue``.
   * - ``--no-evalue``
     - Use explicitly uncalibrated ranking without requesting a background.
   * - ``--no-align``
     - Omit per-hit alignment and CIGAR computation.
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
