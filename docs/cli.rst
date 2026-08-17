Command-line interface
======================

Installing ``vdjmatch`` provides a single console script with two subcommands.

vdjmatch update
---------------

Fetch and cache a VDJdb release (re-used by ``match`` and the Python API).

.. code-block:: bash

   vdjmatch update [--asset {slim,full,default}] [--pin TAG] [--force]

==============  =============================================================
Option          Meaning
==============  =============================================================
``--asset``     which VDJdb table to fetch (default ``slim``)
``--pin``       pin a specific release tag (default: latest)
``--force``     re-download even if already cached
==============  =============================================================

vdjmatch match
--------------

Annotate one or more AIRR rearrangement samples against VDJdb, reporting ranked
hits, control-calibrated E-values and an epitope-enrichment summary.

.. code-block:: bash

   vdjmatch match [options] SAMPLE [SAMPLE ...]

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Option
     - Meaning
   * - ``-o, --output-prefix``
     - output path prefix (default ``vdjmatch_out``)
   * - ``--vdjdb``
     - custom VDJdb table path (default: fetch latest)
   * - ``--asset``
     - VDJdb table to match against (default ``full``)
   * - ``--pin``
     - pin a specific VDJdb release tag
   * - ``--species``
     - species filter (default ``HomoSapiens``)
   * - ``--scope``
     - search budget ``subs,ins,dels,total`` (default ``1,0,0,1``)
   * - ``--matrix``
     - ``vdjam`` (TCR-specific, bundled) or ``none`` (unit cost)
   * - ``--min-score``
     - minimum VDJdb confidence score (default ``0``)
   * - ``--match-v`` / ``--match-j``
     - require the V / J gene to match as well as the CDR3
   * - ``--no-evalue``
     - skip the control-calibrated E-value
   * - ``--no-align``
     - skip the per-hit CIGAR / alignment output
   * - ``--threads``
     - worker threads (``0`` = all cores)

Output
~~~~~~

``match`` writes three tab-separated tables per sample, prefixed with ``-o``:

- ``<prefix>.<sample>.hits.txt`` — every query→VDJdb hit with CDR3 alignment, CIGAR, edit counts and score.
- ``<prefix>.<sample>.calls.txt`` — one predicted epitope per query clonotype with its E-value.
- ``<prefix>.<sample>.summary.txt`` — epitope-level enrichment (unique clonotypes, reads) by MHC class.

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

.. warning::

   Sequences must be **junctions** (Cys104…Phe/Trp118 inclusive), not IMGT CDR3s — VDJdb's column
   is named ``cdr3`` but holds junctions. An anchor-stripped CDR3 scores exactly ``0.0`` with no
   error, so it is dropped and counted in ``n_dropped`` rather than silently scored.

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
     - independent rearrangements, for ``P(>=k precursors)``
   * - ``--min-junctions``
     - skip groups with fewer distinct junctions

Output
~~~~~~

One tab-separated table, one row per group. Key columns:

- ``pgen_log10_span``, ``top1_share``, ``top10_share`` — how far the cognate set's generation
  probabilities spread, and how concentrated the sum is. A mean is the wrong summary here.
- ``naive_sum``, ``union``, ``overlap`` — the per-sequence ball masses added up, the exact union,
  and the share of the naive sum that double-counting would have invented.
- ``retained``, ``F``, ``cells`` — the shell-weighted estimate, the calibrated frequency and the
  expected precursor count at ``--n-cells``.
- ``lambda``, ``p_ge_1``, ``p_ge_10`` — Poisson precursor probabilities (needs ``--n-eff``).
- ``n_ball``, ``n_unseen_ball``, ``unseen_ball_mass`` — the finite census: candidate cognate
  junctions in the neighbourhood that no database has catalogued, and the mass they carry.
- ``n_unseen_ht``, ``unseen_ht_mass``, ``rarity_ratio``, ``richness_reliable`` — the
  Horvitz–Thompson extrapolation. Read the **mass**; the count diverges in the tail and
  ``richness_reliable`` says when that is happening.
