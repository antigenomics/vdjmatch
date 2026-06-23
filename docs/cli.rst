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
