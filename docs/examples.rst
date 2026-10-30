Examples
========

Run three small offline workflows from raw inputs to CLI tables. For longer
lessons with rendered code, figures and results, continue to :doc:`notebooks`.

Setup
-----

Install vdjmatch in an isolated environment and obtain the complete
`examples directory <https://github.com/antigenomics/vdjmatch/tree/master/examples>`_,
keeping ``inputs`` beside the scripts. Activate the environment so
``vdjmatch`` is on ``PATH``. The scripts use four native threads and a supplied
offline reference.

The inputs are synthetic software fixtures. ``DEMO_A`` and ``DEMO_B`` are
demonstration labels, rather than experimental binding labels.

Annotate a cohort
-----------------

From the repository root, choose a new output directory:

.. code-block:: sh

   sh examples/annotate.sh annotation-demo

To use an existing VDJdb build, provide its primary table as the second argument:

.. code-block:: sh

   sh examples/annotate.sh build-demo /path/to/build/tables/vdjdb.parquet

A legacy table is also accepted. AIRR bundles need receptor species metadata
for species-filtered annotation; VDJdb AIRR exports without that field should
use the primary or legacy reference instead. This only annotates
the supplied synthetic queries; it does not rebuild the database or measure
biological accuracy. The manifest records the chosen reference identity.

The cohort contains four single-chain rows across two samples and two paired
receptors, including an incomplete pair. Inspect these outputs:

.. list-table:: Annotation outputs
   :header-rows: 1

   * - File
     - What it shows
   * - ``*.queries.tsv``
     - Source identities and ingestion dispositions
   * - ``*.selection.tsv``
     - Selected evidence, ambiguity and unsupported queries
   * - ``*.manifest.json``
     - Reference identity, scoring options and resources

* :doc:`Read the rendered notebook <notebooks/annotation>`
* :download:`Notebook source <../examples/annotation.ipynb>` · :download:`Shell script <../examples/annotate.sh>`
* :download:`Offline reference <../examples/inputs/reference.tsv>` · :download:`Sample metadata <../examples/inputs/samples.tsv>`
* :download:`First AIRR sample <../examples/inputs/cohort/sample_a.tsv>` · :download:`Second AIRR sample <../examples/inputs/cohort/sample_b.tsv>`
* :download:`Paired sample <../examples/inputs/paired.tsv>`

Compare geometries
------------------

.. code-block:: sh

   sh examples/compare.sh comparison-demo

This calculates fresh, exact-excluded evidence for default stem4 and the
advanced ``vcdr12`` geometry. The notebook joins ``*.enrichment.tsv`` by source
identity and full peptide–MHC key, then displays the score differences.
These are geometry diagnostics: synthetic labels cannot establish biological
accuracy, and score scales do not share calibration.

* :doc:`Read the rendered notebook <notebooks/comparison>`
* :download:`Notebook source <../examples/comparison.ipynb>` · :download:`Shell script <../examples/compare.sh>`

Check database presence
-----------------------

.. code-block:: sh

   sh examples/database.sh database-demo

The supplied three-query table and one synthetic raw control produce these
columns in ``database-demo/presence.database.tsv``:

.. list-table:: Exact CLI output for the synthetic diagnostic
   :header-rows: 1

   * - Query
     - Exact junction
     - Nearest edits
     - Both-direction P
     - Status
   * - ``exact``
     - true
     - 0
     - 1.0
     - matched
   * - ``near``
     - false
     - 1
     - 1.0
     - matched
   * - ``unsupported``
     - false
     - unavailable
     - 1.0
     - no_hit

Pgen is recalculated from the bundled model. This deliberately small control
places its rare query strata near the synthetic reference; exact and near
matches are therefore common under its null, and both directional tails are
one. The unsupported query has both tails one and a zero estimated null match
rate. The example demonstrates source ingestion, stratification, retained exact
overlaps and the distinction between presence, a score-tail P-value and an
any-match rate. The fixture does not validate a biological null.

* :download:`Shell script <../examples/database.sh>`
* :download:`Queries <../examples/inputs/database_queries.tsv>` · :download:`Single synthetic control <../examples/inputs/database_control.tsv>`
* :doc:`Interpret the two directions <control-calibration>`

Continue with biological data
-----------------------------

Compare three retained scores from the terminal:

.. code-block:: sh

   python examples/accuracy.py --benchmark /path/to/airr_benchmark \
     --reference reference.zip --germline-background raw-prior.airr.tsv \
     --sorted-input original-sorted.tsv --threads 4 --output accuracy-run

The raw prior supplies junction/V/J observations for original density scoring.
The optional original sorted table supplies LLW and LLL population labels in
``source_epitope``, ``epitope`` or ``antigen.epitope`` and its original gene calls.
Omitting it uses the mirror's sorted observations and labels that contrast
explicitly. The complete mirror panels remain separate positive/control tasks.

``accuracy.tsv`` gives counts, ROC-AUC and average precision for original,
gap-block and stem4 evidence. ``observations.tsv`` retains all three scores;
``rank-cases.tsv`` identifies individual ranking errors, and
``operating-counts.tsv`` gives TP/FP/FN/TN at observed score thresholds.
Each case's rank loss is the fraction of opposite-class rows that rank in the
wrong order, counting a tie as one half. Its task mean equals 1 minus ROC-AUC;
larger values identify worse rankings. Source labels stay unchanged.
``resources.tsv`` records each complete CLI stage's runtime, sampled process-tree
RAM, occupied cores and live threads. SVG/PDF/PNG figures are exported alongside
the tables. Every invocation rebuilds from raw inputs in a new output directory.
These comparisons rank declared assay targets; they do not calibrate a call.

Score all locally available vaccination repertoires with their metadata:

.. code-block:: sh

   python examples/vaccination.py --dataset /path/to/airr_yfv19 \
     --reference reference.zip --threads 4 --memory-gb 24 --output vaccination-run

This retains every eligible TRB junction/V/J clonotype and technical replicate,
using the author-declared A*02 restriction. ``coverage.tsv`` preserves complete
read denominators and exclusion counts. Valid-amino-acid junctions with
noncanonical anchor residues remain present and flagged; functionality is not
supplied by these inputs. ``observations.tsv`` joins fresh CLI
scores back to each replicate; ``support.tsv`` counts radius/exact-junction
support. Neither support flag is a binding call. ``evidence-coverage.tsv`` ranks
by LLW evidence with all score ties entering together, rather than selecting
the most abundant clonotypes. Raw identities and detailed outputs stay private.

The :doc:`notebooks` gallery covers ROC/PR and true/false positives,
vaccination-repertoire neighbourhood evidence, and precursor-frequency estimates. Its setup
section lists the required local mirrors. See :doc:`performance` for measured
workload sizes, runtime and memory.

Execution notes
---------------

Every run rebuilds CLI tables in a fresh directory. The scripts reject an
existing output directory and propagate CLI failures. Notebooks invoke the
CLI with their kernel's interpreter and a bounded timeout; documentation
builds render saved outputs without executing them.

Use :doc:`how-to` for your own formats, :doc:`explanation` to interpret evidence,
and :doc:`cli`/:doc:`api` for command and library details.
