Examples
========

Run two small offline workflows from raw inputs to CLI tables. For longer
lessons with rendered code, figures and results, continue to :doc:`notebooks`.

Setup
-----

Install vdjmatch in an isolated environment and obtain the complete
`examples directory <https://github.com/antigenomics/vdjmatch/tree/master/examples>`_,
keeping ``inputs`` beside the scripts. Activate the environment so
``vdjmatch`` is on ``PATH``. Both scripts use four native threads and a supplied
offline reference.

The inputs are synthetic software fixtures. ``DEMO_A`` and ``DEMO_B`` are
demonstration labels, rather than experimental binding labels.

Annotate a cohort
-----------------

From the repository root, choose a new output directory:

.. code-block:: sh

   sh examples/annotate.sh annotation-demo

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

Continue with biological data
-----------------------------

The :doc:`notebooks` gallery covers ROC/PR and true/false positives, baseline
vaccination-repertoire evidence, and precursor-frequency estimates. Its setup
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
