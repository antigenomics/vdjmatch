Annotate single and paired receptors
====================================

Run one command, then rank candidate antigens with ``kernel_score``.
A larger score means stronger reference evidence. A match decision additionally
needs the matching calibration described in :doc:`explanation`.

Version 0.6.0 uses the curated manuscript reference population by default.
Pass ``--reference-population all-valid`` for the preserved all-valid baseline.
Exact evidence remains included in ordinary annotation.

Single-chain input
------------------

A tab-separated file contains one receptor per row. ``junction_aa`` includes
the conserved terminal residues; ``v_call`` and ``locus`` identify its V gene
and chain. These are the supplied synthetic example rows:

.. literalinclude:: ../examples/inputs/cohort/sample_a.tsv
   :language: text

Run
~~~

.. code-block:: sh

   vdjmatch annotate examples/inputs/cohort/sample_a.tsv \
     --vdjdb examples/inputs/reference.tsv --locus TRB \
     --threads 4 -o single

Read the output
~~~~~~~~~~~~~~~

``single.selection.tsv`` contains one row per input receptor. Selected columns
from the verified example are shown below; the full table preserves input IDs.

.. list-table:: Single-chain result
   :header-rows: 1

   * - sequence_id
     - selected_epitope
     - kernel_score
     - kernel_competing_score
     - selection_status
   * - q1
     - DEMO_A
     - 0.978264
     - 0.220985
     - ranked_evidence
   * - q2
     - DEMO_A
     - 0.978264
     - 0.219845
     - ranked_evidence

Both receptors rank DEMO_A above the alternative label. The score is evidence
for ranking, not a probability of binding. DEMO_A and DEMO_B are synthetic labels.
Use ``--details`` to export evidence for every candidate peptide–MHC.

Paired input
------------

Put linked alpha and beta chains on the same row, with a stable ``clone_id``.
Do not manufacture pairs from independent single-chain receptors.

.. literalinclude:: ../examples/inputs/paired.tsv
   :language: text

Run
~~~

.. code-block:: sh

   vdjmatch annotate examples/inputs/paired.tsv \
     --vdjdb examples/inputs/reference.tsv --locus paired \
     --threads 4 -o paired

Read the output
~~~~~~~~~~~~~~~

``paired.selection.tsv`` keeps linked evidence and both marginal chains:

.. list-table:: Paired result
   :header-rows: 1

   * - clone_id
     - Linked label / score
     - Alpha label / score
     - Beta label / score
     - kernel_any_axis_score
   * - pair1
     - DEMO_A / 0.962190
     - DEMO_A / 0.739132
     - DEMO_A / 0.978264
     - 0.978264
   * - pair2
     - unavailable
     - DEMO_B / 0.739132
     - unavailable
     - 0.739132

The complete pair supports DEMO_A on all three axes. The second row lacks a
beta chain: its alpha evidence remains available, but it has no linked-pair
score. ``kernel_any_axis_score`` is the maximum of the three selected scores
for the declared any-axis calibration event; it does not invent a joint label.

Parameters that change the question
-----------------------------------

.. list-table:: Main options
   :header-rows: 1

   * - Parameter
     - Meaning
   * - ``--vdjdb reference.tsv``
     - Fix the reference used for all comparisons. Omit to fetch the current release.
   * - ``--locus TRA|TRB|paired``
     - Choose the receptor population; automatic routing is the default.
   * - ``--epitope LLWNGPMAV``
     - Restrict candidate evidence to the requested peptide.
   * - ``--mhc-a 'HLA-A*02' --mhc-match compatible``
     - Restrict to compatible HLA-A*02 annotations at the supplied resolution.
   * - ``--exclude-exact``
     - Remove exact query evidence for a benchmark. Ordinary annotation keeps it.
   * - ``--details``
     - Export all candidate evidence as well as the selected label.
   * - ``--threads 4``
     - Use a budget of four CPU workers. This changes runtime, not the score.

Install and reproduce everything above
--------------------------------------

.. code-block:: sh

   python -m venv .venv
   . .venv/bin/activate
   python -m pip install vdjmatch
   git clone --depth 1 https://github.com/antigenomics/vdjmatch.git
   cd vdjmatch
   sh examples/annotate.sh annotation-demo

The script runs the supplied single-chain cohort and paired inputs from scratch.
An invalid V gene in the second cohort file produces ``no_support`` rather than
losing the row. For real data, start with these same commands and replace the
input and reference paths. See :doc:`explanation` for significance and
:doc:`algorithms` for the scoring steps.
