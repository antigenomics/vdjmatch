Examples gallery
================

Two offline workflows and three biological notebooks show the path from raw inputs to CLI tables. Start with
annotation, then compare detailed evidence. Open :doc:`notebooks` for rendered
previews and notebook source downloads.

.. raw:: html

   <div class="proj-card-grid">
     <a class="proj-card" href="#annotation-tutorial">
       <span class="proj-card__type">Tutorial · four rows and two pairs</span>
       <h3>Annotate a cohort</h3>
       <p>Run metadata and linked-pair inputs. Inspect identities, selections,
       invalid genes and the incomplete pair.</p>
     </a>
     <a class="proj-card" href="#geometry-comparison-how-to">
       <span class="proj-card__type">How-to · matched query identities</span>
       <h3>Compare geometries</h3>
       <p>Export fresh stem4 and vcdr12 evidence and join their full pMHC keys.</p>
     </a>
   </div>

These examples run tiny synthetic raw inputs against a supplied offline reference.
``DEMO_A`` and ``DEMO_B`` are demonstration labels, not experimental binding truth.
Every result is produced by the installed CLI; notebooks analyze only its final tables.
There are no persisted feature inputs, skip/resume paths or reused predictions.

`Download the release source ZIP <https://github.com/antigenomics/vdjmatch/archive/refs/tags/v0.4.0.zip>`_
or clone the repository and use its complete ``examples`` directory, keeping
its ``inputs`` directory beside the scripts and notebooks. The source distribution
also includes these files. Install vdjmatch in an isolated environment and activate
it so ``vdjmatch`` is on ``PATH``. Both scripts default to four native threads.

Annotation tutorial
-------------------

From the repository root, choose a new output directory:

.. code-block:: sh

   sh examples/annotate.sh annotation-demo

The metadata cohort has four single-chain rows across two samples. The wide paired
input has two receptors, including one missing beta chain. Source identities and
missing/invalid dispositions remain in ``*.queries.tsv`` and ``*.selection.tsv``.
``*.manifest.json`` records the reference, scoring and resources.

* :download:`Annotation shell script <../examples/annotate.sh>`
* :download:`Annotation notebook <../examples/annotation.ipynb>`
* :doc:`Read the annotation notebook <notebooks/annotation>`
* :download:`Offline synthetic reference <../examples/inputs/reference.tsv>`
* :download:`Sample metadata <../examples/inputs/samples.tsv>`
* :download:`First AIRR sample <../examples/inputs/cohort/sample_a.tsv>`
* :download:`Second AIRR sample <../examples/inputs/cohort/sample_b.tsv>`
* :download:`Wide paired sample <../examples/inputs/paired.tsv>`

Geometry comparison how-to
--------------------------

.. code-block:: sh

   sh examples/compare.sh comparison-demo

This runs default stem4 and the retained advanced ``vcdr12`` geometry, both with
``--exclude-exact --details`` and separate fresh output prefixes. The notebook
joins detailed ``*.enrichment.tsv`` rows by source identity and complete pMHC key,
then displays their continuous scores. Any arithmetic score delta is a geometry
diagnostic, not an accuracy gain or calibrated comparison of score scales.
It does not calculate scientific accuracy
from demonstration labels or transfer calibration between geometries.

* :download:`Comparison shell script <../examples/compare.sh>`
* :download:`Comparison notebook <../examples/comparison.ipynb>`
* :doc:`Read the comparison notebook <notebooks/comparison>`

Both scripts fail on nonzero CLI exit and reject an existing output directory.
For another fresh run, choose another directory. Notebooks use a new temporary
output directory each time and execute ``sys.executable -m vdjmatch.cli``
with ``check=True`` and
a 120-second timeout. Open them from the repository root or ``examples`` directory.
The notebook kernel and CLI therefore use the same Python environment and do not
require a POSIX shell. Their committed cells contain no execution outputs. Documentation builds never
execute notebooks.

For input formats and operational choices, see :doc:`how-to`; for the meaning of
scores and the distinction from binding probabilities, see :doc:`explanation`.
For all flags and library signatures, see :doc:`cli` and :doc:`api`.

Biological walkthroughs
-----------------------

The :doc:`notebooks` gallery adds longer ROC/PR, vaccination-repertoire and
precursor-frequency lessons. They use the local ``airr_benchmark`` and
``airr_yfv19`` raw-input mirrors, rebuild core calculations through the CLI,
and export publication-size PDF/SVG/PNG figures. The :doc:`performance` page
reports measured workload sizes, elapsed time, process-tree RAM and CPU use.
