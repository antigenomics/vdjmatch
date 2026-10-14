CLI examples and notebooks
==========================

These examples run tiny synthetic raw inputs against a supplied offline reference.
``DEMO_A`` and ``DEMO_B`` are demonstration labels, not experimental binding truth.
Every result is produced by the installed CLI; notebooks analyze only its final tables.
There are no persisted feature inputs, skip/resume paths or reused predictions.

Download or clone the complete ``examples`` directory from the repository, keeping
its ``inputs`` directory beside the scripts and notebooks. The source distribution
also includes these files. Install vdjmatch in an isolated environment and activate
it so ``vdjmatch`` is on ``PATH``. Both scripts default to one native thread.

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
