Notebooks
=========

Start with vX receptor annotation. Each lesson pairs code with its data
requirements and provides the notebook source. The reference and hit policy
are recorded in the annotation manifests.

Start here
----------

.. raw:: html

   <div class="proj-card-grid">
     <a class="proj-card" href="notebooks/vdjdb_annotation.html">
       <span class="proj-card__type">Official VDJdb release</span>
       <h3>Annotate an AIRR query against VDJdb</h3>
       <p>Use the current release or a pinned local ZIP. Read vX evidence,
       query dispositions, and the reference manifest.</p>
     </a>
     <a class="proj-card" href="notebooks/annotation.html">
       <span class="proj-card__type">Offline · custom reference</span>
       <h3>Annotate against your own reference</h3>
       <p>Run vX on a complete peptide–MHC table and look up opaque class
       labels with the direct labeled-reference API.</p>
     </a>
     <a class="proj-card" href="notebooks/comparison.html">
       <span class="proj-card__type">Offline · synthetic inputs</span>
       <h3>Compare exact-hit policies in vX</h3>
       <p>Keep the vX score and reference fixed while changing whether
       exact junction hits can contribute to evidence.</p>
     </a>
   </div>

* :download:`VDJdb annotation source <../examples/vdjdb_annotation.ipynb>`
* :download:`Custom reference source <../examples/annotation.ipynb>`
* :download:`Exact-hit policy source <../examples/comparison.ipynb>`
* :download:`Repertoire evidence report source <../examples/repertoire_evidence.ipynb>`:
  run the shipping report CLI on synthetic final-table fixtures, inspect all
  200 mock endpoints and verify preservation of recorded P and q values.
  This offline lesson performs no biological inference. For real final tables,
  follow :doc:`repertoire-method` and retain the producing manifest.
  Read the :doc:`rendered lesson <notebooks/repertoire_evidence>` alongside its code.

The offline lessons use the supplied ``examples/inputs`` files. See
:doc:`examples` for the shell commands and input downloads.

Biological questions
--------------------

.. raw:: html

   <div class="proj-card-grid">
     <a class="proj-card" href="notebooks/accuracy.html">
       <span class="proj-card__type">Metric tutorial · historical kernel</span>
       <h3>ROC, precision–recall and false positives</h3>
       <p>Compare ROC and precision–recall curves, score distributions and
       true/false-positive counts for 2,670 labelled observations.</p>
     </a>
     <a class="proj-card" href="notebooks/vaccination.html">
       <span class="proj-card__type">Deferred draft · vaccination repertoires</span>
       <h3>LLW evidence across donors</h3>
       <p>Retain all eligible clonotypes and technical replicates, then rank
       LLW evidence with the declared A*02 restriction.</p>
     </a>
     <a class="proj-card" href="notebooks/precursor_estimates.html">
       <span class="proj-card__type">Recombination model · complete target sets</span>
       <h3>From precursor frequency to expected cells</h3>
       <p>Compare GIL, NLV and LLW: probability mass, model frequency,
       expected cells and occupancy.</p>
     </a>
   </div>

* :download:`Accuracy source <../examples/accuracy.ipynb>`
* :download:`Vaccination source <../examples/vaccination.ipynb>`
* :download:`Precursor estimates source <../examples/precursor_estimates.ipynb>`

The historical accuracy lesson uses supplied positive/control labels to
illustrate metric construction. The vaccination workflow is an unexecuted
draft, and the precursor lesson calculates model estimates. See
:doc:`performance` for the scope of measured runtime, RAM and CPU use.

Precursor method and validation
-------------------------------

.. raw:: html

   <div class="proj-card-grid">
     <a class="proj-card" href="notebooks/precursor.html">
       <span class="proj-card__type">Research tutorial · reference data</span>
       <h3>Precursor frequency for an epitope</h3>
       <p>Explore generation probability, union mass, unseen junctions and
       the conversion from repertoire mass to cell counts.</p>
     </a>
     <a class="proj-card" href="compendium.html">
       <h3>Measured precursor frequencies</h3>
       <p>Read the validation data, model conventions and setup for the interactive compendium.</p>
     </a>
   </div>

* :download:`Precursor method source <notebooks/precursor.ipynb>`
* :download:`Interactive compendium source <notebooks/precursor_compendium.py>`

Run a notebook
--------------

Keep the complete `examples directory <https://github.com/antigenomics/vdjmatch/tree/master/examples>`_
together and open it with a Jupyter kernel in your isolated environment.

.. code-block:: sh

   python -m pip install -r examples/requirements-notebooks.txt

Run this installation from the repository root: it installs the checkout whose
CLI the notebooks use. Released packages need examples from the matching release
tag. The gallery checker creates the kernel automatically.

For the biological lessons, set these variables if your local mirrors are
elsewhere:

.. list-table:: Input locations
   :header-rows: 1

   * - Variable
     - Default
     - Used by
   * - ``AIRR_BENCHMARK``
     - ``~/hf/airr_benchmark``
     - Accuracy, precursor estimates and the pinned reference
   * - ``AIRR_YFV19``
     - ``~/hf/airr_yfv19``
     - Vaccination repertoire

Keep :download:`profile_cli.py <../examples/profile_cli.py>` beside the notebooks.
Biological calculations use four native threads and write fresh CLI tables;
reviewed figures and aggregate counts are embedded in these pages. Raw rows,
sample identifiers and complete execution logs remain private. Documentation
builds render saved outputs without running calculations.

Verify the gallery
~~~~~~~~~~~~~~~~~~

The :download:`dependency recipe <../examples/requirements-notebooks.txt>` and
:download:`checker <../examples/check_gallery.py>` provide a clean-kernel check:

.. code-block:: sh

   uv venv .venv-gallery
   uv pip install --python .venv-gallery/bin/python -r examples/requirements-notebooks.txt
   .venv-gallery/bin/python examples/check_gallery.py NEW_OUTPUT_DIRECTORY

The default check runs the offline shell workflows and the custom reference,
exact-hit policy, and repertoire-evidence notebooks. Add
``--vdjdb-reference /path/to/release.zip`` to verify the official VDJdb
workflow offline. Add ``--scientific`` to execute the biological lessons,
including the deferred vaccination draft, if you have the raw mirrors and
network access. The checker saves executed notebooks and cell timings in the
new directory and exits nonzero on failure. The vaccination calculation has a
1,800-second cell budget. Other cells have up to 300 seconds. Documentation
builds execute no cells.

.. toctree::
   :hidden:

   VDJdb annotation <notebooks/vdjdb_annotation>
   Custom reference annotation <notebooks/annotation>
   Exact-hit policy <notebooks/comparison>
   Repertoire evidence report <notebooks/repertoire_evidence>
   ROC and precision–recall <notebooks/accuracy>
   Vaccination repertoire <notebooks/vaccination>
   Precursor estimates <notebooks/precursor_estimates>
   Precursor method <notebooks/precursor>
