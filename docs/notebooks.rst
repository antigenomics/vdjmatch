Notebooks
=========

Read complete notebooks here in the documentation: prose, code and reviewed
figures appear together. Each page also offers its ``.ipynb`` source.
Start with the offline examples, then choose a biological question.

Start here
----------

.. raw:: html

   <div class="proj-card-grid">
     <a class="proj-card" href="notebooks/annotation.html">
       <span class="proj-card__type">Offline · synthetic inputs</span>
       <h3>Raw repertoires to annotation tables</h3>
       <p>Four single-chain rows and two pairs. Follow source identities,
       competing labels and the incomplete pair through the exported tables.</p>
     </a>
     <a class="proj-card" href="notebooks/comparison.html">
       <span class="proj-card__type">Offline · synthetic inputs</span>
       <h3>Compare two CLI geometries</h3>
       <p>Run stem4 and vcdr12 with exact exclusion, join detailed pMHC evidence,
       and inspect score differences on matched query identities.</p>
     </a>
   </div>

* :download:`Cohort annotation source <../examples/annotation.ipynb>`
* :download:`Geometry comparison source <../examples/comparison.ipynb>`

These small synthetic examples need only the supplied ``examples/inputs``
files. See :doc:`examples` for the shell commands and input downloads.

Biological questions
--------------------

.. raw:: html

   <div class="proj-card-grid">
     <a class="proj-card" href="notebooks/accuracy.html">
       <span class="proj-card__type">Accuracy · four TRB targets</span>
       <h3>ROC, precision–recall and false positives</h3>
       <p>Compare ROC and precision–recall curves, score distributions and
       true/false-positive counts for 2,670 labelled observations.</p>
     </a>
     <a class="proj-card" href="notebooks/vaccination.html">
       <span class="proj-card__type">Baseline repertoire · abundance</span>
       <h3>Candidate evidence in a vaccination dataset</h3>
       <p>Explore baseline clone abundance, read-count coverage and
       LLW candidate evidence in a vaccination repertoire.</p>
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

The accuracy lesson uses supplied positive/control labels. The vaccination
lesson describes a baseline repertoire; the precursor lesson calculates model
estimates. Each notebook states its observation unit, assumptions and limits.
See :doc:`performance` for measured runtime, RAM and CPU use.

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

   python -m pip install 'vdjmatch[bench,precursor,notebook]' ipykernel

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
   * - ``YFV_INPUT``
     - First table in sorted filename order
     - Choose a specific baseline repertoire

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

Add ``--scientific`` to include the three biological lessons once their raw
mirrors are available. The checker saves executed notebooks and cell timings
in the new directory, and exits nonzero if a workflow fails.

.. toctree::
   :hidden:

   Cohort annotation <notebooks/annotation>
   Geometry comparison <notebooks/comparison>
   ROC and precision–recall <notebooks/accuracy>
   Vaccination repertoire <notebooks/vaccination>
   Precursor estimates <notebooks/precursor_estimates>
   Precursor method <notebooks/precursor>
