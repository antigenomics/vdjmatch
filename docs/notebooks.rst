Notebooks gallery
=================

Browse the notebooks before downloading them. Annotation and comparison use the
small offline inputs from :doc:`examples`; all tables come from the installed CLI.
The annotation and comparison notebooks have no committed execution outputs.
The precursor tutorial retains its curated illustration outputs; documentation
builds never execute notebooks.

Scientific worked examples
--------------------------

These longer lessons use biological data from local mirrors. Each core
calculation runs through the installed CLI; plots consume its final tables.
Reviewed plots are embedded for browsing. Complete raw rows, donor identifiers
and detailed execution logs remain private. See :doc:`performance` for measured
laptop time, RAM and CPU use.

.. raw:: html

   <div class="proj-card-grid">
     <a class="proj-card" href="notebooks/accuracy.html">
       <span class="proj-card__type">Experimental labels · four TRB targets</span>
       <h3>ROC, precision–recall and false positives</h3>
       <p>Follow 2,670 labelled source observations through exact exclusion,
       score distributions and true/false-positive count trajectories.</p>
     </a>
     <a class="proj-card" href="notebooks/vaccination.html">
       <span class="proj-card__type">Baseline repertoire · abundance</span>
       <h3>Candidate evidence in a vaccination dataset</h3>
       <p>Connect clone abundance, retained read-count coverage and relative
       LLW support without inventing binding labels or expansion.</p>
     </a>
     <a class="proj-card" href="notebooks/precursor_estimates.html">
       <span class="proj-card__type">Recombination model · complete target sets</span>
       <h3>From precursor frequency to expected cells</h3>
       <p>Compare GIL, NLV and LLW while separating probability mass,
       selection assumptions, compartment size and occupancy.</p>
     </a>
   </div>

* :download:`Accuracy notebook <../examples/accuracy.ipynb>`
* :download:`Vaccination notebook <../examples/vaccination.ipynb>`
* :download:`Precursor estimates notebook <../examples/precursor_estimates.ipynb>`
* :download:`CLI process profiler <../examples/profile_cli.py>`

Keep the profiler beside the notebooks. Install
``vdjmatch[bench,precursor,notebook]`` plus ``ipykernel``. Set
``AIRR_BENCHMARK`` and ``AIRR_YFV19`` to your local mirrors; defaults are
``~/hf/airr_benchmark`` and ``~/hf/airr_yfv19``. ``YFV_INPUT`` selects a
specific raw baseline repertoire. Every execution creates a fresh output directory.
The lessons use four native threads and one-thread numerical helper pools.

Annotation and comparison
-------------------------

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

* :download:`Download annotation.ipynb <../examples/annotation.ipynb>`
* :download:`Download comparison.ipynb <../examples/comparison.ipynb>`
* `Get the complete examples directory <https://github.com/antigenomics/vdjmatch/tree/master/examples>`_

Keep the notebooks beside ``inputs`` in ``examples``. Open them from the
repository root or ``examples`` directory, with a kernel that has vdjmatch
installed. Each run creates a fresh temporary output directory and invokes
``sys.executable -m vdjmatch.cli`` with a 120-second timeout. See :doc:`examples`
for raw input downloads and shell commands.

Precursor frequency
-------------------

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

* :download:`Download precursor.ipynb <notebooks/precursor.ipynb>`
* :download:`Download the marimo compendium <notebooks/precursor_compendium.py>`

The precursor workflow is separate from antigen annotation. Its setup and
reference-data requirements are described in :doc:`compendium`.

Run the complete gallery
------------------------

Use an isolated environment containing the published package and notebook tools.
From the repository root:

.. code-block:: sh

   uv venv .venv-gallery
   uv pip install --python .venv-gallery/bin/python --index-url https://pypi.org/simple --no-cache -r examples/requirements-notebooks.txt
   .venv-gallery/bin/python examples/check_gallery.py NEW_GALLERY_OUTPUT_DIRECTORY

This executes both shell examples, the original three Jupyter notebooks with fresh kernels,
and the marimo compendium. It records results and executed outputs in the new
directory and returns a failure status if any workflow fails. The precursor
tutorial uses four native threads and fetches its pinned public Hugging Face
snapshot; the compendium uses one native thread and fetches
its documented tables and reference. Documentation builds still do not execute
notebooks.

Add ``--scientific`` to execute the three biological notebooks too, once their
local raw-input mirrors are available. The checker keeps per-cell timestamps
and all complete execution records in the fresh output directory.

* :download:`Notebook dependency recipe <../examples/requirements-notebooks.txt>`
* :download:`Gallery verification script <../examples/check_gallery.py>`

.. toctree::
   :hidden:

   notebooks/annotation
   notebooks/comparison
   notebooks/precursor
   notebooks/accuracy
   notebooks/vaccination
   notebooks/precursor_estimates
