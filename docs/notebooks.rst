Notebooks gallery
=================

Browse the notebooks before downloading them. Annotation and comparison use the
small offline inputs from :doc:`examples`; all tables come from the installed CLI.
Committed notebooks have no execution outputs, and documentation builds never
execute them.

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

.. toctree::
   :hidden:

   notebooks/annotation
   notebooks/comparison
   notebooks/precursor
