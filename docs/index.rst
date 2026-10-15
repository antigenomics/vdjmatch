vdjmatch
========

.. raw:: html

   <div class="proj-intro">
     <p class="proj-intro__eyebrow">TCR specificity and repertoire annotation</p>
     <p class="proj-intro__lead">From receptor sequences to peptide–MHC evidence.
     Annotate TRA/TRB receptors, linked pairs or whole repertoires while keeping
     source identities, competing labels and unsupported queries visible.</p>
     <p class="proj-intro__links">
       <a href="tutorial.html">Get started</a><span>·</span>
       <a href="cli.html">Command line</a><span>·</span>
       <a href="examples.html">Examples</a><span>·</span>
       <a href="notebooks.html">Notebooks</a>
     </p>
   </div>

Version |release|. Install the published package and annotate your AIRR table:

.. code-block:: sh

   python -m pip install vdjmatch
   vdjmatch annotate sample.tsv --threads 4 -o annotation

The reference is fetched on first use. Choose ``--vdjdb reference.tsv`` for an
offline run, or start with the tiny synthetic inputs in :doc:`examples`.

.. raw:: html

   <div class="proj-card-grid">
     <a class="proj-card" href="tutorial.html">
       <span class="proj-card__type">Tutorials</span>
       <h3>Annotate your first receptors</h3>
       <p>Install, run a small offline cohort and paired sample, then read the output tables.</p>
     </a>
     <a class="proj-card" href="how-to.html">
       <span class="proj-card__type">How-to guides</span>
       <h3>Run your own cohort</h3>
       <p>Choose a reference, prepare AIRR inputs, preserve pairing and set a resource budget.</p>
     </a>
     <a class="proj-card" href="cli.html">
       <span class="proj-card__type">Reference</span>
       <h3>Look up a command</h3>
       <p>Find flags, input contracts and output fields for annotation and calibration.</p>
     </a>
     <a class="proj-card" href="explanation.html">
       <span class="proj-card__type">Explanation</span>
       <h3>Understand the scores</h3>
       <p>Interpret reference support, competing labels, paired evidence and false-call calibration.</p>
     </a>
   </div>

Explore a worked example
------------------------

.. raw:: html

   <div class="proj-card-grid">
     <a class="proj-card" href="examples.html">
       <h3>Examples gallery</h3>
       <p>Two runnable shell workflows with raw inputs: cohort annotation and a geometry comparison.</p>
     </a>
     <a class="proj-card" href="notebooks.html">
       <h3>Notebooks gallery</h3>
       <p>Read worked notebooks with figures, results and downloadable source.</p>
     </a>
   </div>

Use the :doc:`api` to integrate annotation into Python, consult :doc:`releases`
for release changes, and read :doc:`scoring` for the detailed score derivation.
Native search uses `seqtree <https://github.com/antigenomics/seqtree>`_;
ingestion uses `vdjtools <https://github.com/antigenomics/vdjtools>`_.

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Tutorials

   Getting started <tutorial>

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: How-to guides

   Workflows <how-to>
   Examples <examples>
   Notebooks <notebooks>
   Laptop resources <performance>

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Reference

   Command line <cli>
   Python API <api>
   releases

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Explanation

   Evidence and calibration <explanation>
   Two control directions <control-calibration>
   Algorithm diagrams <algorithms>
   Score derivation <scoring>
   pgen
   compendium
