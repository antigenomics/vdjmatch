Annotate your first receptors
=============================

This tutorial starts with raw synthetic repertoire tables, invokes the installed
CLI, and reads its annotation outputs. It uses the pinned offline reference in the
``examples`` directory; no reference download is needed.

Install and run
---------------

Create and activate an isolated environment, install vdjmatch, and obtain the
complete repository ``examples`` directory, including its raw ``inputs`` files.

.. code-block:: sh

   python -m venv .venv
   . .venv/bin/activate
   python -m pip install vdjmatch
   sh examples/annotate.sh annotation-demo

The script submits the metadata cohort and the wide paired population as separate
CLI invocations with ``--threads 1`` and the supplied local reference. Four
single-chain rows retain their sample/source IDs; two paired receptors retain
clone identities, including the incomplete pair. ``DEMO_A`` and ``DEMO_B`` are
synthetic software labels, not claims of antigen binding.

Inspect the CLI tables
----------------------

.. code-block:: python

   import polars as pl

   queries = pl.read_csv("annotation-demo/cohort.queries.tsv", separator="\t")
   selection = pl.read_csv("annotation-demo/cohort.selection.tsv", separator="\t")
   print(queries.select("sample_id", "source_sequence_id"))
   print(selection)

Also inspect ``paired.queries.tsv``, ``paired.selection.tsv`` and both manifests.
Missing chains, invalid input and no-support dispositions stay in the output;
never interpret them as experimental negative labels. Ordinary annotation keeps
exact reference evidence. Relative kernel scores are neither binding probabilities
nor calibrated false-call rates.

Continue with the downloadable :doc:`examples` notebooks, which execute these
CLI commands from a fresh kernel and analyze their final tables. Follow
:doc:`how-to` for your own formats and resource budget; use :doc:`explanation` for
score interpretation and :doc:`cli`/:doc:`api` for reference documentation.
