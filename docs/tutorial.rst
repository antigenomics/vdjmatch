Annotate your first receptors
=============================

This tutorial uses a tiny reference supplied in memory. It needs no download and
shows how individual receptor rows become peptide–MHC candidates and calls.

Install
-------

.. code-block:: bash

   python -m pip install vdjmatch

The package installs ``seqtree`` for native batch search and ``vdjtools`` for
repertoire ingestion. Use a virtual environment for an isolated installation.

Create a reference and queries
------------------------------

.. code-block:: python

   import polars as pl
   from vdjmatch import Annotator

   reference = pl.DataFrame({
       "gene": ["TRB", "TRB"],
       "cdr3": ["CASSIRSSYEQYF", "CASSPGTGYEQFF"],
       "v": ["TRBV19", "TRBV5"],
       "j": ["TRBJ2-7", "TRBJ2-1"],
       "epitope": ["GILGFVFTL", "NLVPMVATV"],
       "mhc_a": ["HLA-A*02:01", "HLA-A*02:01"],
       "mhc_class": ["MHCI", "MHCI"],
       "reference_id": ["study-a", "study-b"],
   })
   queries = pl.DataFrame({
       "sequence_id": ["query-a", "query-b", "query-c"],
       "junction_aa": ["CASSIRSSYEQYF", "CASSIRSSYEQYF", "CWWWWF"],
       "v_call": ["TRBV19", "TRBV19", "TRBV1"],
       "j_call": ["TRBJ2-7", "TRBJ2-7", "TRBJ1-1"],
       "locus": ["TRB", "TRB", "TRB"],
   })
   annotator = Annotator.from_frame(reference)
   calls = annotator.annotate(queries, scope="0", threads=1)
   print(calls.select("sequence_id", "vdjmatch_epitope", "vdjmatch_status"))

The first two rows retain their separate identities and receive the same peptide
candidate. The third has ``no_hit`` and remains in the output. With no background
requested, these calls use uncalibrated ranking evidence.

Inspect alternatives
--------------------

.. code-block:: python

   candidates = annotator.candidates(queries, scope="1", threads=1)
   print(candidates.select(
       "query_id", "epitope", "mhc_a", "rank", "ned_score",
       "n_clonotypes", "n_studies", "n_competing_clonotypes",
   ))

Each candidate is a peptide together with its MHC restriction. Read all candidates
when a query has competing labels. Exact ranking ties produce ``ambiguous`` calls.
The score and its margin describe ranking; neither is a posterior probability.

Use a VDJdb release
-------------------

.. code-block:: python

   annotator = Annotator.latest(species="HomoSapiens")
   calls = annotator.annotate(queries, calibrate=True, species="human", threads=1)

``latest`` resolves the GitHub release by default. Fetching the reference or a
background that is not bundled requires network access. For an offline run,
follow :doc:`how-to` with an existing local reference and explicit background.

The calibrated output adds the expected background match count ``E`` and the
Poisson-tail ``p_enrichment`` for the stated candidate reference set. Their
interpretation and counting units are explained in :doc:`explanation`.
