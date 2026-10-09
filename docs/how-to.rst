Work with releases and repertoires
==================================

Choose or pin a reference
-------------------------

.. code-block:: python

   from vdjmatch import Annotator

   latest = Annotator.latest(asset="default", cache="reference-inputs")
   pinned = Annotator.version("YOUR_RELEASE_TAG", cache="reference-inputs")
   local = Annotator.from_path("reference.zip", asset="default")
   extracted = Annotator.from_path("reference-directory", asset="primary")
   table = Annotator.from_path("reference.tsv")

Local paths resolve without downloading a release. Supported inputs include
legacy tables, release ZIPs and extracted directories, rich records/chains
bundles, and AIRR Rearrangement/Reactivity bundles. Bundle manifests and schema
metadata identify the appropriate assets and joins. Incomplete or incompatible
inputs raise an error rather than silently choosing a different reference.

``default`` selects the release's default supported representation. ``primary``,
``legacy`` and ``airr`` select bundle roles when provided; ``slim`` and ``full``
select the corresponding legacy projections. Request a representation that the
chosen release supplies. A new bundle supported by the resolver is not necessarily
available in the latest public release.

Reference rows with an invalid amino-acid junction are retained with
``reference_valid=False`` and excluded from the search index. The reference report
states retained/excluded counts. Inspect that report when assessing reference coverage.

An HF mirror is selected explicitly with ``source="hf"``; it is distinct from the
GitHub latest/pinned release source. Reuse an ``Annotator`` to keep one explicit
reference index across calls.

Filter a reference before indexing
----------------------------------

``db.load`` and ``Annotator.from_path`` accept ``epitope``, ``mhc_a``, ``mhc_b`` and
``reference_id`` as a string or list of strings. Values within a selector are
alternatives; different selectors must all match. ``exclude_reference_ids``
removes specified studies while retaining rows without a study identifier.

.. code-block:: python

   restricted = Annotator.from_path(
       "reference.zip", epitope=["GILGFVFTL", "NLVPMVATV"],
       mhc_a="HLA-A*02:01", exclude_reference_ids=["study-to-exclude"],
   )

``evidence_type`` selects records with any attached evidence observation of a
requested type, retaining every selected source chain and its complete evidence
metadata. It accepts the same string/list convention. For a rich release, for
example, ``evidence_type="structure_native"`` requires that type in the attached
evidence. A requested evidence predicate raises when the reference does not supply
an ``evidence_type`` field; it never silently ignores an unavailable predicate.
These selectors retain source order and do not deduplicate reference observations.

Read AIRR and custom tables
---------------------------

.. code-block:: python

   from vdjmatch.io import read_rearrangement

   queries, report = read_rearrangement(
       "repertoire.tsv.gz", source="airr", return_report=True,
   )
   calls = local.annotate(queries, threads=1)

AIRR ``junction_aa`` includes the conserved anchors and is preferred even when
``cdr3_aa`` is also present. AIRR ``cdr3_aa`` excludes them; a true AIRR table that
has only that field is rejected. The reader never fabricates anchors.

Legacy VDJdb ``cdr3`` and legacy producer ``cdr3_aa`` can already mean a junction.
For an ambiguous custom table, state the convention:

.. code-block:: python

   queries = read_rearrangement(
       "custom.csv.gz", source="custom", sequence_convention="junction",
   )

Both plain and gzipped CSV/TSV tables are supported. Standard repertoire formats
use published ``vdjtools`` readers. If a converter collapses or removes source rows
so that original identities cannot be preserved, ingestion raises explicitly.

By default, retained source rows keep their order, ``query_id``, supplied
``sequence_id``, V/J calls, counts and pairing linkage. Missing or invalid sequences
are counted in the optional report. A supplied count must be a nonnegative integer;
malformed, missing, negative and fractional counts raise an error.

To deliberately collapse matching keys, pass ``dedup=True``. The key includes
junction, V, J, locus and pair linkage; counts are summed, and ``query_ids`` and
``sequence_ids`` retain the contributing source identities. This is a different
counting unit from the original row table.

Pair chains by a cell identifier
--------------------------------

.. code-block:: python

   from vdjmatch.io import read_cell

   cells, report = read_cell(
       "cell-rearrangements.tsv", link="cell_id", return_report=True,
   )
   calls = local.annotate_paired(cells, cdr3a="cdr3a", cdr3b="cdr3b", threads=1)

``link`` can name a custom source column. Every retained chain must have a linkage
value and a TRA/TRB locus. Multiple retained alpha or beta rows for a cell raise by
default. ``ambiguity="first"`` explicitly selects the first row in source order.
Missing chains have typed null columns and a descriptive ingestion ``pair_status``.

Paired annotation requires both chains to match the same reference complex with
the same peptide–MHC identity. Agreement between two independent chain labels does
not establish this paired-reference evidence. Request paired calibration explicitly:

.. code-block:: python

   calls = local.annotate_paired(
       cells, cdr3a="cdr3a", cdr3b="cdr3b", calibrate=True,
       species="human", threads=1,
   )

For user-supplied paired controls, pass ``control={"TRA": alpha_index,
"TRB": beta_index}``; both indexes are required. The returned calibration label
``paired_independent_fixed_ball`` states the independent-chain background assumption.

Supply a background explicitly
------------------------------

.. code-block:: python

   from seqtree import Index

   controls = {"TRB": Index.build(["CWWWWF", "CYYYYF"], alphabet="aa")}
   candidates = local.candidates(queries, control=controls, threads=1)

This tiny background illustrates the interface; use an appropriate versioned
background for an application. The API uses uncalibrated ranking unless
``calibrate=True`` or a control index is supplied. A locus-to-index mapping allows
separate alpha and beta backgrounds. Empty backgrounds and requested background
loading failures raise errors.

Hard V/J matching can be requested for uncalibrated annotation with
``match_v=True`` or ``match_j=True``. Calibrated requests with these restrictions
are rejected because sequence-only control indexes cannot apply the same V/J
predicate.

Run the command line
--------------------

.. code-block:: bash

   vdjmatch update --asset default --cache reference-inputs
   vdjmatch match --vdjdb reference.zip --input-format airr --threads 1 \
       --no-evalue -o results/run repertoire.tsv
   vdjmatch match --vdjdb reference-directory --asset primary \
       --input-format custom --sequence-convention junction --no-evalue \
       -o results/custom custom.csv.gz
   vdjmatch match --vdjdb reference.zip --paired --link cell_id \
       --input-format airr -o results/cells cell-rearrangements.tsv

CLI matching requests background calibration by default; ``--no-evalue`` selects
uncalibrated ranking. The command writes detailed hits, all candidates, calls,
descriptive count summaries, ingestion counts and a manifest. See :doc:`cli` for
options and :doc:`explanation` before interpreting the summaries.
