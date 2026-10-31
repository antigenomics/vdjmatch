API reference
=============

See :doc:`tutorial` for a runnable example, :doc:`how-to` for source conventions
and offline inputs, and :doc:`explanation` for score definitions.

The public API is re-exported from the top-level :mod:`vdjmatch` package —
:class:`vdjmatch.api.Annotator` and the module-level :func:`vdjmatch.annotate` shortcut. The
submodules below expose the building blocks (VDJdb access, the search index, E-values, I/O).

Saved reference assets
----------------------

``Annotator.from_path(path, threads=4)`` prepares or reopens an adjacent native
reference bundle. ``reference_profiles`` exposes its flat pooled-epitope/full-pMHC
class table; ``reference_index_manifest`` records the bundle path and identity.
Source/filter changes invalidate reuse. ``Annotator.latest`` uses the same path.
``Annotator.from_frame`` keeps anonymous references in memory.
See :doc:`cli` for the descriptive decision rules and the ``vdjmatch index`` command.

Batch repertoire annotation
---------------------------

``annotate_repertoires`` accepts a table, directory, list of paths or a metadata
TSV with ``path`` and optional ``sample_id``. It preserves sample/source identities
and batches single or linked populations. Compact kernel selection is the default;
``selection_only=False`` also exports full candidate/per-target evidence.
``threads`` also bounds the Polars pool. Set ``POLARS_MAX_THREADS`` before importing
Polars if that library has already initialized its pool; a larger existing pool
raises an error before input loading. Generated-field collisions retain input
metadata under ``input_`` prefixes.

.. autofunction:: vdjmatch.annotate_repertoires

On-demand control calibration
-----------------------------

See :doc:`control-calibration` for the two null hypotheses, event units and
complete-score selection requirements. These helpers use transient tables;
the user supplies one raw repertoire, with no calibration-file lifecycle.

.. automodule:: vdjmatch.evalue.mock_reference
   :members:

.. automodule:: vdjmatch.evalue.directional
   :members:

vdjmatch.api
------------

The original ``Annotator`` API's automatic calibration derives the background
organism from the queried reference loci; mixed species or an explicit mismatch
raises an error, and missing metadata requires ``species="human"`` or
``species="mouse"``. Supplied control mappings must cover every active locus
(both TRA and TRB for paired queries); their organism is caller-owned.

.. automodule:: vdjmatch.api
   :members:
   :undoc-members:
   :show-inheritance:

Bounded native dependency
~~~~~~~~~~~~~~~~~~~~~~~~~

Candidate-only single-chain reduction requires the native
``seqtree.Index.search_batch_consume`` method, supplied by the declared
``seqtree>=1.3.2`` dependency. The API also checks the installed native method
and raises an explicit error when it is absent. Requests that retain full hit
tables continue to use the existing batch API.

Released-dependency CI verifies the native single/paired consumers, empty-axis
support and hit-list-only exact exclusion in the installed environment.
The bounded candidate remains an explicit experimental scoring path; the
dependency release does not promote it to the default.

Detailed paired output
~~~~~~~~~~~~~~~~~~~~~~

``paired_candidates(..., return_hits=True)`` returns ``(detail, candidates)`` from
one reference search per chain. Detailed rows retain shared query/complex/pMHC keys
and chain metadata under ``alpha_`` and ``beta_`` prefixes. ``align=True`` adds each
chain's alignment/CIGAR fields. ``score_scale`` must be finite and positive; its
default is 400 for the bundled scoring convention. Both arguments are also supported
by ``annotate_paired``.

.. code-block:: python

   detail, candidates = annotator.paired_candidates(
       cells, cdr3a="cdr3a", cdr3b="cdr3b", align=True,
       score_scale=400.0, return_hits=True, threads=1,
   )

vdjmatch.db
-----------

.. automodule:: vdjmatch.db
   :members:
   :undoc-members:
   :show-inheritance:

vdjmatch.match
--------------

.. automodule:: vdjmatch.match
   :members:
   :undoc-members:
   :show-inheritance:

vdjmatch.evalue
---------------

.. automodule:: vdjmatch.evalue
   :members:
   :undoc-members:
   :show-inheritance:

vdjmatch.aggregate
------------------

.. automodule:: vdjmatch.aggregate
   :members:
   :undoc-members:
   :show-inheritance:

vdjmatch.io
-----------

.. automodule:: vdjmatch.io
   :members:
   :undoc-members:
   :show-inheritance:

vdjmatch.cluster
----------------

.. automodule:: vdjmatch.cluster
   :members:
   :undoc-members:
   :show-inheritance:

vdjmatch.precursor
------------------

.. note::

   Precursor estimation is a separate API. Its optional dependencies are available
   with ``pip install 'vdjmatch[precursor]'``; annotation does not run precursor estimation.

.. automodule:: vdjmatch.precursor
   :members:
   :undoc-members:
   :show-inheritance:

Candidate evidence
------------------

.. automodule:: vdjmatch.aggregate.candidates
   :members:
   :undoc-members:
   :show-inheritance:

Query normalization
-------------------

.. automodule:: vdjmatch.io.columns
   :members:
   :undoc-members:
   :show-inheritance:

Legacy paired first-hit helpers
-------------------------------

These helpers retain their adaptive first-hit contract separately from the
Annotator's fixed-ball paired-reference interface. ``paired_scan`` returns
``(hits, alpha_control_costs, beta_control_costs)``. By default each hit is the
historical two-tuple ``(radius, epitope)`` and ``pvalue`` counts observation
associations; its target size ``N`` must use that same unit.

With ``include_identity=True``, a hit is the three-tuple
``(radius, epitope, (alpha_junction, beta_junction))``. ``pvalue`` then counts
distinct junction pairs at the nearest paired radius; pass a distinct-pair target
size ``N``. When selecting an epitope, compute ``N`` for that epitope's reference
set. Mixing two-tuple and three-tuple hits is rejected.

.. automodule:: vdjmatch.evalue.paired
   :members: build_paired_ref, paired_scan, pvalue
   :undoc-members:
   :show-inheritance:

Experimental local-rank evidence
--------------------------------

.. autofunction:: vdjmatch.evalue.local_rank.local_rank_evidence

Curated source replication
--------------------------

:func:`vdjmatch.db.replicated` retains source provenance and complete pMHC
identity when selecting corroborated records.

Experimental V-conditioned motif evidence
-----------------------------------------

.. autofunction:: vdjmatch.match.motifs.masked_motif_scores

Targeted kernel evidence
------------------------

.. automodule:: vdjmatch.evalue.kernel
   :members: ladder_cutoffs, kernel_ladder, select_kernel

.. autofunction:: vdjmatch.match.geometry.precompute_v_loop_distances

Legacy calculation reproduction
-------------------------------

.. automodule:: vdjmatch.match.scoring_audit
   :members: matrix_audit, loo_audit

Reference graph morphology
--------------------------

.. autofunction:: vdjmatch.aggregate.structure.graph_profiles

.. autofunction:: vdjmatch.aggregate.structure.gene_selection_information

.. autofunction:: vdjmatch.aggregate.motif_information.motif_population

.. autofunction:: vdjmatch.aggregate.motif_information.motif_selection_information

.. autofunction:: vdjmatch.aggregate.motif_information.motif_bootstrap

.. autofunction:: vdjmatch.aggregate.motif_position.position_contrast

.. autofunction:: vdjmatch.aggregate.motif_position.generative_position_contrast

.. autofunction:: vdjmatch.aggregate.length_profile.length_histogram

.. autofunction:: vdjmatch.aggregate.length_profile.length_control

.. autofunction:: vdjmatch.aggregate.length_profile.length_summary

.. autofunction:: vdjmatch.aggregate.reference_coverage.reference_graph

.. autofunction:: vdjmatch.aggregate.reference_coverage.leiden_taxa

.. autofunction:: vdjmatch.aggregate.reference_coverage.abundance_coverage

.. autofunction:: vdjmatch.aggregate.germline_information.germline_projection_map

.. autofunction:: vdjmatch.aggregate.germline_information.iter_germline_information

.. autofunction:: vdjmatch.aggregate.paired_labels.paired_labels

Reference-null repertoire inference
-----------------------------------

See :func:`vdjmatch.api.repertoire_enrichment` above.

.. autofunction:: vdjmatch.aggregate.repertoire_enrichment.reference_enrichment

See :func:`vdjmatch.evalue.mock_reference.iter_epitope_mocks` above.

TCRnet-filtered repertoire counts
---------------------------------

.. autoclass:: vdjmatch.aggregate.repertoire_tcrnet.TcrnetMatcher
   :members: match

The owned reference and control indexes may be shared across samples within a
fresh cohort invocation. Filtering evidence is available as
``reference_filtering``. ``match`` returns summary counts, real hit links and
sample-filter evidence; it does not return repertoire P-values.


Bounded query ingestion
-----------------------

``vdjmatch.io.iter_rearrangement``

Use one complete sample at a time when its measured footprint is modest.
``iter_rearrangement`` provides within-sample batches for unusually large files;
it does not perform chunk-local gene imputation or deduplication. Consumers must
combine repeated clonotypes and paired identities across batch boundaries.
TCRnet requires sample-wide neighbourhood information.


Exact bulk reference-panel primitives
-------------------------------------

The ``repertoire-exact`` CLI uses checked model-specific selection
factors from vdjtools by default; ``--legacy`` preserves the compiled
uncalibrated convention. It composes these normalized-table primitives:

:func:`vdjmatch.evalue.mock_reference.compiled_reference_strata` constructs the probability strata and V/J audit.
Its default count-table path preserves the compiled Q=1 convention. The optional
``model_vj_probabilities`` argument instead accepts exact selected ``v,j,p``
marginals summing to one, with positive support, and replaces model counts.
The CLI ``repertoire-exact --selection-profile`` applies the coarse shape to joint mass
before the AA maximum and passes these selected marginals to the residual
within-V composition correction. Modern CLI defaults use calibrated global Q;
legacy reproduction is explicit.

:func:`vdjmatch.evalue.mock_reference.sample_stratified_panels` draws fresh panels.

.. autofunction:: vdjmatch.aggregate.repertoire_enrichment.exact_reference_counts
