API reference
=============

See :doc:`tutorial` for a runnable example, :doc:`how-to` for source conventions
and offline inputs, and :doc:`explanation` for score definitions.

The public API is re-exported from the top-level :mod:`vdjmatch` package —
:class:`vdjmatch.api.Annotator` and the module-level :func:`vdjmatch.annotate` shortcut. The
submodules below expose the building blocks (VDJdb access, the search index, E-values, I/O).

vdjmatch.api
------------

.. automodule:: vdjmatch.api
   :members:
   :undoc-members:
   :show-inheritance:

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
