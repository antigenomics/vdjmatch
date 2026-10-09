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
