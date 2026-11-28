Validation data: measured precursor frequencies
===============================================

The precursor estimator ``F(e)`` describes the fraction of a naive repertoire
that recognizes epitope ``e``. The compendium collects published measurements
from tetramer and dextramer enrichment of unexposed donors, umbilical cord blood
and naive mice. Compare receptor annotation scores and null tails using the
separate :doc:`vx-formalism` workflow.

Two loaders fetch pinned measurement and prediction tables.

.. list-table::
   :header-rows: 1
   :widths: 22 16 62

   * - loader
     - provenance
     - what it is
   * - :func:`~vdjmatch.precursor.load_compendium`
     - **experimental**
     - 147 published rows in 14 source publications, tracing to 31 originating laboratories
   * - :func:`~vdjmatch.precursor.load_estimates`
     - **derived**
     - saved estimates for measured epitopes with VDJdb cognate sets

Both are mirrored at ``isalgo/airr_benchmark`` under
``vdjmatch/precursor_freq/`` and cached alongside the VDJdb releases.

.. code-block:: python

   from vdjmatch.precursor import load_compendium, load_estimates

   comp = load_compendium()                    # measured, naive rows only
   est = load_estimates(species="human", chain="TRB", model="arda")

   # both axes in cells per million of the naive compartment
   est.select("epitope_seq", (est["occupancy"] * 1e6).alias("predicted"), "measured_per_1e6")

Denominators and repeated measurements
--------------------------------------

The adult total-CD8 conversion assumes a naive fraction of one half. Frequencies
originally reported per total CD8 cells are doubled for ``freq_per_1e6_naive``.
``freq_per_1e6_published`` and ``denominator`` retain the original value and unit.

The harmonization treats cord blood as a naive compartment and leaves its
frequencies on the published scale.

Ten review rows restate primary measurements under a different denominator.
``load_compendium(dedup_retabulations=True)`` removes those rows. The default
preserves the published tabulations.

Historical generation models
----------------------------

The saved comparison scored human under **both historical** bundled sets on purpose.
Its mouse models used the historical ``arda`` family, so a
human-``olga`` against mouse-``arda`` comparison would confound species with model family.
With vdjtools 5.4.0, ``arda`` and ``learned`` alias the current ``isalgo`` fits.
Reproducing the saved comparison requires its original pinned parameter files.

The historical models differ in how they define V deletion support. OLGA stores
``P(delV | V)`` on one deletion-bin grid per locus, sized
for the longest allele, so a shorter allele carries probability on trims it cannot reach. **5.7% of
human TRA junctions score** ``Pgen`` **exactly zero** under ``olga`` — contributing nothing to any
mass and raising no error — against 19 junctions under ``arda``. The ``n_zero_pgen`` column of
:func:`~vdjmatch.precursor.summarise` exists to surface that loss.

Beyond a point estimate
-----------------------

The same two inputs — the ``Pgen`` spectrum and how many junctions were catalogued — answer several
questions that a single frequency does not.

Coverage correction uses a junction's model probability ``p`` and capture
probability ``π = 1 − exp(−Nθp)``, where ``N`` is sampling effort and
``θ`` is the fitted capture scale. The mass weight ``p/π`` approaches
``1/(Nθ)`` as ``p`` tends to zero. The corresponding richness
weight ``1/π`` diverges. :func:`~vdjmatch.precursor.unseen_junctions` marks the
richness estimate unreliable when the rarest observed junction is captured too
seldom. The sequence census below counts candidates in a specified finite
Hamming ball.

**The ball census.** ``union_mass(..., count_members=True)`` returns
``n_union`` beside the mass. ``n_union - n_seqs`` counts candidate junctions
in the union outside the supplied cognate set.

.. code-block:: python

   from vdjmatch.precursor import load_model, union_mass, observed_mass

   model = load_model("TRB", source="arda", organism="human")
   u = union_mass(model, junctions, r=1, count_members=True)

   u["n_union"] - u["n_seqs"]      # uncatalogued neighbours — exact given the ball
   observed_mass(model, junctions) / u["union"]   # share of the mass already catalogued
   u["overlap"]                    # what naive per-sequence summing would have double-counted

When cognate junctions share neighbors, summing individual ball masses counts
that shared mass repeatedly. :func:`~vdjmatch.precursor.union_mass` counts each
member of the model-defined union once. Its ``overlap`` field reached
0.41–0.61 for the best-sampled epitopes in the saved comparison.

**A sum over sampled masses depends on sample size.** If each of ``n``
independent draws selects receptor ``i`` in proportion to its mass ``p_i``,
the expected sum of the drawn masses is ``n·Σ(p_i²)/F``, where
``F = Σp_i`` is the full set mass. Use
:func:`~vdjmatch.precursor.occupancy` for expected distinct clonotype counts.

The notebook
------------

An interactive `marimo <https://marimo.io>`_ notebook walks the whole comparison — the scatter
against measurement, the ``Pgen`` spectrum inside one cognate set, the ball census, and the step from
a mass to a number of cells:

.. code-block:: bash

   pip install 'vdjmatch[precursor,notebook]'
   marimo edit docs/notebooks/precursor_compendium.py

Both comparison tables come from the mirror. The interactive ball census and
``Pgen`` diagnostics compute freshly from reference junctions with one native
thread; allow additional time for those cells and the first reference download.

Where the numbers come from
---------------------------

Every row of the compendium carries its own ``pmid`` and ``doi``: cite the laboratory that made the
measurement, not the mirror. Everything that computed a number lives in
`repseq/2026-precursor-freq <https://github.com/repseq/2026-precursor-freq>`_; the mirror's
``SOURCES.md`` gives per-file provenance, the denominator traps, and the arithmetic re-checks that
make the visually-read tables falsifiable.
