Validation data: measured precursor frequencies
===============================================

``F(e)`` estimates a **naive repertoire mass** — the fraction of a naive repertoire that can see an
epitope. Other laboratories have measured that quantity directly, by tetramer and dextramer
enrichment of unexposed donors, umbilical cord blood and naive mice, so the estimator can be checked
against the estimand itself rather than against a proxy.

Two tables ship for that comparison, and they are deliberately kept apart:

.. list-table::
   :header-rows: 1
   :widths: 22 16 62

   * - loader
     - provenance
     - what it is
   * - :func:`~vdjmatch.precursor.load_compendium`
     - **experimental**
     - 147 published measurements, 14 studies, 31 originating laboratories
   * - :func:`~vdjmatch.precursor.load_estimates`
     - **derived**
     - what ``vdjmatch.precursor`` predicts for the epitopes VDJdb can supply a cognate set for

Merging them into one column would conflate a measurement with a prediction of it, which is the
whole point of the comparison. Both are mirrored at ``isalgo/airr_benchmark`` under
``vdjmatch/precursor_freq/`` and cached alongside the VDJdb releases.

.. code-block:: python

   from vdjmatch.precursor import load_compendium, load_estimates

   comp = load_compendium()                    # measured, naive rows only
   est = load_estimates(species="human", chain="TRB", model="arda")

   # both axes in cells per million of the naive compartment
   est.select("epitope_seq", (est["occupancy"] * 1e6).alias("predicted"), "measured_per_1e6")

Three things that will bite you
-------------------------------

**Total CD8 is not naive CD8.** About half of adult circulating CD8 T cells are naive-phenotype, so
a frequency "per CD8" and one "per naive CD8" differ by 2×. ``freq_per_1e6_naive`` is the harmonized
column and the conversion is already applied; ``freq_per_1e6_published`` and ``denominator`` keep it
reversible. Do not apply it twice.

**Cord blood is not adult blood.** Umbilical cord blood is essentially the naive compartment in
full, which is exactly why it is used to measure precursors. The adult 0.5 factor must not be
applied there, and is not.

**A review row may restate a primary row.** Ten do — the same measurement under a different
denominator. ``load_compendium(dedup_retabulations=True)`` drops them; it is off by default because
a per-row median legitimately weights by how many donors stand behind a row. Turn it on when you are
counting *measurements* rather than weighting them.

The model set is not interchangeable
------------------------------------

Human is scored under **both** bundled sets on purpose. Mouse exists only under ``arda``, so a
human-``olga`` against mouse-``arda`` comparison would confound species with model family.

They differ where it matters: OLGA stores ``P(delV | V)`` on one deletion-bin grid per locus, sized
for the longest allele, so a shorter allele carries probability on trims it cannot reach. **5.7% of
human TRA junctions score** ``Pgen`` **exactly zero** under ``olga`` — contributing nothing to any
mass and raising no error — against 19 junctions under ``arda``. The ``n_zero_pgen`` column of
:func:`~vdjmatch.precursor.summarise` exists to surface that loss.

Beyond a point estimate
-----------------------

The same two inputs — the ``Pgen`` spectrum and how many junctions were catalogued — answer several
questions that a single frequency does not.

**Unseen mass and unseen count are different questions.** The mass converges: the Horvitz–Thompson
weight ``p/π → 1/N`` as ``p → 0``, which is what lets
:func:`~vdjmatch.precursor.coverage_corrected_mass` work without guessing the shape of the tail. The
count does not: ``1/π`` diverges, and a richness extrapolation on real TCR data returns 10¹¹–10¹⁴
unseen junctions, which is not a number. So two answers ship —
:func:`~vdjmatch.precursor.unseen_junctions` marks ``richness_reliable`` false when the rarest
observed junction was captured too seldom, and the *ball census* below is the finite, exact
alternative.

**The ball census.** ``union_mass(..., count_members=True)`` returns ``n_union`` beside the mass;
``n_union - n_seqs`` counts candidate cognate junctions no database holds.

.. code-block:: python

   from vdjmatch.precursor import load_model, union_mass, observed_mass

   model = load_model("TRB", source="arda", organism="human")
   u = union_mass(model, junctions, r=1, count_members=True)

   u["n_union"] - u["n_seqs"]      # uncatalogued neighbours — exact given the ball
   observed_mass(model, junctions) / u["union"]   # share of the mass already catalogued
   u["overlap"]                    # what naive per-sequence summing would have double-counted

**The overlap is a measurement, not an inconvenience.** Cognate TCRs are near-duplicates by
construction, so their balls overlap and the *sum* of per-sequence ball masses double-counts.
:func:`~vdjmatch.precursor.union_mass` is exact without inclusion–exclusion, hence carries no
truncation error, and ``overlap`` quantifies how tight the specificity group is — it reaches
0.41–0.61 for the best-sampled epitopes.

**A set total is the wrong functional.** Under the size bias that governs which receptors get
catalogued, a summed mass has expectation ``n·(Σπ²)/F`` — proportional to how many receptors were
catalogued and *inversely* proportional to the quantity it is meant to estimate. Count distinct
clonotypes through the saturating :func:`~vdjmatch.precursor.occupancy` model for anything
count-shaped.

The notebook
------------

An interactive `marimo <https://marimo.io>`_ notebook walks the whole comparison — the scatter
against measurement, the ``Pgen`` spectrum inside one cognate set, the ball census, and the step from
a mass to a number of cells:

.. code-block:: bash

   pip install 'vdjmatch[precursor,notebook]'
   marimo edit docs/notebooks/precursor_compendium.py

It recomputes no ``Pgen`` for the comparison itself — both tables come from the mirror — so it opens
in a second.

Where the numbers come from
---------------------------

Every row of the compendium carries its own ``pmid`` and ``doi``: cite the laboratory that made the
measurement, not the mirror. Everything that computed a number lives in
`repseq/2026-precursor-freq <https://github.com/repseq/2026-precursor-freq>`_; the mirror's
``SOURCES.md`` gives per-file provenance, the denominator traps, and the arithmetic re-checks that
make the visually-read tables falsifiable.
