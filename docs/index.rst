vdjmatch
========

Annotate selected T-cell receptors and repertoires with VDJdb peptide–MHC evidence.
The package preserves receptor identity, reports competing candidate labels and
supports explicit background calibration through a Python API and command line.
Native sequence search uses ``seqtree``; repertoire ingestion reuses ``vdjtools``.

Start with :doc:`tutorial`, follow :doc:`how-to` for release and repertoire workflows,
and read :doc:`explanation` to interpret ranking and calibration. Options and public
interfaces are documented in :doc:`cli` and :doc:`api`.

.. toctree::
   :maxdepth: 2

   tutorial
   how-to
   explanation
   releases
   cli
   api
   compendium
   notebooks/precursor

Research context
----------------

The following recorded scoring studies describe their stated historical reference
and evaluation design. They do not identify every production score as the same
estimator or make a query-level statistic a specificity posterior.

**What scoring actually buys.** An empirical study on VDJdb (the scoring appendix; 2026-06-11-ZENODO
release, composition-controlled) finds that **Hamming distance 1 is the signal:noise optimum** (macro
purity 0.49 → 0.07 over edit distance 1–5; the original VDJdb observation), that **central substitutions
carry the specificity signal** (P(same) ≈ 0.31 in the core vs ≈ 0.75 near the anchors), and that **no
amino-acid matrix clearly beats BLOSUM62** — a genetic-code null (VDJAMr) even ties it, so CDR3
substitution structure is generative, not chemical. **Reweighting BLOSUM62 by the central-position
significance** (a native ``seqtree`` ``PositionalMatrix``, end-anchored) is the one change that does beat
it, in 7/8 held-out epitopes. The substitution alphabet is second-order, position is the first-order
matrix lever, and the overall first-order statistic is the control-calibrated E-value. Finally, the
**V gene is a strong near-binary prior** (same-V neighbours share an epitope up to ~7× more than
cross-V); loose CDR1/CDR2 similarity does not recover it, and near-exact germline-loop identity recovers
only about **half** of it (the rest is gene-identity-specific).
