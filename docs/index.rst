vdjmatch
========

Fast, control-calibrated annotation of **T-cell receptor antigen specificity**.

``vdjmatch`` annotates clonotypes in large AIRR repertoires against
`VDJdb <https://github.com/antigenomics/vdjdb-db>`_ by fuzzy CDR3 search, reporting a
**control-calibrated E-value** (BLAST-style significance against a background repertoire) and
enriched antigen-specificity labels. It is a Python rewrite of the legacy Java/Groovy vdjmatch,
built on the `seqtree <https://github.com/antigenomics/seqtree>`_ search core.

.. note::

   ``vdjmatch`` 0.0.1 is an early release. The single-chain annotator (VDJdb fetch, AIRR I/O,
   fuzzy search, E-values, epitope-enrichment summaries, CLI) is in place; paired α/β scoring, a
   re-derived segment-aware substitution matrix (VDJAM), and tool comparisons are in progress.

Installation
------------

.. code-block:: bash

   pip install vdjmatch

``seqtree`` (the search engine) is installed as a dependency. For development:

.. code-block:: bash

   python -m venv .venv && source .venv/bin/activate
   pip install -e ".[test,bench]"

Quickstart
----------

Fetch the latest VDJdb release and annotate an AIRR rearrangement sample:

.. code-block:: bash

   vdjmatch update                                   # cache the latest VDJdb release
   vdjmatch match --species HomoSapiens --scope 1,0,0,1 -o out sample_airr.tsv

This writes three tab-separated tables per sample:

- ``out.<sample>.hits.txt`` — every query→VDJdb hit with CDR3 alignment, CIGAR, edit counts and score.
- ``out.<sample>.calls.txt`` — one predicted epitope per query clonotype with its E-value.
- ``out.<sample>.summary.txt`` — epitope-level enrichment (unique clonotypes, reads) by MHC class
  and antigen species.

Key ideas
---------

**Control-calibrated E-value.** Immune repertoires are biologically redundant (convergent
recombination, public clones), so a naive i.i.d. null massively over-calls. ``vdjmatch`` counts a
query's VDJdb neighbours within a fixed search scope and compares to the count expected from a
matched **background control** repertoire; the Poisson-tail ``p_enrichment`` is significant only
when a clonotype has *more* VDJdb neighbours than the generative process predicts — the hallmark of
antigen-driven selection. The theory is derived in the ``seqtree`` appendix.

**Scope / budget.** ``--scope s,i,d,t`` sets the maximum substitutions, insertions, deletions and
total edits of the CDR3 search ball.

**VDJAM.** A TCR-specific amino-acid substitution matrix (bundled), with optional region-aware
weighting that emphasises the antigen-contacting NDN core over the germline-fixed V/J flanks
(germline-retention profiles derived from the OLGA model via ``mirpy``).

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

.. toctree::
   :hidden:

   self
