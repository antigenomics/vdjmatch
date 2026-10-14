Read scores and decide whether to call a match
==============================================

Use ``kernel_score`` to rank evidence and a matching ``p_null`` to judge how
unusual the selected score is under the declared null. ``call_status`` applies
the frozen operating cutoff. These quantities answer different questions.

Which output should I read?
---------------------------

.. list-table:: Current annotation output
   :header-rows: 1

   * - Column
     - Read it as
   * - ``selected_epitope``
     - Highest-scoring peptide for a single-chain query; inspect ties before assigning it.
   * - ``kernel_score``
     - Relative reference evidence. Larger values rank higher.
   * - ``kernel_competing_score``
     - Strongest alternative peptide's score, on the same scale.
   * - ``tied_peptides``
     - Number of peptide labels sharing the best score. More than one means ambiguity.
   * - ``selection_status``
     - Whether ranking evidence is available; ``ranked_evidence`` is not a calibrated call.
   * - ``kernel_support_share``
     - Best score divided by best plus strongest competing score. Not binding probability.
   * - ``p_null``
     - Empirical upper-tail probability of the selected score under matching null queries.
   * - ``call_status``
     - ``accepted`` or ``no_supported_match`` after applying the matching calibration.

The :doc:`tutorial` shows actual single-chain and paired input/output rows.
The :doc:`algorithms` page explains how each route produces its evidence.

Obtain a calibrated decision
----------------------------

Supply the frozen calibration and its original null-case table for the same
reference and complete scoring pipeline. The input selection table must have
one unique ``query_id`` per query and the score column recorded in the recipe.
For a recipe fitted with ``--score-column kernel_score``:

.. code-block:: sh

   vdjmatch apply-call-calibration annotation.selection.tsv \
     --calibration calibration.calibration.json \
     --input-manifest annotation.manifest.json \
     --null-cases calibration-cases.tsv \
     --output calls.tsv

``calls.tsv`` preserves the ranking columns and adds ``p_null``, ``E_null``,
``call_score`` and ``call_status``. The adjacent ``.null-tails.tsv`` gives the
numerator and denominator separately for every null source. The command checks
producer identity for bound recipes and the exact null-case file hash.

A recipe fitted to a column named ``score`` requires that column; a recipe fitted
to ``kernel_score`` requires ``kernel_score``. Do not rename columns to bypass a
mismatch in the underlying scorer, reference or selection procedure.

What does the P-value test?
---------------------------

For one null source, let :math:`N` be all null queries, including those with no
support, and :math:`k` the number with selected scores at least as high as the
query. Ties count against the query:

.. math::

   p_{\mathrm{null}} = \frac{k+1}{N+1}.

This is a conservative rank P-value when query and null scores are exchangeable
under the same complete procedure. With several null sources, the output takes
the largest tail probability. ``E_null = k/N`` is the empirical false-call rate
per query at that score, not an upper confidence bound. Unsupported queries
receive ``p_null = 1`` and no accepted call.

A small P-value means the selected evidence is unusual under the declared null.
It does not establish experimental binding or resolve an exact tie between
peptides. The configured cutoff controls the measured any-call rate in the
supplied calibration populations; an arbitrary score cutoff has no such meaning.

Single-chain and paired decisions
---------------------------------

For a single chain, the complete selection maximizes over restrictions within
each peptide, then over peptides. Its null must include that entire search.
Testing the chosen peptide as though it were specified beforehand changes the
question and understates selection.

For a pair, ``selected_linked_epitope``, ``selected_alpha_epitope`` and
``selected_beta_epitope`` remain separate. Linked evidence uses the same
reference pair. Marginal evidence includes independent same-target single-chain
references. A missing chain does not erase the available chain's evidence.

``kernel_any_axis_score`` is the maximum of those three selected scores. Its
calibration therefore includes selection across all three axes. A single-chain
or linked-only P-value cannot be substituted for that paired any-axis P-value.

Choose the matching recipe
--------------------------

A calibration belongs to its species, chain population, reference identity,
reference filtering, geometry, scorer, exact-match policy and label-selection
procedure. Ordinary annotation retains exact matches. A benchmark using
``--exclude-exact`` needs its own calibration. Restricting the candidate peptide
or HLA changes the selection question too.

If a matching calibration is unavailable, the output still supports ranking,
but it cannot supply a validated match decision. Do not turn ``kernel_score``
into a P-value by rescaling it, or substitute a component probability.

Retained methods
----------------

``historical-density`` and ``match`` expose earlier scoring models. Their
reference-count, density and probability fields have different definitions.
See :doc:`legacy-explanation` for those interfaces and :doc:`scoring` for the
current kernel and retained alternative formulas.

.. toctree::
   :hidden:

   legacy-explanation
