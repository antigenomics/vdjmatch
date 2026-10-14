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

Is this junction present or unusually close to the database?
-------------------------------------------------------------

Literal presence is a lookup. ``database-evidence`` separately reports the
significance of similarity to *any* valid database junction, pooling specificity
labels within the declared species and chain:

.. code-block:: sh

   vdjmatch database-evidence query.airr.tsv --vdjdb reference.zip \
     --species human --locus TRB --null-sample null-junctions.airr.tsv \
     --threads 4 --output-prefix presence

Read ``presence.database.tsv``:

.. list-table:: Database-level diagnostics
   :header-rows: 1

   * - Column
     - Question answered
   * - ``database_exact_junction``
     - Does this complete anchor-inclusive junction occur in the database? A boolean lookup, ignoring V/J and antigen labels.
   * - ``database_exact_gene_clonotype``
     - Does the same junction/V-gene/J-gene combination occur? Alleles are pooled; missing or ambiguous gene calls remain unavailable.
   * - ``database_nearest_edits``
     - How many unit edits separate it from the closest database junction within the declared search limits?
   * - ``database_p``
     - Would a matching null receptor be at least as close to any database junction? Smaller values indicate more unusual resemblance.
   * - ``database_status``
     - ``matched`` means a neighbour exists within the edit limit; it is not a calibrated call. Other values report no-hit, missing/invalid sequence or species/chain mismatch.

For a fixed reference, let :math:`d(x)` be a query's nearest admissible junction
edit distance. For :math:`N` supplied null draws, :math:`k` counts draws with
:math:`d(z) \leq d(x)`. Then :math:`p_{\mathrm{database}}=(k+1)/(N+1)`.
Exact ties count against the query; repeated draws and null no-hits remain in
the denominator. A query with no hit receives P = 1. Invalid query rows retain
unavailable diagnostics. Every null row must be valid for the declared species
and chain; the command rejects a null table that would silently lose rows.
With multiple ``--null-sample`` files, it reports the largest source-specific
tail. ``presence.null.tsv`` and ``presence.null-tails.tsv`` retain the null
distances and each source's numerator/denominator.
For a query no-hit, no finite radius is tested: its tail numerator is unavailable
and P = 1 is the conservative no-evidence convention.

This calculation retains exact records, applies no study-support cutoff, and
does not select an epitope. It uses the same complete database search for queries
and nulls, so selection over database junctions is included in its null statistic.
Its P-value assumes query/null exchangeability under the declared null and a
fixed reference. A common exact junction can have a large P-value; an unusual
neighbour can have a small one. Neither identifies the correct antigen.
``--max-edits`` sets the full-junction Levenshtein search limit (default 5);
substitutions, insertions and deletions each cost one. The diagnostic is separate
from ``kernel_score`` and the scoring pipeline's
``p_null``. TRA and TRB diagnostics remain separate; they cannot supply a
paired-receptor P-value. Its P-value ignores V/J; the gene-clonotype lookup is
a separate exact-presence flag.

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
