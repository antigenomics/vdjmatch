:orphan:

Marginal junction generation probability
========================================

The ``pgen`` command computes the amino-acid generation probability of each
full, anchor-inclusive AIRR ``junction_aa`` under an explicitly selected raw
``vdjtools`` model directory::

   vdjmatch pgen input.airr.tsv --model-path /path/to/raw/model \
       --species human --locus TRA --threads 12 --output-prefix result

Use a model whose manifest matches the declared species and locus. The command
loads and validates that directory directly; it neither selects a bundled model
implicitly nor collapses its alleles. ``vdjtools.load_model`` applies its documented
anchor repairs. The manifest records the consumed model-file hashes, model metadata,
loader and native binary hashes, options, input hash, row dispositions and timing.
No model download or new model fitting occurs.

This is marginal ``P(junction_aa)``: all V/J alleles and recombination scenarios in
the model contribute, with zero allowed amino-acid mismatches. Input V/J annotations
are preserved but do not condition the probability. This is neither an antigen
specificity probability nor a probability calibrated from the annotation score.
The model's generative mass does not include post-recombination selection merely
because the input repertoire was observed or productive. Paired input remains one
result per chain row; the command does not infer a joint pair probability.

``result.scores.tsv`` preserves all raw rows, metadata and supplied unique,
non-null ``query_id`` values. If absent, ``query_id`` is the zero-based raw row
ordinal. ``junction_aa_normalized`` records uppercase, whitespace-stripped input;
anchors are never fabricated from bare ``cdr3_aa``. Canonical amino-acid junctions
of the selected species/locus are sent to one native batch. Missing locus metadata
uses the explicit ``--locus``; missing V/J calls do not prevent marginal Pgen.
Explicitly different species/loci, missing junctions and noncanonical sequences
remain in the table with an explanatory ``status`` and null probabilities.

``pgen`` is the marginal probability and ``log10_pgen`` is its base-10 logarithm
when positive. ``available`` is true for a valid numeric result, including zero.
``zero_mass`` records a returned zero, with null logarithm; it does not distinguish
impossible model mass from numerical underflow. Negative, nonfinite or greater-than-one
native probabilities fail the command instead of being clipped or exported.
``result.manifest.json`` records the exact computation contract. Positive
``--threads`` caps native workers and the CLI's table-thread budget; no outer worker
pool or computed-intermediate cache is used.
