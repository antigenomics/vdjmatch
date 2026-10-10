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

VDJdb reference observations
----------------------------------------

Use ``--vdjdb`` to read a frozen local legacy or rich reference ZIP through the
same database loader used for annotation::

   vdjmatch pgen reference.zip --vdjdb --model-path /path/to/raw/model \
       --species human --locus TRA --threads 12 --output-prefix reference-tra

The loader selects the declared receptor species and locus. Every normalized
reference observation in that selection is retained, including repeated junctions
and its ``reference_id``, peptide, MHC and other metadata. The command projects
internal ``cdr3/v/j/gene`` fields to ``junction_aa/v_call/j_call/locus`` without
renaming or dropping the original columns. Legacy VDJdb ``cdr3`` is a full
anchor-inclusive junction under the database loader's source contract.
Nested rich metadata is exported as JSON cells using the existing flat-table
exporter. ``query_id`` is the selected loader row ordinal when not supplied.

The manifest declares ``input_format=vdjdb`` and records ``db.provenance`` together
with the reference loader/normalizer hashes. Row counts refer to the selected
normalized reference observations, rather than all records in the archive.
Pgen remains marginal over the model's V/J alleles; database V/J annotations,
epitopes and study metadata do not condition it. Local reference tables and
release directories are also supported by the existing loader. AIRR remains the
default when ``--vdjdb`` is absent.
