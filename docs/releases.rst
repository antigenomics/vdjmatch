Release notes
=============

0.4.0
-----

Standalone reference resolution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GitHub latest/pinned releases and user-supplied tables, ZIPs and directories share
one resolver. Historical legacy references, richer primary bundles and AIRR
Rearrangement/Reactivity exports preserve independent observations and metadata.
Manifest roles, checksums, atomic publication and writer locks protect fetched inputs.
An explicit primary or AIRR selection fails when the release lacks that format.

Selected-receptor evidence
~~~~~~~~~~~~~~~~~~~~~~~~~~

Query identities, loci, V/J calls and pairing survive matching and join-back.
AIRR uses anchor-inclusive ``junction_aa``; ambiguous custom inputs require an
explicit convention. Supplied counts must be nonnegative integers.

CLI/API candidate reduction reports peptide–MHC alternatives, NED-v1 ranking,
reference/study/clonotype support and competing evidence. Optional fixed-ball
E/P calibration is distinct from posterior specificity confidence. Paired
reference evidence requires both chains in the same complex and reports detailed
chain identities. Finite paired controls have bounded rule-of-three corrections.

Compatibility and limits
~~~~~~~~~~~~~~~~~~~~~~~~

CLI output now includes hits, candidates, calls, descriptive summaries, ingestion
diagnostics and a run manifest. Calls preserve input rows; exact ranking ties
abstain. Default search threads are one. Samples run sequentially and native
search uses complete query batches per locus. VDJAM is the CLI/API default scoring matrix.

Requires ``seqtree>=1.0.2``, ``vdjtools>=5.0.0``, ``polars>=1.41.2,<2`` and
``scipy>=1.9`` for exact finite-control bounds and binomial tails.
HF benchmark snapshots are explicit, separate from GitHub release resolution.
Hard V/J calibration requires metadata-aware controls and fails explicitly.
Paired CLI V/J restrictions are currently unsupported and fail explicitly.
``first-hit`` exposes global adaptive first-hit enrichment through the CLI,
with fresh raw control indexing, original query rows and a provenance manifest.
``match --fresh-control`` and repeatable ``--control LOCUS=TABLE`` support fresh
raw-control manuscript runs. Single-chain ``--exclude-exact`` applies the same
puncture to target and control evidence. Paired exclusion removes only exact/exact
junction pairs and retains one-exact/other-neighbour evidence. Its all-count
finite-control method uses four marginal category Clopper–Pearson upper bounds,
fixed total failure budget ``1e-6``, and a binomial target tail plus that budget.
Empirical Cartesian counts, raw E and raw Poisson tails remain diagnostics.
The bound is explicitly conditional on IID category/target sampling and
independent background chains; deduplicated controls do not establish those
assumptions or unconditional biological calibration. Default paired Poisson
calibration remains unchanged.
Declared paired reference links survive string IDs, independent chain record IDs
and repeated normalization; conflicting links fail explicitly.
``search`` exposes fresh native batch matching, positional row identities,
counts-only output, same-key filtering and vectorized positional matrix scores.
Its manifest separates native build/search and positional rescoring times.
Repertoire-level statistical inference remains future work.

Graded neighbourhood annotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``match --search-mode ball`` exposes fixed radii 1–5 without a global neighbour
cap. Unit-cost balls share one minimum-edit search; weighted balls use a complete
batch per radius. Fixed paired annotation retains its explicit joint control
calibration. Existing fixed-mode behaviour remains available.

The experimental fixed-K CLI/sample-runner mode and ``--top-k`` argument were
removed before publication. Historical ranked research APIs remain available for
reproducing earlier diagnostics. Their calibration and label ranks must not be
presented as the production neighbourhood scorer. Publishing seqtree 1.0.2
remains a dependency prerequisite for the pending release.

Repeated query junctions are searched once within each native batch, including
single-chain and paired marginal controls. Hits are expanded back to every query
row before evidence reduction, preserving distinct V/J calls, counts and pairing.
No computed search results are persisted between runs.

V-gene matching and germline model lookup preserve legitimate IMGT slash-containing
names (for example TRAV23/DV6 and TRBV24/OR9-2) while stripping allele suffixes.
These names must not be treated as decorations or collapsed to a different gene.
Historical benchmark exports that stripped slash suffixes retain their original
provenance and require explicit reconciliation before a fresh scientific comparison.
