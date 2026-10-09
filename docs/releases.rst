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
search uses one batch per locus. VDJAM is the CLI/API default scoring matrix.

Requires ``seqtree>=1.0.1``, ``vdjtools>=5.0.0`` and ``polars>=1.41.2,<2``.
HF benchmark snapshots are explicit, separate from GitHub release resolution.
Hard V/J calibration requires metadata-aware controls and fails explicitly.
Paired CLI V/J restrictions are currently unsupported and fail explicitly.
``first-hit`` exposes global adaptive first-hit enrichment through the CLI,
with fresh raw control indexing, original query rows and a provenance manifest.
Repertoire-level statistical inference remains future work.
