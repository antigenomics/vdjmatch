# CLI-first synthetic examples

Run from a fresh installation with `vdjmatch` on PATH:

```sh
sh examples/annotate.sh NEW_ANNOTATION_DIRECTORY
sh examples/compare.sh NEW_COMPARISON_DIRECTORY
```

Keep this complete directory together. The scripts find their raw inputs beside
themselves, use the offline reference, default to one native thread and fail if
the output directory exists. Every invocation calculates fresh CLI tables;
there are no feature caches, predictions as inputs, skip or resume paths.

The annotation notebook covers a four-row metadata cohort and two wide pairs,
including an incomplete pair. The comparison notebook calls the CLI with exact
exclusion and detailed default stem4 / advanced vcdr12 settings, then joins its
continuous evidence on source identity and complete pMHC keys. Open either
notebook from the repository root or this directory. No outputs are committed.

All inputs were written specifically as deterministic software fixtures, not
sampled from a biological dataset. `DEMO_A`/`DEMO_B` are demonstration labels,
not experimental binding truth. Literal junctions include conserved anchors.
There is no RNG or fit and no accuracy or calibration claim.

Offline reference SHA256: `3ed6c1baafe2ec51bac658a02fa806d4a866ffcffc1ff7e1a391188050bcb76d`.
The Git version of this directory pins the raw inputs; CLI manifests record the
reference and scoring configuration of each fresh invocation.
