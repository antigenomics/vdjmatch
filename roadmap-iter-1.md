# vdjmatch — roadmap, iteration 1

Living plan for the current iteration. The cross-package design contract lives in
`seqtree/ROADMAP.md` §2; this file tracks what is built and what iteration 1 adds.

## Done (Phase A MVP — single-chain annotator, released as 0.0.1)

- **Repo / gitflow.** Legacy Java frozen on `legacy-java` (tags `1.1.4`–`1.3.1`); the Python rewrite
  is `dev` → `master`. GPL-3.0-or-later. CI (`tests.yml`), docs (`docs.yml` → GitHub Pages),
  PyPI publish (`publish.yml`, OIDC trusted publishing) mirror the antigenomics/tcren setup.
- **`db/`** — fetch + cache the latest VDJdb release ZIP; parse slim/full to a canonical schema.
- **`io/`** — polars readers: AIRR rearrangement (single chain) and paired (long `cell_id`/`clone_id`,
  wide TCRvdb), alias map, dedup to unique clonotypes.
- **`match/`** — `VdjdbIndex` (per-gene seqtree indices) + vectorized hit expansion; scope→SearchParams;
  bundled VDJAM matrix; CIGAR from alignment ops.
- **`evalue/` + `aggregate/`** — control-calibrated E-values (wrap `seqtree.evalues`); epitope
  enrichment + per-query best call.
- **`cli/` + `runner/` + `cluster/`** — `vdjmatch update` / `vdjmatch match`; multi-sample runner;
  pairwise sample overlap. Verified end-to-end (sample1 → NLVPMVATV top epitope; canonical TCRs
  map to GILGFVFTL/NLVPMVATV/GLCTLVAML, p≈0; junk → p=1.0). 8 unit tests pass.

Dependency note: `seqtree` 0.0.3 is not yet on PyPI; CI/publish install it from
`git+https://github.com/antigenomics/seqtree`. `pip install vdjmatch` resolves cleanly once
seqtree 0.0.3 is published.

## Iteration 1 — Task #8: paired α/β + per-region V/NDN/J + V-gene / control normalization

### Paired α/β E-values (ROADMAP §2.4)
- Joint null factorizes under chain independence: `π0^{αβ} ≈ π0^α · π0^β`; paired enrichment = Poisson
  tail on references matching **both** chains within budget (intersect per-chain `complex_id` hit sets;
  needs the fat/full VDJdb that carries pairing). Report Fisher-combined p and the joint-count p.
- Single-chain-with-partner-rarity estimate (§2.5) deferred within this task.

### Per-region V / NDN / J scoring — **and how VDJAM is derived (key decision, do not forget)**
- **VDJAM is learned for the NDN core, not for V/J.** The specificity-relevant substitution signal
  lives in the non-template (N-D-N) insert; that is where the data-driven VDJAM log-odds matrix is
  estimated (same-antigen Hamming-1 pairs, NDN-debiased — see Task #9 / `seqtree/ROADMAP.md` §2.2).
- **V- and J-region substitution behaviour is auto-computed, not learned.** For each V (and J) gene,
  take the **germline amino-acid sequence**, apply the observed **3′ V-trimming / 5′ J-trimming**
  distribution, and **fill the trimmed positions with random amino acids** (the N additions). The
  expected per-position residue distribution that results *is* the V/J segment model: it defines how
  conservative each V/J-derived CDR3 position is (germline-fixed positions are near-invariant; trimmed
  positions revert to the random-insertion background). So:
    - V/J region "matrices"/weights are a deterministic function of (germline V/J seq, trimming
      profile, random-N background) — computed, not fit.
    - Only the NDN matrix consumes VDJdb substitution data.
  This keeps the learned parameter count tiny (one NDN matrix per chain) and makes the V/J treatment
  principled and transferable across genes.
- Implementation: segment the aligned CDR3 (`seqtree.Alignment.ops`) into V-germline prefix / NDN core /
  J-germline suffix by trimming offsets; score the NDN core with the learned VDJAM and the V/J flanks
  with the computed germline+trim+random model. seqtree `PositionalMatrix` is fixed-width, so v1 does
  Python-side region rescoring of candidate hits (escalate to a native seqtree per-region scorer only
  if a benchmark proves it the bottleneck).

### V-gene incorporation + control normalization (the subtle research piece)
- Use **V + CDR3 jointly** (V+CDR3 ≈ the true motif; V carries biological bias). Compare: V as hard
  filter vs a soft S(V) term vs per-(V-cluster) null. Test on epitope data.
- **V clustering by CDR1/CDR2/CDR2.5 similarity** (so e.g. TRBV5-1 ≈ TRBV5-8): extract CDRs from
  germline and compute V-pairwise distances with seqtree (port the relevant bits of `mirpy`, whose
  `tcrtrie` is slated to become seqtree). Bundle precomputed CDR sequences + V distances.
- **Control V-J normalization tension:** a real control captures V/J bias automatically *only if its
  V-J usage matches the query repertoire*; ALICE-style P(V,J)-conditioned scaling renormalizes V-J
  usage but **risks erasing genuine V-gene biological bias** that is itself part of the specificity
  signal. Deliver **both** nulls — (i) matched real control, (ii) V-J-conditioned scaled null — and an
  experiment on held-out epitopes measuring whether renormalization helps or hurts. Decide after data.

### Verification for iteration 1
- Paired E < min single-chain E for true pairs (TCRvdb sample6).
- Region rescoring reproduces the intended V/NDN/J weighting; NDN-only-learned + computed-V/J matches
  or beats a single global matrix on VDJdb-vs-VDJdb retrieval.
- V+CDR3 fusion and the two control nulls compared on epitope PR/AUC; pick the winner.

## Later (subsequent iterations)
- Task #9: VDJAM extractor + leave-one-out-by-epitope (the NDN matrix; feeds the manuscript).
- Task #7: golden tests vs DETECT `.ods`, property tests, 300k-scale benchmark harness.
- Task #10: NAR manuscript skeleton + transfer of the seqtree E-value appendix.
