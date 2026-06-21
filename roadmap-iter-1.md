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

### Result so far (bench/gen_vdjam.py — position-debiased CDR3 substitution log-odds)
Per-region correlation of the learned matrix with BLOSUM62, on VDJdb same-antigen Hamming-1 pairs
(human; min 30 CDR3s/epitope), **empirically confirms the NDN-vs-V/J split**:

| chain | region | events | r(BLOSUM62) |
|-------|--------|--------|-------------|
| TRB | NDN | 6584 | **0.42** |
| TRB | V   |  749 | 0.08 |
| TRB | J   | 2797 | 0.22 |
| TRA | NDN | 5531 | **0.35** |
| TRA | V   | 8083 | 0.07 |

NDN carries the BLOSUM-like free-substitution signal; the V region is germline-dominated (≈no free
substitution chemistry) → confirms VDJAM should be **learned on NDN** and V/J **auto-computed from
germline + trimming**.

**DONE — germline+trimming V/J model + region-aware scoring.** Added `mir.basic.trimming`
(OLGA-derived per-gene germline-retention profiles + germline-match-aware `PgenLite`, mirpy v1.3.0);
emitted the retention profiles to `src/vdjmatch/resources/trimming/human_vj_retention.tsv`. New
`vdjmatch.match.regions`: per-position weight `1 − P(germline-retained)` (NDN core → ~1, V/J flanks →
~0) + VDJAM penalties + region-aware substitution/aligned scoring; wired into `engine.annotate`
(`region_aware=True` → `region_score` column). Leave-one-out-by-epitope (TRB, `bench/loo_vdjam.py`,
4 arms): **region weighting beats flat VDJAM (+0.004 @scope2, +0.007 @scope3, never worse, grows with
scope) but BLOSUM62 stays strongest** at this narrow-scope retrieval — CDR3 variation is already
NDN-concentrated, so flank down-weighting moves few hits. Documented honestly in `appendix/scoring.tex`
(retention figure + 4-arm LOO figure). **Next lever:** the V-gene dimension (V clustering by
CDR1/2/2.5; cross-V comparisons, where flank weighting should matter) + control V–J normalization.

### Verification for iteration 1
- Paired E < min single-chain E for true pairs (TCRvdb sample6).
- Region rescoring reproduces the intended V/NDN/J weighting; NDN-only-learned + computed-V/J matches
  or beats a single global matrix on VDJdb-vs-VDJdb retrieval.
- V+CDR3 fusion and the two control nulls compared on epitope PR/AUC; pick the winner.

**DONE — central-significance PSSM beats BLOSUM62.** Encoded experiment-2 (central CDR3 mismatches
carry the specificity signal) as a native seqtree `PositionalMatrix.from_weights(BLOSUM62, ω)` (centre
~1.4×, V/J borders ~0.7×; `regions.significance_pssm`). Beats flat BLOSUM62 in **8/8** held-out epitopes
at edit distance ≤2 (0.333→0.356) and ≤4 (0.218→0.236). The substitution *alphabet* is second-order;
*position* is the first-order matrix lever.

**DONE — refuted tcrBLOSUM (Postovskaya et al., Brief Bioinform 2024, doi:10.1093/bib/bbae602).** Their
claim that a TCR-specific matrix beats BLOSUM62 does not survive held-out-epitope retrieval
(`bench/tcrblosum_refute.py`): tcrBLOSUMb scores **below** BLOSUM62 (0.303 vs 0.330 @≤2, 0/8 wins;
0.197 vs 0.211 @≤4), while their bundled BLOSUM62 reproduces ours exactly. Their counts use no
redundancy clustering and are validated on the same same-epitope relation they are fit to → in-sample.

**DONE — V-gene dimension is a strong but near-binary prior; CDR1/2 clustering refuted.** Comprehensive
scan over {human,mouse}×{TRA,TRB}×{MHC-I,MHC-II} (`bench/vgene_scan.py`, `vgene_strat.py`;
`match/vgene.py` ships germline-loop similarity + `v_clusters`): a V match is a **1.4–17× co-specificity
prior** in every cell, but germline CDR1+CDR2 similarity barely predicts cross-V co-specificity
(point-biserial r ∈ [−0.06, +0.13] in every well-powered cell). Soft V-clustering by CDR1/2 does **not**
recover the prior — V is gene identity, not loop chemistry. Open work: model V as a strong prior in a
joint V+CDR3 E-value (not soften it); per-chain TRA/TRB; control V–J normalization.

## Later (subsequent iterations)
- **Hard vs easy (featured/featureless) epitopes.** Per-epitope PR-AUC varies enormously (convergent
  CMV NLV ≫ diffuse influenza GIL); this is biology, not noise. Build an *a-priori* "annotability"
  score from repertoire convergence (and, where available, pMHC structural features) to report
  calibrated per-epitope confidence and triage hard cases. Featured pMHC select focused motif-rich
  repertoires; featureless pMHC give structurally diverse TCRs and no tight motif — clustering-based
  inference succeeds/fails accordingly. Refs (verified via PubMed):
  - Turner et al. 2005, *Nat Immunol* 6(4):382–389, doi:10.1038/ni1175 (PMID 15735650) — featureless
    pMHC → limited/public repertoire; mutating one feature collapses diversity.
  - Yang et al. 2017, *J Biol Chem* 292(45):18618–18627, doi:10.1074/jbc.M117.810382 (PMID 28931605) —
    GIL/HLA-A2 recognised by structurally distinct TCRs in different docking modes.
  - Song et al. 2017, *Nat Struct Mol Biol* 24(4):395–406, doi:10.1038/nsmb.3383 (PMC5383516) — broad
    repertoire / diverse structural solutions for the featureless M1/HLA-A2 epitope.
  - Peng et al. 2014, *PNAS* 111(26):E2656–E2665, doi:10.1073/pnas.1401131111 (PMC4084487) — biophysics
    of recognising "featureless" surfaces (antibody; transferable principle).
  - Hudson et al. 2024, *Immunoinformatics* 13:100033, doi:10.1016/j.immuno.2024.100033 (PMC10955519) —
    clustering models for TCR specificity succeed/fail by epitope.
- Task #9: VDJAM extractor + leave-one-out-by-epitope (the NDN matrix; feeds the manuscript).
- Task #7: golden tests vs DETECT `.ods`, property tests, 300k-scale benchmark harness.
- Task #10: NAR manuscript skeleton + transfer of the seqtree E-value appendix.
