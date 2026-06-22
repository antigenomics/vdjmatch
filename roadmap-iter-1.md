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

**DONE — FR/CDR exact-match decomposition (which region carries the V prior? none).**
`bench/vregion_decompose.py` (human TRB MHC-I, germline FR1/CDR1/FR2/CDR2/FR3 from mirpy, bundled
`resources/vgene/human_v_regions.tsv`): among cross-V neighbour pairs, **no** framework or CDR region —
identical or not, and not even CDR1∧CDR2 identical (4.5%, n=22) — lifts co-specificity toward same-V
(70%); all strata sit at/below the cross-V floor (~22%), and CDR1+CDR2 edit distance shows a
flat-to-negative gradient (closer ≠ more co-specific). Decisive feature = exact gene identity (CDR3
germline anchor + integrated geometry), not transferable contacting-loop similarity. Robust at subs 1&2.

**DONE — benchmark = latest VDJdb release; gold-standard shortlist.** Standing benchmark is the latest
VDJdb release (`db.load`/`fetch_latest`; bench paths via `$VDJDB_SAMPLE`, no hardcoding). High-confidence
shortlist = clonotype–epitope pairs in **≥2 distinct references** (`db.replicated`, min_refs=2; 1214 TRB
+ 839 TRA human). Accuracy (`bench/shortlist_accuracy.py`, LOO NN vote, subs=1, exact self excluded):
**~53% top-1 on the shortlist vs ~16% on single-reference controls** (~3×); recall ~58–62% vs ~20–29%.
V-restriction helps more at wider scope (TRB subs=2: 54→58%). Replicated = public/learnable; single-ref
= private singletons. (Cf. the mhcmatch shortlist; noted in tests + memory.)

**DONE — regenerated everything on the correct release (2026-06-11-ZENODO) + composition controls.**
We had been benchmarking on a stale local export (`sample3`, 116k rows), NOT the release (`vdjdb.txt`,
284,546 records). Fixed: `db.fetch_hf` pulls the pinned release from the `isalgo/airr_benchmark` HF
dataset (uploaded gzipped); all bench scripts default to it via `_bench.source`. The dense 2026 release
is dominated by 10x mega-studies (SLLMWITQV ~30k), so `_bench.long_list` keeps epitopes ≥30 clonotypes,
caps each to a random 3000, and drops **spike studies** (one study dumping a length-skewed mega-set into
one epitope — only PMID:40498839, the SLL len-14 spike; spectratype detector). Metrics are
imbalance-robust (`bench/metrics.py`: ROC-AUC + balanced PR/F1, π0-renormalised); all tests exclude
exact self-hits. Re-run results (balanced PR-AUC unless noted):
- Hamming-1 optimum **holds** (macro purity 0.49→0.07, lift 56×→2.5×); central substitutions **sharper**
  (P(same) 0.31 core vs ~0.75 near-anchor; PSSM ~2× centre).
- No amino-acid matrix beats BLOSUM62; **BLOSUM+possig wins 7/8** (0.564→0.598 @≤2, 6/8 @≤4). tcrBLOSUM
  still ≤ BLOSUM (0.555 vs 0.567 @≤2; tied @≤4).
- **NEW — VDJAMr**: a genetic-code-only matrix (`loo_vdjam.codon_dissim`, mutational accessibility, no
  chemistry) **ties BLOSUM62** (0.581 vs 0.564) and beats data-derived VDJAM → CDR3 substitution
  structure is generative, not chemical.
- **PSSM v2** (`regions.significance_weights`): end-anchored (offset from V/J anchors, not relative),
  Beta-Binomial-smoothed; BLOSUM severity matters only in the core (P(same) 0.50 conservative vs 0.38
  radical), flat at anchors → validates the multiplicative ω·BLOSUM form.
- **V-mechanism (tempered after a deterministic re-run + tiebreaker)**: loose CDR1/2 similarity doesn't
  recover the prior (r≈0), but near-exact germline CDR1+CDR2 identity gives a real ~3× lift (edit-0 ~32%
  vs cross-V floor ~11%) — yet only ~**half** the same-V level (32% vs 64%); the rest is
  gene-identity-specific. Whole-loop, not a sparse pseudosequence (per-position flat except CDR2 pos5
  1.56×). NB the edit-0 bin is n=25 and bounced 32–60% across held-out sets; the monotone gradient over
  the well-powered bins is the robust part. (`bench/vregion_decompose.py`, `bench/vpseudo.py`)
- Shortlist accuracy ~44–47% top-1, **modestly** above single-ref controls — the earlier 3× gap was a
  sparse-export artefact (dense release makes controls findable too).
- Reproducibility: deterministic held-out-epitope selection (tiebreak by epitope name); seqtree bumped
  to **0.1.0** (`SubstitutionMatrix.penalty` API).

### OLGA spurious rate is GENUINE coincidental overlap (a same-run control gives a false 0%)
The OLGA "spurious-hit" rate estimates `1 - P(no hit)` = the probability a Pgen-drawn TCR coincides
(within first-hit scope) with the epitope reference — a coincidental-collision quantity set by Pgen,
not method noise (`bench/olga_overlap_limit.py`). On 5k **TRB** OLGA draws (sample4, chain-consistent):
raw overlap `1-P(no hit)` is scope-dependent and **large** (15.3% @1 edit → 74% @3 → 94% @5 edits;
Λ=0.52 @1 → 2546 @5); the Poisson model `P_overlap ≈ 1-exp(-Λ)`, Λ = Σ_r Pgen-mass of r's edit-ball, is
already overdispersed at k=1 (0.153 vs 0.405) → reference neighbours are **heavily clustered, not
Poisson** (convergent-recombination signal). (The old TRA-query numbers were ~100× smaller — a TRA query
rarely coincides with a TRB reference; chain matters enormously.) vdjmatch significant rate 13.2%
(TRB, M=250k).
**Control test (`bench/control_pgen_test.py`, chain-consistent — FINAL).** Critical: **`sample5` is
ALPHA chain (300k × TRAV), `sample4` is BETA (300k × TRBV)** — earlier runs used sample5 (TRA) as
queries against the TRB reference, so the "0% with a matched control" was a same-CHAIN artifact (TRA
queries vs a TRA same-run control absorbed every hit). Re-run **chain-consistent** (sample4 TRB queries,
TRB reference, all TRB controls, M=50k): **real 29.8%, same-run 30.3%, independent fresh-OLGA 28.8%** —
all three agree. The spurious rate is **genuine coincidental overlap, robust to control choice** (real /
same-run / independent); there is no control-mismatch escape hatch. (Rate scales with control size M;
compare like-for-like.) The independent control is a freshly generated TRB set
(`olga-generate_sequences --humanTRB`, different seed; `bench/out/fresh_trb_olga.tsv`), resolving the
fixed-seed caveat. **Lesson:** all TRB benchmarks must use sample4 (TRB) for OLGA negatives, not sample5
(TRA). Experiments default to 5000 cases.

**DONE — paired α/β E-value (Task #2.4) + comparison datasets (temporal holdout, TCRvdb).**
`vdjmatch.evalue.paired`: joint first-hit E-value, chain-factorized null `π0^{αβ} ≈ π0^α·π0^β`. A paired
hit needs BOTH chains within budget (radius `R = max(cost_α, cost_β)`); `E = N·(n_ctrl_α(R)/Ma)·
(n_ctrl_β(R)/Mb)`, `p = Poisson P(X ≥ n_pair(R)|E)`. `build_paired_ref` groups VDJdb by `complex_id`
(2026 release: 85,922 paired complexes). New `bench/compare.py` datasets + `--impl {vote,evalue}` arm:
- **temporal**: 2025 release (`2025-12-29`) reference, test = 2026 clonotypes new vs 2025 (held out by
  time), exact matches KEPT (cross-release match is a legit annotation). Held epitopes must have
  ≥`min_epi` records in the 2025 reference (only annotate what predates the cutoff). The honest
  real-world generalization benchmark — and it is *hard*: TRA/TRB PR-AUC 0.61/0.64, retention 16%/10%
  with the E-value arm (the subs=1 vote is near-useless here, retention ~1.7% → use the ≤5-edit
  first-hit E-value for novel TCRs).
- **tcrvdb / tcrvdb-paired**: validated TCRvdb pairs (padj<1e-5; 2 epitopes) annotated vs the 2026
  reference. **Paired beats single chain**: f1 0.89 / retention 0.82 (paired) > 0.86 / 0.77 (TRB) >
  0.75 / 0.63 (TRA), purity ~0.98 throughout — the joint α+β requirement lifts recall at high purity.

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
