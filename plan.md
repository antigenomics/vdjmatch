# vdjmatch 2.0 (Python on seqtree) + Nucleic Acids Research manuscript

## Context

The legacy `vdjmatch` (Java/Groovy, `~/vcs/code/vdjmatch`, latest tag `1.3.1`) annotates TCR clonotypes
by antigen specificity against VDJdb using a cloglog GLM over V/CDR1/CDR2/J/CDR3 similarities, a custom
TCR substitution matrix (VDJAM), and an informativeness weight for "public" TCRs. It was never properly
published and is slow/heavy. We are **rewriting it in Python on top of `seqtree`** (0.0.3b1, on PyPI,
GPL-3.0-or-later) — which already provides fast fuzzy search (`Index`/`KmerIndex`), control-calibrated
E-values (`evalues`/`load_control`, full theory in `seqtree/appendix/evalue.tex`), substitution &
positional matrices, alignment with M/S/I/D ops, and batched multithreaded search — and writing an NAR
paper. Goal: extremely fast, light, highly tested annotation of huge AIRR samples; E-values + enriched
epitope labels; single- and paired-chain (α/β) scoring; later a vdjdb-web backend.

**Decisions locked (this session):** (1) License = **GPL-3.0-or-later** (required: imports GPL-3 seqtree).
(2) Headline statistic = **control-calibrated E-value primary**, with the legacy GLM kept as an optional
similarity/ranking score, **and a dedicated V-gene workstream** (V+CDR3 as a joint motif; V-clustering
by CDR1/CDR2/CDR2.5 similarity; the control V-J-usage normalization problem). (3) First milestone = **all
of**: git migration + software (single **and** paired chain) + manuscript skeleton + VDJAM extractor +
porting/expanding the seqtree appendix theory into the manuscript.

Owner hardware: Apple M3, 32 GB, repo-local `.venv`. Gitflow: feature → `dev` → `master`; commit trailer
`Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`; no PyPI publish without explicit release.

---

## A. Git migration (repo `~/vcs/code/vdjmatch`) — execute on approval

Latest Java tag `1.3.1`; branches `master` + experimental. Steps:
1. `git switch -c legacy-java master`; push. Old tags `1.1.4…1.3.1` keep history + releases accessible.
2. **Fix the old release on `legacy-java`**: gradle wrapper is `gradle-5.0-rc-1` (broken). Bump the wrapper
   to a working Gradle, pin a JDK it builds under, get `./gradlew build` green, tag a final `1.3.2`
   (or `java-final`). This satisfies "fix old release."
3. `git switch -c dev master`; do the rewrite here. First rewrite commit `git rm -r src build.gradle
   gradlew* gradle jitpack.yaml .travis.yml benchmark` (legacy deleted from `dev`, recoverable from
   `legacy-java`/tags). Replace `LICENSE`/`README.md` for the Python tool (GPL-3-or-later).
4. When the first milestone is stable: `git switch master && git merge --no-ff dev`; tag `2.0.0a1`.
   PyPI release deferred to explicit request.

---

## B. Software: package layout (on `dev`)

```
vdjmatch/
  pyproject.toml          # hatchling; deps: seqtree>=0.0.3, polars>=1.0; extras: test/bench/control/ods
  src/vdjmatch/
    db/      vdjdb.py cache.py schema.py        # fetch latest GitHub release, cache by tag+ETag, polars parse
    io/      airr.py columns.py writer.py        # polars readers: rearrangement + cell(paired); alias map; writers
    match/   engine.py scope.py scoring.py regions.py vgene.py cigar.py hits.py
    evalue/  single.py paired.py partner_prior.py control.py
    aggregate/ enrichment.py glm.py
    cluster/ overlap.py
    runner/  multisample.py
    cli/     __main__.py (match, cluster, update, annotate)
    resources/ vdjam.txt score_coef.txt segm_score.txt cdr_seqs/ VERSIONS.md
  tests/ unit/ golden/ property/ conftest.py
  bench/ bench_annotate.py gen_vdjam.py harness.py
```
Pure Python by default. C++ (scikit-build-core + pybind11, or a seqtree bump) only if a benchmark proves
a hot path is the bottleneck (§G).

### Modules (reuse seqtree; do not reimplement search/E-values/matrices/alignment)
- **db/vdjdb.py** — `fetch_latest(asset="slim"|"fat", pin=tag)` via stdlib `urllib` against
  `api.github.com/repos/antigenomics/vdjdb-db/releases/latest`; cache keyed by tag (reproducible `pin`).
  polars parse to canonical cols (`cdr3, v.segm, j.segm, antigen.epitope, antigen.species, antigen.gene,
  mhc.a/mhc.b, mhc.class, vdjdb.score, complex.id`). `complex.id` carries α/β pairing (fat/full only) —
  required for paired mode.
- **io/airr.py + columns.py** — polars `scan_csv`/`read_csv`; alias map handles AIRR-canonical
  (`junction_aa`,`v_call`,`j_call`,`locus`,`duplicate_count`,`cell_id`) **and** test-data short forms
  (`v_gene`,`j_gene`,`cdr3_*_aa`,`TRBV`,`clone_id`). `read_cell()` pivots α/β by `cell_id`; ingest
  `sample6_TCRvdb.csv`. **Dedup to unique clonotypes first** (E-values require it), keep multiplicity
  for read-weighted summaries.
- **match/engine.py** — per (species,gene) partition: `Index.build(unique_cdr3s,"aa")`; `ref_id`→DB row
  via positional polars columns (no per-hit dicts). `search_batch(queries,params,threads=0)` (GIL-released).
  For 100k–1M query sets: `KmerIndex.seed_and_gather` candidate gen → `search_batch`/`align` only
  candidates. **Grouped index per (V,J)** (mirpy `tcrtrie` pattern, ~1.5–2×) = planned scaling lever.
- **match/scope.py** — legacy `s,i,d,t` → `SearchParams(max_subs,max_ins,max_dels,max_total_edits, engine,
  matrix, pos_matrix, max_penalty)`. Default `seqtm`; `seqtrie`+matrix for similarity-scored mode.
- **match/cigar.py** — RLE `Alignment.ops` (M/S/I/D) → CIGAR (`=`/`X`/`I`/`D`) + aligned match/gap strings
  + `Alignment.score`. Pure Python, unit-tested.
- **match/regions.py** — per-region V/NDN/J matrices. seqtree `PositionalMatrix` is **fixed-width** but
  CDR3 is variable-length, so v1 uses **Python-side region rescoring**: segment the aligned CDR3 by
  IMGT-style offsets from the conserved C…F/W ends (V-germline prefix, J-germline suffix, NDN middle),
  walk `Alignment.ops` columns, rescore each region with its matrix, sum into S(CDR3). `{region: matrix}`
  config (one matrix may serve all three in v1). Native C++ switching deferred to §G if it's the hot path.
- **match/scoring.py + aggregate/glm.py** — load `vdjam.txt` via `SubstitutionMatrix.from_similarity`
  (penalty `s_aa+s_bb−2s_ab` = legacy "subtract self-score"). Optional cloglog GLM from `score_coef.txt`
  over S(V),S(CDR1),S(CDR2),S(J),S(CDR3),indels. **E-value is primary; GLM is an optional ranking score.**
- **evalue/** — `single.py` wraps `seqtree.evalues` (reuse already-computed target hits; search only the
  control; `exclude_exact=True` for self/DB scans). `control.py` wraps `load_control` (bundled
  `human_trb_aa`; HF `human_tra_aa`/`mouse_*`). `paired.py` (α/β): joint null `π0^{αβ}≈π0^α·π0^β`, paired
  Poisson tail on references matching **both** chains (intersect per-chain `complex.id` hit sets); report
  Fisher-combined p and joint-count p; surface appendix `b2` co-occupancy as a dependence diagnostic.
  `partner_prior.py` (single observed chain weighted by partner rarity: OLGA `Pgen` from samples 4/5, or
  empirical abundance) — pluggable, calibration TBD.
- **aggregate/enrichment.py** — port `ClonotypeSearchSummary`: per (epitope, mhc.class, antigen.species,
  antigen.gene) aggregate matched unique/reads/freq + "not found" complement (polars `group_by().agg()`).
- **cluster/overlap.py** — port `cluster`: `seqtree.pairwise_batch(A,B,params,"aa",threads)` for
  sample-vs-sample fuzzy overlap; emit counts / weighted Jaccard / Morisita.
- **runner/multisample.py** — build VDJdb `Index` once; iterate samples on a shared thread pool; stream
  outputs to bound RAM.
- **cli/** — `match`, `cluster`, `update` (refresh VDJdb cache), `annotate`; flags mirror legacy
  (`--scope s,i,d,t`, `--scoring`, `--matrix`, `--species`, `--gene`, `--match-v/--match-j`, `--paired`,
  `--control`, `--vdjdb-conf`, `--min-epi-size`, `--top-k`, `--threads`, `-o`).
- **Reporting schema** (per-hit rows): query/db ids + cdr3/v/j/locus + epitope/mhc/vdjdb.score +
  n_subs/n_ins/n_dels + alignment_score + cigar + match_string + (glm_score) + informativeness_weight +
  E/p_any/p_enrichment/rule_of_three (+ paired: cdr3a/cdr3b/E_paired/p_paired). Stable column order =
  the future vdjdb-web contract.

---

## C. V-gene + control-normalization workstream (the hard, owner-flagged problem)

A first-class design+experiment track (not a silent default). Three coupled questions, each tested on
epitope data:
1. **Using V + CDR3 jointly** (V+CDR3 is the true motif; V carries biological bias). Options: (a) V as hard
   filter (legacy); (b) V as a soft S(V) term added to score and/or the joint E-value; (c) V-restricted
   neighbour counting with a per-(V-cluster) null. Benchmark which maximizes epitope PR/AUC.
2. **V-gene clustering by CDR1/CDR2/CDR2.5 similarity** (allow TRBV5-1≈TRBV5-8). Extract CDR1/2/2.5 from
   germline (port from `mirpy`, https://github.com/antigenomics/mirpy; its `tcrtrie` is slated to become
   seqtree); compute V-pairwise distances with seqtree; cluster; fuzzy V-match at a tunable cut. Bundle
   precomputed CDR seqs + V distances in `resources/cdr_seqs/`.
3. **Control V-J normalization** (subtle): a real control captures V/J bias automatically **only if its
   V-J usage matches the query repertoire**; ALICE-style P(V,J)-conditioned scaling (`λ=N·Pgen` scaled by
   P_OLGA(V,J)) renormalizes V-J usage but **risks erasing genuine V-gene biological bias** that is part
   of the specificity signal. Plan: implement **both** nulls — (i) matched real control, (ii)
   V-J-conditioned scaled null — and an experiment measuring, on held-out epitopes, whether renormalization
   helps or hurts. Deliver machinery + experiment, not a premature default. V clustering / per-(V,J)
   grouped tries are candidates to port into seqtree (§G).

---

## D. VDJAM re-derivation (research; feeds the manuscript)

A BLOSUM/PAM-analogue for unalignable, over-conserved CDR3s. Extractor + estimator:
- **Observation unit:** same-antigen **Hamming-1** CDR3 pairs from VDJdb (equal length, one substitution),
  per chain (TRA/TRB) and per segment (V/NDN/J via germline-trim offsets).
- **Debias by NDN-origin probability** `w_NDN(p|L,c)` (empirical V/J-trim coverage; OLGA/IGoR posterior as
  robustness) entering **both** foreground and a **position-specific background** `Pr_bg(·|p)` — driving
  log-odds at over-conserved columns (e.g. `CASS`) to ≈0 (solves "A near start ≠ A in middle").
- **Score:** `VDJAM_{c,seg}(a,b)=(1/λ) log2[f/e]` with Laplace pseudocount; load via
  `SubstitutionMatrix.from_similarity`. Over-conserved/contact handling also via
  `PositionalMatrix.from_weights(vdjam, weights)` with weights from `w_NDN` and the Shcherbinin 2023
  contact profile (PMID 37649481: tolerated mismatches are non-contacting / loop-preserving).
- **Leave-one-out-by-epitope:** derive `VDJAM^{(−E*)}`, test VDJdb-vs-VDJdb retrieval on held-out `E*`
  (`exclude_exact=True`) vs BLOSUM62 and unit cost; aggregate per-epitope ΔPR-AUC (Wilcoxon). Ships as a
  reproducible script + `bench/gen_vdjam.py` (analogue of seqtree `bench/gen_pam50.py`).

---

## E. Manuscript (`~/vcs/manuscripts/2026-vdjmatch`, OUP/NAR template, no commits yet)

Full NAR **Article** in `oup-authoring-template/oup-authoring-template.tex`
(`\documentclass[unnumsec,webpdf,contemporary,large]{oup-authoring-template}`, `\journaltitle{Nucleic
Acids Research}`, numbered refs). `git init` + first commit. Structure:
- **Intro** (lift framing from `doc/2018-vdjmatch-old.docx`; refresh to 2025 VDJdb).
- **Results** = 3 deliverables: (1) control-calibrated algorithm + E-value (Fig 1); (2) VDJAM matrix
  (Fig 2: heatmap/MDS/segment matrices/BLOSUM62 corr/indel/contact overlay); (3) VDJdb epitope verification
  PR/ROC/purity/retention stratified human×mouse × MHC-I×II (Fig 3); (4) tool comparison accuracy+speed/RAM
  (Fig 4).
- **Methods**: VDJAM derivation (§D), benchmark protocol (§F), comparison protocol; **transfer/expand the
  relevant seqtree appendix theory** — condense the E-value/Poisson/Chen–Stein/self-match-exclusion
  derivation from `seqtree/appendix/evalue.tex` into NAR Methods, and ship the fuller derivation as
  **Supplementary Note 1** (`evalue.pdf`) or separate supplementary appendices in the manuscript repo.
- **Data/Code availability**, author/funding/COI (owner; declare any ImmuneWatch relationship).
- Lift from 2018: algorithm narrative + VDJAM MDS interpretation + benchmark framing. Write fresh: all
  statistics, segment/chain matrices, contact integration, tool comparison.

### F. Benchmark designs
- **Deliverable 1 (verification, nothing trained):** target=VDJdb, control=OLGA (per species/chain);
  `exclude_exact=True` (punctured null = rigorous form of the legacy informativeness rule). Metrics:
  precision, ROC-AUC, PR-AUC, **purity** (homogeneity of a query's hit-epitope set) and **retention**
  (fraction of labelled clonotypes keeping ≥1 significant same-epitope hit after the E-value filter) —
  exact defs to confirm with owner. Stratify {human,mouse}×{MHC-I,MHC-II}, TRA/TRB; shuffled-label
  negative → AUC→0.5. Reuse `seqtree/bench/bench_epitope.py` machinery.
- **Deliverable 3 (comparison, data-gated):** vs TCRMatch (PMID 33777034, sandbox build+run), ImmuneWatch
  DETECT (commercial → use provided `test_data/immunedetect_results/*.ods`, accuracy only), ERGO/ERGO-II
  (immrep23), immrep25 (owner data). Datasets: samples 1–2 (3 epitopes single-chain), 4–5 (300k OLGA
  random = feasibility frontier), 6 (TCRvdb paired). Metrics: accuracy (PR/ROC) + throughput (TCRs/s) +
  peak RSS vs input size; find where competitors OOM/time out. Mark steps awaiting owner data.

---

## G. Likely seqtree extensions (bump version only when a benchmark justifies it)
- Native **per-region matrix** scoring (variable-width region→matrix map) — replaces §B regions.py rescoring if hot.
- Native **paired** `seed_and_gather` keyed on a paired id (§B paired E-value).
- **V-gene CDR extraction + V-distance/clustering** and per-(V,J) **grouped tries** ported from `mirpy`'s
  tcrtrie into seqtree (§C) — general-purpose, belongs upstream.
Default: pure Python in vdjmatch first; escalate to seqtree C++ only behind a proven gain.

---

## H. Testing & benchmarking
- **Unit:** parsers/alias map, scope→params, CIGAR round-trip, GLM vs reference numbers, E-value formulas
  vs `_poisson_sf`, region segmentation, V-distance, VDJAM estimator on a toy set.
- **Golden:** pinned VDJdb tag; annotate `test_data` samples 1–6; cross-check vs DETECT `.ods` (membership
  exact, scores within tolerance).
- **Property (hypothesis):** dedup idempotent; E monotone in n_control; p_enrichment monotone in n_target;
  paired E ≤ min(E_α,E_β); CIGAR length == alignment length; wider scope never drops a hit.
- **Benchmark harness:** `RUN_BENCHMARK=1`-gated; wall-time + peak RSS; sweep 1k→10k→100k→300k on
  `test_data/sample4/5_olga_airr.txt` (300k) and `sample3_airr.txt` (116k). Targets: annotate 300k vs
  VDJdb in seconds–low-minutes, peak RSS < 8 GB on 32 GB.

---

## I. Verification (end-to-end)
1. `pip install -e .` in repo `.venv` (seqtree from PyPI). `vdjmatch update` caches VDJdb;
   `vdjmatch match --scope 1,0,0,1 test_data/sample1_airr.txt out` → per-hit + summary files.
2. `pytest tests/` green (unit + golden + property); golden diff vs DETECT `.ods` within tolerance.
3. `RUN_BENCHMARK=1 python bench/bench_annotate.py` hits speed/RAM targets on 300k OLGA.
4. Paired: `vdjmatch match --paired test_data/sample6_vdb_airr.txt out` → paired E < min single-chain E for
   true pairs.
5. VDJAM: `python bench/gen_vdjam.py` emits TRA/TRB NDN matrices; LOO-epitope script shows ΔPR-AUC > 0 vs
   BLOSUM62 and unit cost.
6. Manuscript: `latexmk -pdf` compiles the OUP skeleton with Fig/Table placeholders + transferred E-value
   Methods/supplementary.

---

## J. Open questions to iterate with owner (do not block the first milestone)
- **V-gene/control normalization (§C):** central research choice — which V+CDR3 fusion and which null
  (matched real control vs V-J-conditioned scaling) wins on held-out epitopes. Plan delivers machinery +
  experiment; owner decides after seeing results.
- **purity/retention** exact definitions (§F).
- **VDJAM granularity** (per-chain NDN only vs also V/J standalone) and matrix scale λ (legacy integer vs fresh half-bit).
- **Authorship/affiliations/COI/funding**; ImmuneWatch relationship disclosure.
- **immrep25 + extra comparison datasets**, mouse/TRA OLGA controls — owner-provided later.
- **Citations:** verify every DOI via PubMed/arXiv before adding (Shcherbinin 2023 PMID 37649481, TCRMatch
  33777034, ERGO 32983088 / ERGO-II 33981311, VDJdb, GLIPH, TCRdist, OLGA, BLOSUM, AIRR standard, etc.).

---

## Critical files
- seqtree API to build on: `~/vcs/code/seqtree/python/seqtree/__init__.py`, `evalue.py`, `layout.py`; matrix penalty path `src/substitution_matrix.cpp` + `_bindings.cpp`.
- E-value theory to transfer to the paper: `~/vcs/code/seqtree/appendix/evalue.tex`.
- Legacy semantics + resources to port: `~/vcs/code/vdjmatch/src/main/groovy/.../ClonotypeDatabase.groovy`, `src/main/resources/{vdjam.txt,score_coef.txt,segm_score.txt}`.
- Manuscript skeleton: `~/vcs/manuscripts/2026-vdjmatch/oup-authoring-template/oup-authoring-template.tex`; reuse `doc/2018-vdjmatch-old.docx`; data in `test_data/` (+ `immunedetect_results/*.ods`).
- V-gene/CDR + V-J normalization reference to port: https://github.com/antigenomics/mirpy (tcrtrie → seqtree).
