# vdjmatch <-> vdjtools: division of labour, and what vdjtools 3.8 adds

Written 2026-08-14, from the vdjtools side, while `vdjtools.motif` was being built. Nothing here
changes vdjmatch; it records what is now available upstream, what vdjmatch should **not** take on,
and the one problem that genuinely belongs here rather than there.

## The line

**vdjmatch searches across samples and against databases. vdjtools discovers and describes within
one sample.** That split is not bureaucratic — the two need different objects:

| | vdjmatch | vdjtools |
|---|---|---|
| question | "what does this clonotype match?" | "what is convergent in this repertoire?" |
| reference | VDJdb, another sample, a prototype library | the sample's own neighbourhood |
| object | a hit list with an E-value | a metaclonotype and its PWM |
| output | `vdjmatch_{epitope,score,n_hits}` | `members` / `summary` / logo |

vdjtools 3.8 delegates *every* sequence operation to vdjmatch/seqtree and reimplements none of it:
`vdjmatch.cluster.overlap` builds every graph edge, `vdjmatch.evalue.query_evalues` computes every
enrichment, `vdjmatch.db.load` reads VDJdb, `seqtree.Index` does every neighbour search, and
`seqtree.embed_in_frame` defines the column frame. If any of that starts being duplicated in
vdjtools, that is a bug in vdjtools.

## What landed in vdjtools 3.8 that vdjmatch can use

### `vdjtools.overlap.tcrnet` grew a `grouping=` and a reusable control

Two defects were fixed that also affect anything built on `evalue.query_evalues`:

1. **Grouping.** `E = (N/M)*n_control` scales the control degree by the ratio of index sizes, so
   both sides must be restricted to the same clonotype group. Under a substitution-only scope the
   correct scale is `N_L/M_L`, and that ratio moves **2.2-fold** between junction length 11 and 18;
   ungrouped, a "core set" is substantially a short-CDR3 selection. Measured over 149,224 real
   Hamming-1 edges, neighbours are **23.5x** enriched for sharing the same (V, J). `tcrnet` now
   defaults to `grouping="vj"` with backoff `vjl -> vj -> v -> l -> none` and reports which level it
   used per row.
2. **The empty control ball.** With `n_control == 0` the old code produced `E = 0` and
   `poisson.sf(k-1, 0) == 0` exactly, so every such clonotype came back maximally significant.
   `E` is now floored at the rule-of-three bound `3*N/M`, matching what `seqtree.evalues` already
   does internally. **If `vdjmatch.evalue.single.query_evalues` callers recompute `E` themselves,
   they have the same hole** — worth a check; `rule_of_three` is already in the output frame but is
   easy to ignore.

`vdjtools.overlap.tcrnet.build_grouped_control(control_frame, grouping="vj")` returns a
`GroupedControl` that builds one `seqtree.Index` per group **lazily and once**, so a 287-sample
cohort indexes its control a single time. If vdjmatch ever scores a cohort against a partitioned
background, this is the shape to copy (or to move down here and have vdjtools import).

### The control repertoire vdjmatch's `evalue.background` cannot reach

`seqtree.control.load_control` returns a sequence-only `Index`, so no V/J-matched cell is reachable
through it. The V/J-annotated build is
`isalgo/airr_control` → `vdjdb-web-control/{HomoSapiens,MusMusculus}.{TRA,TRB}.txt.gz`: 1e6 uniform
random productive clonotypes per species/chain carrying `v`, `j`, `VEnd`, `JStart`. vdjtools reads
it with `vdjtools.motif.pipeline.load_control(species=, locus=)`. Provenance is recorded in
vdjtools' `SOURCES.md`; **seqtree's own `SOURCES.md` is stale** — it names
`antigenomics/vdjdb-cdr3-frequencies` and says the cache is keyed by name, and neither is true of
`seqtree/control.py` any more.

### VDJdb's shipped motif artifacts, and a defect in them

`vdjtools.motif.vdjdb` reads `motif_pwms.txt` and `cluster_members.txt` from the release ZIP.
Three things worth knowing if vdjmatch ever surfaces these (its `Motifs.scala` consumer does):

- **`I.norm` is a cross-entropy, not a divergence.** `H(f,q) = H(f) + D_KL(f||q)` and the
  observed-entropy term is never subtracted, so `I.norm` scores high whenever letters are merely
  *rare* in the background. Recover the divergence as `D_KL [bits] = log2(20)*(2*I.norm + I - 1)`;
  verified on all 13 columns of `H.B.GILGFVFTL.1`.
- **`freq` does not sum to 1 per `(cid, pos)`** — the R generator inner-joins the foreground onto its
  background pool, dropping zero-background cells from 9.68% of columns. Do **not** renormalise: the
  dropped cells are the maximally enriched letters.
- **`v.end` / `j.start` in `cluster_members.txt` are amino-acid positions**, unlike the nucleotide
  `VEnd` / `JStart` of the control repertoires. Same names, different units, one file apart.

### `aggregate.enrichment` counts but does not test — vdjtools 3.8 supplies the test

`vdjmatch.aggregate.epitope_summary` tabulates unique matched query clonotypes, reads and best score
per `(epitope, mhc_class, antigen_species)`, and stops there. That table is not yet a result: public,
high-Pgen junctions match VDJdb in every repertoire, and the epitopes with the most records collect
the most matches for that reason alone. `vdjtools.motif.enrich` adds the missing layer and is happy
to consume a vdjmatch hit list (`unit=` and `by=` are plain column names).

The design decision worth copying, and the one that is easy to get backwards:

- **The null is a control repertoire annotated identically, not the database's own composition.**
  `expected = (records for epitope e / all records) * matches` looks obvious and is wrong.
  `vdjdb-web` documents why at `ControlRepertoires.scala:40-46`: reachability and record count are
  different things — `YLEPGPVTA` holds 11 records and is reached by far more of a random repertoire
  than `EPLPQGQLTAY` with 39, because its motif sits where V(D)J recombination puts mass.
- Run the identical search, with the identical scope and the identical HLA/confidence restriction,
  against a control repertoire; compare compositions **among matched clonotypes** on both sides, so
  depth and the overall annotation rate cancel. Beta-Binomial upper tail with the control's counts
  as a Jeffreys prior (`alpha = control_hit + 0.5`, `beta = max(0, control_matched - control_hit)
  + 0.5`), because the control proportion is itself an estimate.
- **Restricting VDJdb to `cluster_members.txt` sharpens it.** Measured on yellow-fever day-15
  paratope motifs: against all records LLWNGPMAV is 10/38 vs 479/10,743 control, log2 OR +2.98,
  q = 8.0e-5; against the motif clusters, 10/11 vs 192/2,740, log2 OR +6.53, q = 7.3e-11.
- **A convergent family is one observation.** `collapse_families` + `unit="family_id"` applied to
  both sides is the conservative read.

**One guard to copy verbatim.** `vdjmatch.db.load(species="human")` returns **zero rows** — it wants
`"HomoSapiens"` — and every count downstream is then a well-formed zero that reads as a clean null.
vdjtools' `annotate` now raises on an empty database and names the mistake in the message. `db.load`
itself would be a better place for that check.

HLA matching (`parse_hla` / `hla_matches` / `hla_expr`) is ported from `vdjdb-web`'s
`HlaAllele.scala`: two-field resolution, prefix compatibility in both directions, free-text split,
an **OR across `mhc.a` and `mhc.b`** (no single-column filter expresses it), and non-HLA cells
(`B2M`, `H-2Db`) never matching — without that last rule every class-I hit is "HLA-matched" via
`B2M`. If vdjmatch grows an `--hla` filter, this is the semantics to match.

## The open problem, and why it belongs here

**Aligning omega loops of unequal length so a shared motif stays in one column.**

`seqtree.embed_in_frame(seq, width, c)` puts one contiguous gap block at frame column `c`. vdjtools
now sweeps `c` over its whole range and picks the minimum-column-entropy frame
(`vdjtools.motif.pwm.frame_scan`). Measured on synthetic metaclonotypes over lengths 12-16, 300
members, with a two-residue motif that cannot arise by chance:

| motif placement | argmin `c` | members with the motif in one column |
|---|---|---|
| V-anchored at offset 5 | 7 | **100%** at `c*`; **20%** at the fixed `c = 3` or `4` |
| J-anchored, 6 from the end | 5 | 100% at `c*`; also 100% at `c = 3` or `4` |
| offset varies per member | 6 | 30.7% at best — no single-block frame aligns it |

So the three rules one might fix in advance — "open the block after residue 3-4", "before the last
3-4", "pad the middle" — are each a **single point of this curve**, and choosing one in advance is
exactly what smears a V-anchored motif across five columns. Sweeping recovers it, because the
minimum lands *after* the motif; that is the same statement as "the motif folds to the left".

That is as far as a **display** frame needs to go, and it is where vdjtools stops. What is left is
the **search** problem, and it is yours:

1. **A frame is not an aligner.** vdjtools picks one frame per metaclonotype from members it already
   has. Searching a repertoire against a prototype needs a *scoring* alignment for a query of
   arbitrary length, which is `seqtree.gapblock` territory: `gapblock_score`, `GapBlockIndex`, and
   the `frame_prior` / `positions_prior` / `central_prior` family. seqtree's own measurement is that
   free placement is worse than a hard pin (precision 0.156 vs 0.414 at matched FPR), which is the
   search-side analogue of the table above.
2. **The unsolved case is the third row, and a better `c` cannot fix it.** A motif at a
   member-varying offset carries no column-frame information at all — under any single-block frame
   it is indistinguishable from no motif (`frame_depth` 0.232 vs 0.215). Representing it needs a
   different object: a profile HMM with a real insert state, or a gapped-k-mer representation. That
   is a seqtree engine question, not a vdjtools plotting question.
3. **Two blocks are cheap to define and hard to justify.** Any deterministic length -> gap-columns
   rule gives a consistent column index (the frame *is* the alignment), so multi-block frames are
   not ruled out by transitivity — seqtree's transitivity argument is about reproducing the
   *pairwise-optimal single-block* alignment, which a display frame does not need. The real question
   is whether a second block ever buys enough conserved columns to pay for the search cost. That is
   measurable with `frame_scan` extended to a 2-D `(c1, c2)` sweep, on real VDJdb clusters with
   length spread. Nobody has measured it.

**`IslandProfile` is the natural home for the result.** It already is a per-island PWM over a fixed
frame with a calibrated cutoff (`thetas_from_scores`), measured at +10.9 pp recall over
min-over-members in the repertoire-annotation regime. A vdjtools metaclonotype is exactly an island;
if vdjmatch wants to annotate a repertoire against a library of discovered motifs, the path is
`vdjtools.motif.discover` -> members -> `seqtree.IslandProfile.fit` -> `vdjmatch` search, and the
only missing piece is the library format. vdjtools deliberately did **not** build
`vdjtools motif annotate` for this reason.

## What vdjtools will not do

Cross-sample and cross-database search stays here. vdjtools has no annotator, no epitope caller, no
prototype library, and no paired-chain matcher, and should not grow them — it calls
`vdjmatch.db.load` and `vdjmatch.cluster.overlap` when it needs either.
