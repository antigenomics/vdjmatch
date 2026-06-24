export const meta = {
  name: 'glc-lll-feature-discovery',
  description: 'Exhaustively probe V/J/length/pairing/structure features to beat GLC and LLL detection',
  phases: [
    { title: 'Probe', detail: 'one agent per feature axis; ROC-alone + combined for GLC/LLL/YLQ' },
    { title: 'Combine', detail: 'build best combined score from discriminative features' },
    { title: 'Verify', detail: 'adversarial: label-shuffle null + no-harm on NLV/YLQ' },
  ],
}

// ---- shared preamble every probe agent receives (agents do NOT see the parent conversation) -----
const PRE = `
CONTEXT. You are probing TCR-epitope DETECTION features in the vdjmatch repo to beat two HARD tasks.
Working dir: /Users/mikesh/vcs/code/vdjmatch  (use ./.venv/bin/python ALWAYS).
ABSOLUTE DATA RULE: read the manuscript test_data ONLY at runtime via the harness below; NEVER copy
any sequence data into the repo or into your answer — report ONLY aggregate numbers/distributions.

TASKS (all TRB unless noted), with the score to BEAT and the current vdjmatch baseline:
  GLC (GLCTLVAML, A*02): pos=107 padj<1e-5, neg=42 padj>=1e-5. BEAT imw-detect ROC 0.679. vdjmatch now 0.591.
        -> GLC is NOT separable by CDR3 distance: pos & neg are EQUIDISTANT from the dense ref (6592).
           tcrdist also = 0.592. The signal must come from germline/structure features, NOT CDR3 distance.
  LLL (LLLGIGILV, A*02): pos=232, neg=86 (LLW as control). BEAT tcrdist ROC 0.615. vdjmatch now 0.599.
           LLL reference is SPARSE (233 records) — few neighbours.
  YLQ (YLQPRTFLL, A*02): pos=149 neg=254. Already won (vdjmatch 0.842). Use as a SANITY check that a
           feature does not HURT an already-separable task.

HARNESS (bench/_feat_probe.py — import it; it reads test_data at runtime, nothing copied):
  import sys; sys.path.insert(0,'bench'); sys.path.insert(0,'src')
  from _feat_probe import task_table, ref_table, baseline_scores
  from metrics import roc_auc           # roc_auc(list_of_(label,score)); ties grouped (0.5 credit)
  import benchmark as B
  d   = task_table('GLC')               # polars df: query_id,label(1/0),cdr3,v,j,length [,a_cdr3,a_v,a_j]
  ref = ref_table('GLC')                # epitope's A*02 ref: cdr3,v,j (allele-stripped gene names)
  base= baseline_scores('GLC')          # {cdr3 -> vdjmatch caldens score}  (the committed detection score)
  # background OLGA repertoire for V/J/length priors:
  from vdjmatch.evalue import background; ctrl_trb = background('TRB')   # seqtree Index of OLGA TRB
  # NOTE for a V/J/length BACKGROUND distribution use sample4 (OLGA TRB) raw if needed:
  #   import polars as pl; from compare import TESTDATA
  #   olga = pl.read_csv(TESTDATA/'sample4_olga_airr.txt', separator='\\t')  # AIRR: v_call,j_call,junction_aa
  v25 = None  # if you need J/V of the reference it is already in ref_table (ref['j'], ref['v']).

ROC CONVENTION: higher score = more likely POSITIVE. roc_auc>0.5 means the feature ranks pos above neg.
If a feature is anti-correlated, report |its| direction. Always compute on GLC and LLL (and YLQ sanity).
Combined score: try base[cdr3] * exp(feature_logodds) or base + lambda*feature; report best combined ROC.

Return ONLY the structured object. Keep any text you write free of raw sequences.`

const PROBE_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['feature','glc_roc_alone','lll_roc_alone','ylq_roc_alone','glc_combined','lll_combined',
             'discriminative','signal','recommend'],
  properties: {
    feature: { type: 'string' },
    glc_roc_alone: { type: 'number', description: 'ROC of this feature alone on GLC TRB' },
    lll_roc_alone: { type: 'number' },
    ylq_roc_alone: { type: 'number' },
    glc_combined: { type: 'number', description: 'best ROC of baseline caldens COMBINED with this feature, GLC' },
    lll_combined: { type: 'number' },
    ylq_combined: { type: 'number', description: 'YLQ combined (sanity: must not drop below 0.84)' },
    discriminative: { type: 'boolean', description: 'true if it adds real signal on GLC or LLL beyond noise' },
    signal: { type: 'string', description: '1-3 sentences: how pos vs neg differ (numbers only, no sequences)' },
    how_combined: { type: 'string', description: 'exact formula used to combine with baseline' },
    recommend: { type: 'string', description: 'wire it? how? caveats?' },
  },
}

const FEATURES = [
  { key: 'length', task: 'Probe', prompt: `${PRE}
FEATURE = CDR3 LENGTH PRIOR. Hypothesis: true binders have epitope-typical CDR3 length; bystanders do not.
Build P(L | reference) from ref['length'] (=ref cdr3 len) smoothed, and P(L | background OLGA) from sample4.
Score each query by log[P(L|ref)/P(L|bg)]. Report ROC alone + combined (base * exp(score)). Also just test
raw length and |length - ref_modal_length| as discriminators. Cover GLC, LLL, YLQ.` },
  { key: 'jgene', task: 'Probe', prompt: `${PRE}
FEATURE = J-GENE PRIOR. Hypothesis: GLC/LLL binders are J-restricted; bystanders use random J.
Build P(J | reference) from ref['j'] and P(J | background) from sample4 j_call (OLGA). Score query by
log[P(J|ref)/P(J|bg)] (add-0.5 smoothing). ROC alone + combined. Report which J genes are enriched in
pos vs neg (gene names + frequencies only). Cover GLC, LLL, YLQ.` },
  { key: 'vgene', task: 'Probe', prompt: `${PRE}
FEATURE = V-GENE PRIOR + V-LOOP. Two parts: (a) V-usage prior log[P(V|ref)/P(V|bg)] like the J probe;
(b) germline CDR1+CDR2 loop similarity of the query V to the reference V set — use
    from vdjmatch.match import vgene as vg ; vg.vsim(qV, refV) ; vg.load_v_regions() gives FR1-3+CDR1/2.
Score query by max/mean vsim to the reference V's, and by the V-usage prior. ROC alone + combined.
Report top enriched V genes in pos vs neg (names+freqs). Cover GLC, LLL, YLQ.` },
  { key: 'vjl_prior', task: 'Probe', prompt: `${PRE}
FEATURE = JOINT V+J+LENGTH GERMLINE PRIOR (the KEY hypothesis for GLC). Model the reference germline
configuration: score = log P(V|ref)+log P(J|ref)+log P(L|ref) - [same under OLGA background]. This asks
"is this TCR's germline (V,J,length) typical of the epitope's binders?" — orthogonal to CDR3 distance.
Use add-0.5 smoothing; background from sample4 (OLGA TRB: v_call,j_call,junction_aa length). Also test the
PRODUCT with baseline caldens: final = base * exp(beta*prior), sweep beta in {0.5,1,2}. This is the most
important probe — be thorough. ROC alone + combined on GLC, LLL, YLQ. Report each component's individual ROC.` },
  { key: 'pairing', task: 'Probe', prompt: `${PRE}
FEATURE = ALPHA/BETA PAIRING (GLC & YLQ only; LLL is TRB-only — report lll as null/0.5). sample6 is paired:
task_table('GLC') has a_cdr3,a_v,a_j (the paired ALPHA chain of each beta query). Compute a TRA baseline
  baseTRA = baseline_scores('GLC','TRA')   # {alpha_cdr3 -> caldens TRA score}
then per query combine beta score base[cdr3] with alpha score baseTRA[a_cdr3] (lookup; missing->0):
  test sum, product, max, and geometric mean of the two; also (beta_rank+alpha_rank). Report best combined
ROC for GLC and YLQ vs the TRB-only baseline (GLC 0.591, YLQ 0.842). Does pairing lift GLC above 0.679?` },
  { key: 'physchem', task: 'Probe', prompt: `${PRE}
FEATURE = CDR3 PHYSICOCHEMISTRY at the hypervariable apex. For the middle CDR3 residues (positions ~5..L-4,
the loop apex that contacts peptide), compute mean charge, hydrophobicity (Kyte-Doolittle), volume, and
aromatic content. Compare pos vs neg distributions; build a logistic-free linear discriminant (just
z-score each, sum the ones that separate) and report ROC alone + combined. Cover GLC, LLL, YLQ.` },
  { key: 'structure', task: 'Probe', prompt: `${PRE}
FEATURE = STRUCTURAL CONTACT POSITIONS (GLC focus; try LLL if a template exists). Goal: find WHICH CDR3beta
positions contact the GLC peptide, then weight CDR3 distance by those positions.
STEPS:
 1. Find a GLCTLVAML/HLA-A*02 TCR-pMHC structure in /Users/mikesh/vcs/code/tcren-ms/data/Native2026/*.pdb.gz
    (zcat | grep for the peptide 'GLCTLVAML' in SEQRES, or check a known PDB list). Report the PDB id(s).
 2. Preferred: use the tcren env to get the CDR3beta<->peptide contact map:
      cd /Users/mikesh/vcs/code/tcren-ms ; conda run -n tcren python -c "..."
      from tcren.structure.io import import_structure; from tcren.annotation import classify_chains
      from tcren.mhc import annotate_mhc; from tcren.project2d import residue_markup_table, region_pair_contacts
      -> region_pair_contacts(s, kind='closest') filtered to trb/cdr3 <-> peptide; gives the contacting
      CDR3beta residue indices. Convert to RELATIVE positions within the CDR3 (fraction along the loop).
    Budget this to ~10 min; if the tcren env is unavailable/slow, FALL BACK to the canonical assumption
    (CDR3 apex positions 5..L-4 contact peptide) and SAY you used the fallback.
 3. Build a position-weight vector over CDR3 (contacts weighted high) and re-score GLC CDR3 distance with it
    (substitutions at contact positions cost more). Test ROC alone + combined vs baseline.
Report the contact positions found (the reusable structural fact) and whether contact-weighting helps GLC/LLL.` },
]

phase('Probe')
const probes = await parallel(FEATURES.map(f => () =>
  agent(f.prompt, { label: `probe:${f.key}`, phase: 'Probe', schema: PROBE_SCHEMA })
))
const found = probes.filter(Boolean)
log(`probes done: ${found.length}/${FEATURES.length}; discriminative=${found.filter(p=>p.discriminative).map(p=>p.feature).join(', ')}`)

phase('Combine')
const COMBINE_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['glc_best_roc','lll_best_roc','beats_glc_target','beats_lll_target','recipe','code_sketch','no_harm'],
  properties: {
    glc_best_roc: { type: 'number' },
    lll_best_roc: { type: 'number' },
    beats_glc_target: { type: 'boolean', description: 'glc_best_roc >= 0.679' },
    beats_lll_target: { type: 'boolean', description: 'lll_best_roc >= 0.615' },
    recipe: { type: 'string', description: 'the exact combined score that wins, in words' },
    code_sketch: { type: 'string', description: 'minimal python to reproduce the winning GLC & LLL ROC' },
    no_harm: { type: 'string', description: 'NLV and YLQ ROC under the same combined recipe (must stay >= baseline)' },
    notes: { type: 'string' },
  },
}
const combine = await agent(`${PRE}
You are GIVEN the single-feature probe results (JSON). Build the BEST COMBINED detection score for GLC and
LLL using the discriminative features, starting from the baseline caldens (baseline_scores). Actually RUN
the combination in python and measure ROC. Targets: GLC>=0.679, LLL>=0.615. CRUCIAL: also measure NLV and
YLQ under the SAME recipe — the combination must NOT drop NLV below 0.624 or YLQ below 0.842 (no per-task
cherry-picking; one unified recipe). Prefer principled products of likelihood-ratios (germline prior) with
the caldens distance score. Report the winning recipe, code, and the no-harm NLV/YLQ numbers.

PROBE RESULTS JSON:
${JSON.stringify(found, null, 1)}`,
  { label: 'combine', phase: 'Combine', schema: COMBINE_SCHEMA })

phase('Verify')
const VERIFY_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['holds','shuffle_roc','verdict'],
  properties: {
    holds: { type: 'boolean', description: 'true if the recipe survives the null + no-harm checks' },
    shuffle_roc: { type: 'number', description: 'GLC ROC of the winning recipe under a label-shuffle (should be ~0.5)' },
    verdict: { type: 'string', description: 'is the GLC/LLL gain real and fair, or an artifact? be adversarial' },
    fairness: { type: 'string', description: 'does the recipe look principled & unquestionable for a paper?' },
  },
}
const verify = await agent(`${PRE}
ADVERSARIALLY VERIFY this winning combined recipe for GLC/LLL. Try to REFUTE it:
 1. Label-shuffle null: randomly permute GLC labels (fixed seed via index parity, NOT Math.random) and
    recompute the recipe ROC — it MUST collapse to ~0.5. If it stays high, the recipe is leaking/overfit.
 2. No-harm: re-confirm NLV>=0.624 and YLQ>=0.842 under the recipe.
 3. Principle check: is each feature biologically justifiable (germline restriction, structural contacts)
    and free of test-label leakage? Would a skeptical NAR reviewer accept it as fair?
Default to holds=false if anything is uncertain.

WINNING RECIPE:
${JSON.stringify(combine, null, 1)}`,
  { label: 'verify', phase: 'Verify', schema: VERIFY_SCHEMA })

return { probes: found, combine, verify }
