export const meta = {
  name: 'competitor-coverage',
  description: 'tcrdist3 alpha/paired detection + tcrdist3 & GLIPH2 clustering vs vdjmatch',
  phases: [
    { title: 'Run', detail: 'tcrdist3 detection (alpha/paired), tcrdist3 clustering, GLIPH2 clustering' },
    { title: 'Synthesize', detail: 'assemble competitor tables vs vdjmatch' },
  ],
}

const PRE = `
CONTEXT. Extend the vdjmatch benchmark with competitor coverage. Be rigorous; report ONLY numbers.
Repos/envs:
- vdjmatch: /Users/mikesh/vcs/code/vdjmatch  (use ./.venv/bin/python; harness in bench/).
- tcrdist3: conda env "cmp-tcrdist"  (conda run -n cmp-tcrdist python SCRIPT.py ; conda run does NOT
  forward heredoc stdin -> write scripts to /tmp and run by path).
- GLIPH2 binary: /Users/mikesh/vcs/manuscripts/2026-vdjmatch/test_data/gliph2/irtools.osx
- Shared clustering subset (VDJdb public, safe): /tmp/cluster_subset.tsv  (cols: cdr3  v  j  epitope; TRB).
DATA RULE: GLC/YLQ detection queries come from sample6 (TCRvdb) — write temp files only under /tmp,
NEVER under either git repo; report numbers only. Predictions with sample6 sequences (if any) stay in
/tmp or vdjmatch/bench/predictions (already gitignored), never the manuscript repo.

ROC metric (use vdjmatch's, ties grouped):
  import sys; sys.path.insert(0,'/Users/mikesh/vcs/code/vdjmatch/bench'); from metrics import roc_auc, pr_auc_balanced
  roc_auc(list_of_(label0or1, score))   # higher score = more likely positive

vdjmatch numbers to compare against:
  clustering (VDJdb2026 TRB shortlist, n=2270, 11 epitopes): macro_purity=0.989, retention=0.425
  detection GLC: alpha 0.675, beta 0.591, paired 0.682   |  YLQ: alpha 0.865, beta 0.842, paired 0.889
`

const DET_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['tool', 'success', 'results', 'notes'],
  properties: {
    tool: { type: 'string' },
    success: { type: 'boolean' },
    results: { type: 'array', items: {
      type: 'object', additionalProperties: false,
      required: ['task', 'chain', 'roc_auc'],
      properties: { task: {type:'string'}, chain: {type:'string'}, roc_auc: {type:'number'},
                    pr_auc: {type:'number'}, n_pos: {type:'number'}, n_neg: {type:'number'} } } },
    notes: { type: 'string' },
  },
}
const CLUST_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['tool', 'success', 'n_clonotypes', 'n_clusters', 'macro_purity', 'retention', 'notes'],
  properties: {
    tool: { type: 'string' }, success: { type: 'boolean' },
    n_clonotypes: { type: 'number' }, n_clusters: { type: 'number' },
    macro_purity: { type: 'number' }, retention: { type: 'number' },
    params: { type: 'string' }, notes: { type: 'string' },
  },
}

phase('Run')
const tasks = [
  () => agent(`${PRE}
TASK: tcrdist3 ALPHA and PAIRED detection for GLC and YLQ.
The beta-only tcrdist3 detection already exists; you add alpha and paired.
Template (BETA): /Users/mikesh/vcs/code/vdjmatch/bench/_tcrdist_compute.py — single-chain 1-NN distance
of each query to a per-epitope reference. Adapt it for the ALPHA chain (rename to cdr3_a_aa / v_a_gene /
j_a_gene, chains=["alpha"], read trR.rw_alpha) and run for GLC and YLQ.
Data (vdjmatch harness, ./.venv/bin/python):
  from _feat_probe import task_table, ref_table     # sys.path bench+src
  d   = task_table('GLC')      # cols: cdr3(beta),v,j, a_cdr3,a_v,a_j, label(1/0)
  refA= ref_table('GLC','TRA') # alpha reference cdr3,v,j ; refB=ref_table('GLC','TRB')
  # alpha query = (a_cdr3,a_v,a_j); beta query = (cdr3,v,j). Write /tmp ref/query TSVs for tcrdist3.
For each task: compute the 1-NN tcrdist to the ALPHA reference (score = -nearest_alpha_dist, or
radius - dist) and to the BETA reference; report ROC_alpha, ROC_beta, and ROC_paired where paired =
combine the two per-query distances (sum the distances, score = -(dist_alpha+dist_beta)). Score sign:
SMALLER tcrdist = more likely positive, so negate before roc_auc. Drop exact CDR3 self-matches (mirror
vdjmatch exclude_exact). Report ROC for chain in {alpha, beta, paired} for GLC and YLQ. If a gene is
unknown to tcrdist3, drop that clone (note how many). Return DET_SCHEMA.`,
    { label: 'tcrdist:detection', phase: 'Run', schema: DET_SCHEMA, agentType: 'general-purpose' }),

  () => agent(`${PRE}
TASK: tcrdist3 CLUSTERING on /tmp/cluster_subset.tsv (TRB), to compare purity/retention with vdjmatch.
In conda env cmp-tcrdist: build a TCRrep (beta) from the subset (cols cdr3->cdr3_b_aa, v->v_b_gene,
j->j_b_gene, count=1; drop genes tcrdist3 doesn't know), compute the pairwise TCRdist
(tr.compute_distances() gives tr.pw_beta, or compute_sparse_rect_distances self-vs-self within a radius).
Cluster by SINGLE-LINKAGE: connect two clonotypes when tcrdist <= R; union-find -> communities. Sweep R
in {12,24,36,50} and pick the R giving the best purity/retention trade-off (report which R). Using the
KNOWN epitope labels in the subset, compute (same definitions as vdjmatch):
  - macro_purity = size-weighted mean over multi-member clusters of (dominant-epitope count / cluster size)
  - retention = fraction of all clonotypes in a cluster of size>=2 whose dominant-epitope fraction >= 0.7
Return CLUST_SCHEMA (params = chosen R). Compare to vdjmatch 0.989/0.425.`,
    { label: 'tcrdist:clustering', phase: 'Run', schema: CLUST_SCHEMA, agentType: 'general-purpose' }),

  () => agent(`${PRE}
TASK: GLIPH2 CLUSTERING on /tmp/cluster_subset.tsv (TRB), to compare purity/retention with vdjmatch.
GLIPH2 binary: /Users/mikesh/vcs/manuscripts/2026-vdjmatch/test_data/gliph2/irtools.osx (chmod +x first).
Format the input cdr3_file (tab-separated, one row per clonotype):
  CDR3b<TAB>TRBV<TAB>TRBJ<TAB>CDR3a<TAB>subject:condition<TAB>count
Use cdr3/v/j from the subset, CDR3a="NA", subject:condition="s1:c1" (single subject is fine), count=1.
Write a parameter file (see the gliph2 readme.txt in that dir) with out_prefix, cdr3_file=, hla_file=NA,
algorithm=GLIPH2, all_aa_interchangeable=1, local_min_pvalue=0.001, p_depth=1000, kmer_min_depth=3, and
omit the optional refer/v_usage/cdr3_length freq files if GLIPH2 runs without them (else use its bundled
defaults). Run "./irtools.osx -c param_file" from a /tmp working dir. Parse the *_cluster.csv (or the
clone-membership output): map each CDR3b to its GLIPH2 cluster, then using the KNOWN epitope labels from
the subset compute macro_purity and retention with the SAME definitions as the tcrdist agent
(dominant-epitope fraction; pure>=0.7; retention=fraction in pure multi-member clusters). A CDR3 may
appear in multiple GLIPH2 clusters — assign it to its largest/most-significant cluster. Return
CLUST_SCHEMA. If GLIPH2 cannot run (missing ref files), say exactly why in notes and set success=false.`,
    { label: 'gliph2:clustering', phase: 'Run', schema: CLUST_SCHEMA, agentType: 'general-purpose' }),
]
const out = (await parallel(tasks)).filter(Boolean)
log(`competitor runs: ${out.map(o => o.tool + (o.success ? ' ok' : ' FAILED')).join(', ')}`)

phase('Synthesize')
const SYNTH = {
  type: 'object', additionalProperties: false,
  required: ['detection_table', 'clustering_table', 'verdict'],
  properties: {
    detection_table: { type: 'string', description: 'GLC/YLQ alpha/beta/paired ROC: vdjmatch vs tcrdist3' },
    clustering_table: { type: 'string', description: 'purity/retention: vdjmatch vs tcrdist3 vs GLIPH2' },
    verdict: { type: 'string', description: 'where vdjmatch wins/ties/loses; any caveats' },
  },
}
const synth = await agent(`${PRE}
Assemble the competitor-coverage results (JSON below) into two compact comparison tables (detection
paired/alpha; clustering purity/retention) vs vdjmatch, and a short verdict on where vdjmatch wins,
ties, or loses, with caveats. Numbers only.

RESULTS JSON:
${JSON.stringify(out, null, 1)}`, { label: 'synthesize', phase: 'Synthesize', schema: SYNTH })

return { runs: out, synthesis: synth }
