export const meta = {
  name: 'epitope-structure-hotspots',
  description: 'Per-epitope TCR-pMHC contact hotspots + structure-informed CDR3 motifs (NLV/YLQ/LLW/LLL/GLC)',
  phases: [
    { title: 'Structure', detail: 'one agent per epitope: contacts + hotspots + structure-informed motif' },
    { title: 'Synthesize', detail: 'cross-epitope pattern, easy/hard heuristic, detection recipe' },
  ],
}

const PRE = `
You are extracting TCR-pMHC structural CONTACT HOTSPOTS and CDR3 MOTIFS for a TCR-epitope detection paper.
Two repos / two envs (orchestrate both from bash):
  STRUCTURE: /Users/mikesh/vcs/code/tcren-ms  -- run via:  conda run -n tcren python ...
  REFERENCE: /Users/mikesh/vcs/code/vdjmatch   -- run via:  ./.venv/bin/python ...
DATA RULE: read structures + the manuscript reference at runtime; report ONLY positions/motifs/numbers,
never dump large sequence lists.

STRUCTURE RECIPE (proven; tcren env). For a given PDB id (file data/Native2026/<id>.pdb.gz):
  cd /Users/mikesh/vcs/code/tcren-ms
  conda run -n tcren python - <<'PY'
  from pathlib import Path
  from tcren.structure.io import import_structure
  from tcren.annotation import classify_chains
  from tcren.mhc import annotate_mhc
  from tcren.project2d import residue_markup_table, region_pair_contacts
  s = import_structure(Path('data/Native2026/3o4l.pdb.gz'), pdb_id='3o4l')
  classify_chains(s, organism='human'); annotate_mhc(s)
  markup = residue_markup_table(s)   # rows: complex_chain (tra/trb/peptide/mhca..), complex_region (cdr1/2/3,fr..),
                                     #       aa, aa_index (position within region), residue_index
  rp = region_pair_contacts(s, kind='closest')   # contacting residue pairs across chains/regions
  # Keep CDR3<->peptide pairs (both orientations): complex_chain in {tra,trb} & region==cdr3 vs peptide.
  # For each contact, record the CDR3 residue's aa_index (position within CDR3) and the peptide position.
  PY
  If tcren env fails, fall back to Biopython (already in tcren-ms env or vdjmatch .venv): parse the .pdb.gz,
  the peptide is the short 8-11mer chain; the two TCR V-domains are the chains with a Cys...FGXG CDR3 loop;
  compute heavy-atom contacts < 4.5 A between CDR3 residues and peptide residues. SAY if you used fallback.

CDR3 POSITION CONVENTION: report contacts as RELATIVE position within the CDR3 (e.g. "residues 5-9 of a
13-mer", or as fraction 0..1 along the loop, and IMGT number if tcren gives it). The CDR3 apex (middle)
typically contacts the peptide; the Cys/FGXG ends are germline.

STRUCTURE-INFORMED MOTIF (combine structure + reference): once you know WHICH CDR3 positions contact the
peptide, look at the epitope's reference CDR3 set to see WHAT residues sit there:
  cd /Users/mikesh/vcs/code/vdjmatch
  ./.venv/bin/python - <<'PY'
  import sys; sys.path.insert(0,'bench'); sys.path.insert(0,'src')
  from _feat_probe import ref_table
  r = ref_table('GLC')   # cdr3,v,j of the epitope's A*02 reference (TRB). len(r) rows.
  # For the structurally-contacting positions, tabulate residue frequencies across reference CDR3s
  # (align by the apex: e.g. center the CDR3, or use position-from-Cys and position-from-FGXG).
  PY
Report: is there a FOCUSED motif (a few conserved residues at the hotspot = easy to detect) or is it
DEGENERATE (no conservation at contacts = hard, e.g. GLC)? This is the easy/hard heuristic.`

const STRUCT_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['epitope','peptide','pdb_ids','cdr3b_contacts','peptide_hotspots','structural_motif',
             'motif_focus','easy_or_hard','detection_use'],
  properties: {
    epitope: { type: 'string' },
    peptide: { type: 'string' },
    pdb_ids: { type: 'array', items: { type: 'string' }, description: 'structures actually analyzed (peptide-verified)' },
    used_fallback: { type: 'boolean' },
    cdr3b_contacts: { type: 'string', description: 'CDR3beta positions contacting peptide (relative + IMGT if avail)' },
    cdr3a_contacts: { type: 'string', description: 'CDR3alpha positions contacting peptide, or "n/a"' },
    peptide_hotspots: { type: 'string', description: 'which peptide positions (P1..P9) the TCR contacts most' },
    structural_motif: { type: 'string', description: 'residues conserved at the contacting CDR3 positions in the reference (numbers/letters only)' },
    motif_focus: { type: 'string', enum: ['focused','partial','degenerate'], description: 'conservation at the hotspot' },
    easy_or_hard: { type: 'string', enum: ['easy','medium','hard'] },
    detection_use: { type: 'string', description: 'how to use these contacts/motif to score detection for THIS epitope' },
    notes: { type: 'string' },
  },
}

const EPIS = [
  { tag: 'NLV', pep: 'NLVPMVATV', pdb: ['3gsn','5d2l','5d2n'] },
  { tag: 'YLQ', pep: 'YLQPRTFLL', pdb: ['7n1f','7n6e','7pbe','7rtr'] },
  { tag: 'GLC', pep: 'GLCTLVAML', pdb: ['3o4l'] },
  { tag: 'LLL', pep: 'LLLGIGILV', pdb: ['7q9a'] },
  { tag: 'LLW', pep: 'LLWNGPMAV', pdb: [] },
]

phase('Structure')
const results = await parallel(EPIS.map(e => () => agent(`${PRE}

YOUR EPITOPE: ${e.tag} = ${e.pep} (HLA-A*02:01).
Structures found in tcren Native2026 (peptide sequence already verified in ATOM records): ${e.pdb.length ? e.pdb.join(', ') : 'NONE'}.
${e.pdb.length ? 'Analyze each (or a representative if several agree). Report consensus contacts/hotspots.' :
  'NO structure in Native2026. Search RCSB for a '+e.pep+'/HLA-A2 TCR-pMHC complex (verify the peptide '+
  'sequence before trusting any PDB id; do NOT invent a PDB id from memory). If none exists, set pdb_ids=[] '+
  'and motif_focus from the reference CDR3 conservation alone; say no structure was available.'}
Extract CDR3beta (and CDR3alpha) residues contacting the peptide, the peptide hotspot positions, and the
structure-informed motif (reference conservation at the contacting positions). Use ref_table('${e.tag}').
Return the structured object only.`,
  { label: `struct:${e.tag}`, phase: 'Structure', schema: STRUCT_SCHEMA })))

const ok = results.filter(Boolean)
log(`structures done: ${ok.map(r=>r.epitope+'='+r.motif_focus).join(', ')}`)

phase('Synthesize')
const SYNTH_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['contact_pattern','easy_hard_rule','glc_recipe','lll_recipe','detection_recommendation'],
  properties: {
    contact_pattern: { type: 'string', description: 'the shared CDR3<->peptide contact geometry across the 5 epitopes' },
    easy_hard_rule: { type: 'string', description: 'a concrete heuristic: which structural/motif feature predicts an easy vs hard epitope' },
    glc_recipe: { type: 'string', description: 'specific structure/motif-based way to lift GLC detection above 0.679' },
    lll_recipe: { type: 'string', description: 'specific way to lift LLL detection above 0.615' },
    detection_recommendation: { type: 'string', description: 'how to fold contact-position weighting + structural motif into vdjmatch scoring, fairly' },
    per_epitope_table: { type: 'string', description: 'compact table: epitope | hotspot CDR3 positions | peptide hotspots | motif | easy/hard' },
  },
}
const synth = await agent(`${PRE}

You are GIVEN the per-epitope structural results (JSON). Synthesize for the detection paper:
 1. The shared contact geometry (which CDR3 positions / peptide positions dominate across NLV/YLQ/LLW/LLL/GLC).
 2. A concrete EASY-vs-HARD heuristic: what structural/motif property predicts an epitope is easy (e.g. focused
    CDR3 motif at a few hotspots, like the GIL RS motif / YLQ) vs hard (degenerate contacts, like GLC)?
 3. A specific, FAIR (no label leakage) recipe to lift GLC above 0.679 and LLL above 0.615 using contact-position
    weighting and/or the structure-informed motif PSSM. Be honest if structure does NOT help a given task.
 4. How to fold this into vdjmatch scoring as a clean, unquestionable feature.

PER-EPITOPE RESULTS JSON:
${JSON.stringify(ok, null, 1)}`,
  { label: 'synthesize', phase: 'Synthesize', schema: SYNTH_SCHEMA })

return { per_epitope: ok, synthesis: synth }
