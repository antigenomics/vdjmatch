"""Fresh exhaustive single-chain reference calls and batched peptide-withholding nulls."""
from pathlib import Path


def register(subparsers):
    p = subparsers.add_parser('reference-calls', help='benchmark full-reference winning labels or peptide-withheld nulls')
    p.add_argument('--vdjdb', required=True)
    p.add_argument('--species', default='HomoSapiens')
    p.add_argument('--locus', choices=['TRA', 'TRB'], required=True)
    p.add_argument('--control', help='matching raw repertoire; human TRB has a bundled input')
    p.add_argument('--scope', default='1,0,0,1')
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--calibration-score', choices=['E', 'p_enrichment', 'score'], default='score',
                   help='optional recommendation statistic; never filters continuous outputs')
    p.add_argument('--min-target-records', type=int, default=30,
                   help='benchmark query targets only; the searchable reference remains complete')
    p.add_argument('--leave-epitope-out', action='store_true')
    p.add_argument('--output-prefix', required=True)
    p.set_defaults(func=main)
    return p


def main(args):
    import importlib.metadata
    import json
    import time
    import polars as pl
    from .. import db
    from ..api import Annotator, _prepare, _append_calls
    from ..aggregate.candidates import PMHC
    from ..evalue.control import raw_background
    from ..match.scope import search_params
    from ..match.scoring import load_vdjam, DEFAULT_SCALE

    started = time.monotonic()
    if args.threads < 1:
        raise ValueError('threads must be positive')
    if args.min_target_records < 1:
        raise ValueError('min-target-records must be positive')
    reference = db.load(args.vdjdb, species=args.species, gene=args.locus)
    valid = reference.filter(pl.col('reference_valid'))
    support = valid.group_by('epitope').len().rename({'len': 'reference_records'})
    eligible = support.filter(pl.col('reference_records') >= args.min_target_records)
    benchmark = valid.join(eligible.select('epitope'), on='epitope', how='semi')
    keys = ['cdr3', 'v', 'j']
    if args.leave_epitope_out:
        keys.append('epitope')
    raw = benchmark.select(keys).unique().sort(keys).with_row_index('query_id').with_columns(
        pl.lit(args.locus).alias('locus'))
    ann = Annotator.from_frame(valid)
    _, q = _prepare(raw)
    control, control_provenance = raw_background(args.locus, args.species, args.control)
    exclusions = raw.select('query_id', 'epitope') if args.leave_epitope_out else None
    params = search_params(args.scope, matrix=load_vdjam(), gap_open=DEFAULT_SCALE, gap_extend=DEFAULT_SCALE)
    _, candidates = ann._evidence(q, params, control=control,
        threads=args.threads, exclude_exact=True, exclude_epitopes=exclusions)
    calls = _append_calls(raw, q, candidates, ann.loci, 'vdjmatch_')
    score = pl.when(pl.col('vdjmatch_status') == 'matched').then(
        pl.col('vdjmatch_' + args.calibration_score)).otherwise(None).alias('score')
    if args.leave_epitope_out:
        cases = calls.select('query_id', score,
            pl.lit('null').alias('outcome'), pl.lit('withheld').alias('null_source'))
    else:
        truth = raw.join(valid, on=['cdr3', 'v', 'j'], nulls_equal=True, validate='1:m').select(
            'query_id', *PMHC).unique()
        predicted = calls.select('query_id', *[pl.col('vdjmatch_'+x).alias(x) for x in PMHC])
        correct = predicted.join(truth, on=['query_id', *PMHC], how='semi', nulls_equal=True)['query_id']
        cases = calls.select('query_id', score,
            pl.when(pl.col('query_id').is_in(correct.implode()) & (pl.col('vdjmatch_status') == 'matched'))
            .then(pl.lit('correct')).otherwise(pl.lit('incorrect')).alias('outcome'),
            pl.lit(None, dtype=pl.String).alias('null_source'))
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    calls.write_csv(str(prefix)+'.calls.tsv', separator='\t')
    candidates.write_csv(str(prefix)+'.candidates.tsv', separator='\t')
    cases.write_csv(str(prefix)+'.calibration-cases.tsv', separator='\t')
    support.with_columns((pl.col('reference_records') >= args.min_target_records).alias('benchmark_eligible')).write_csv(
        str(prefix)+'.target-support.tsv', separator='\t')
    raw.join(valid, on=keys, nulls_equal=True, validate='1:m').select('query_id', *PMHC).unique().write_csv(
        str(prefix)+'.truth.tsv', separator='\t')
    provenance = db.provenance(args.vdjdb)
    manifest = {'reference': provenance, 'reference_sha256': provenance.get('input_sha256'),
        'control': control_provenance,
        'software': {x: importlib.metadata.version(x) for x in ['vdjmatch', 'seqtree', 'vdjtools']},
        'options': vars(args) | {'func': 'reference-calls'}, 'source_records': reference.height,
        'valid_reference_records': valid.height, 'invalid_reference_records': reference.height-valid.height,
        'eligible_benchmark_records': benchmark.height, 'ineligible_benchmark_records': valid.height-benchmark.height,
        'eligible_benchmark_peptides': eligible.height, 'all_reference_peptides': support.height,
        'target_support_rule': 'valid source records per species/locus/peptide, before self-exclusion; pool MHC restrictions',
        'search_reference': 'entire valid species/locus reference, including sparse peptides',
        'queries': raw.height, 'candidate_rows': candidates.height, 'wall_seconds': time.monotonic()-started,
        'self_exclusion': 'benchmark-only full-junction puncture of hits and population sizes',
        'observation_unit': 'distinct full-junction/V/J'+('/truth peptide' if args.leave_epitope_out else ''),
        'peptide_withholding': 'all restrictions removed before reduction' if args.leave_epitope_out else None,
        'scoring': 'production default vdjam matrix, score_scale400, soft V weight0.25',
        'limitations': 'Fixed-ball enrichment is evidence, not a posterior for the winning label. No paired calibration.'}
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return 0
