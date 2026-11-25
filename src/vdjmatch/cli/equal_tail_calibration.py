"""Fit or apply empirical any-call calibration for the complete equal-tail search."""
import json
from importlib.resources import files
from pathlib import Path

import polars as pl

from .. import db
from ..db.cache import sha256
from ..evalue.equal_tail_selection import equal_tail_selection
from ..evalue.paired import build_paired_ref
from .calibration import apply_calibration, null_tail_probabilities, operating_points
from .equal_tail_compatibility import compatible_binding


def null_path(recipe):
    path = Path(recipe)
    if not path.name.endswith('.calibration.json'):
        raise ValueError('equal-tail recipe must end in .calibration.json')
    return path.with_name(path.name.removesuffix('.calibration.json') + '.null-cases.tsv')


def _binding(args, species):
    from seqtree import _core, gapblock
    from ..evalue import directional, equal_tail, equal_tail_selection as selection, kernel, paired
    from ..match import geometry
    from ..aggregate import candidates
    from ..db import schema
    from ..io import airr
    from . import calibration, historical, equal_tail as input_owner
    modules = [equal_tail, selection, kernel, geometry, paired, candidates, schema, historical, input_owner, airr, calibration, directional, gapblock, _core]
    return dict(event='maximum supported equal-tail score over the complete species TRA/TRB compatible-pMHC catalogue',
                axes=['paired','TRA','TRB'] if args.locus=='paired' else [args.locus],
                paired='same-reference linked pair; full marginal references, including unlinked observations',
                ties='exact ratio product_numerator/tail_denominator**2; distinct tied pMHCs ambiguous; float64 negative log of selected ratio is the calibrated scalar',
                no_support='null scalar; retained in null query denominator',
                reference_sha256=db.provenance(args.vdjdb)['input_sha256'], reference_population='curated-benchmark',
                species=species, locus=args.locus, exclude_exact=args.exclude_exact, distance_geometry=args.distance_geometry,
                background_scope=args.background_scope, kernel='stem4 v-bounded, trim3/3',
                exact_policy='full anchored junction on one chain, both anchored junctions on the linked axis; fixed model',
                model_sha256=sha256(Path(str(files('vdjmatch.resources')/'vloops'/f'{species}_v_loops.tsv'))),
                source_sha256={m.__name__:sha256(Path(m.__file__)) for m in modules} | {'calibration':sha256(Path(__file__))})


def _withheld(reference, locus, minimum):
    # Match the benchmark query unit of reference-calls, before any hit exclusion.
    if locus=='paired':
        counts=reference.filter(pl.col('gene').is_in(['TRA','TRB'])).group_by('gene','epitope').len()
        raw=build_paired_ref(reference).select(
            'epitope',pl.col('alpha').alias('cdr3a'),pl.col('beta').alias('cdr3b'),
            pl.col('v_alpha').alias('va'),pl.col('v_beta').alias('vb'),'j_alpha','j_beta')
    else:
        c='a' if locus=='TRA' else 'b'
        ref=reference.filter(pl.col('gene')==locus)
        counts=ref.group_by('epitope').len()
        raw=ref.select('epitope',pl.col('cdr3').alias('cdr3'+c),pl.col('v').alias('v'+c),'j')
    eligible=counts.filter(pl.col('len')>=minimum).select('epitope').unique()
    return raw.join(eligible,on='epitope',how='semi').unique().sort(raw.columns).with_row_index('query_id')


def run_database(args, queries, reference, species, prefix, ingestion):
    binding=_binding(args,species)
    options=dict(locus=args.locus,species=species,threads=args.threads,exclude_exact=args.exclude_exact,
                 distance_geometry=args.distance_geometry,background_scope=args.background_scope,
                 kernel_score_budget_bytes=args.kernel_score_budget_bytes)
    prefix.parent.mkdir(parents=True,exist_ok=True)
    meta=dict(command='equal-tail', binding=binding, kernel_score_budget_bytes=args.kernel_score_budget_bytes,
              input_sha256=sha256(Path(args.sample)), ingestion=ingestion)
    recipe=None
    if args.apply_any_call:
        recipe=json.loads(Path(args.apply_any_call).read_text())
        if not compatible_binding(recipe.get('equal_tail_binding'),binding):
            raise ValueError('any-call recipe does not match this equal-tail search')
        source=null_path(args.apply_any_call)
        if sha256(source)!=recipe['null_cases_sha256']:
            raise ValueError('null-case file differs from the fitted recipe')
        nulls=pl.read_csv(source,separator='\t',schema_overrides={'score':pl.Float64})
    if args.fit_any_call is not None:
        withheld=_withheld(reference,args.locus,args.min_null_target_records)
        if not len(withheld):
            raise ValueError('no qualifying peptide-withheld reference queries')
        # Score both populations together once. The withheld label applies only
        # to the corresponding reference query, never to a generated control.
        parts=[]
        for source, q in [('withheld',withheld.rename({'epitope':'_withheld_peptide'})),('generated',queries)]:
            parts.append(q.with_columns(pl.concat_str(pl.lit(source),pl.col('query_id'),separator=':').alias('query_id'),
                                        pl.lit(source).alias('null_source')))
        null_queries=pl.concat(parts,how='diagonal_relaxed')
        # Null generated inputs never supply a withheld peptide. A null label is the sentinel.
        selected=equal_tail_selection(null_queries,reference,**options,withheld_column='_withheld_peptide')
        nulls=selected.select('query_id','score').join(null_queries.select('query_id','null_source'),on='query_id',validate='1:1').with_columns(pl.lit('null').alias('outcome'))
        point,curve=operating_points(nulls,args.fit_any_call)
        path=Path(str(prefix)+'.calibration.json')
        nulls.write_csv(null_path(path),separator='\t')
        recipe=dict(operating_point=point,equal_tail_binding=binding,null_cases_sha256=sha256(null_path(path)),score_column='score',
                      generated_input_sha256=meta['input_sha256'],min_withheld_target_records=args.min_null_target_records,
                      withheld='exclude the annotated peptide from selection at every restriction and axis; fixed reference models',
                      generated='caller-supplied pinned generated control receptors; paired controls require independent chains')
        curve.write_csv(str(prefix)+'.operating_points.tsv',separator='\t')
        path.write_text(json.dumps(recipe,indent=2,allow_nan=False)+'\n')
        meta.update(null_query_counts=nulls.group_by('null_source').len().to_dicts())
    else:
        selected=equal_tail_selection(queries,reference,**options)
        if recipe is not None:
            selected=apply_calibration(selected,recipe['operating_point'])
            selected,tails=null_tail_probabilities(selected,nulls)
            tails.write_csv(str(prefix)+'.null-tails.tsv',separator='\t')
            meta.update(calibration_sha256=sha256(Path(args.apply_any_call)),null_cases_sha256=recipe['null_cases_sha256'])
        selected.write_csv(str(prefix)+'.selection.tsv',separator='\t')
        meta.update(queries=len(selected),status_counts=selected.group_by('selection_status').len().to_dicts())
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(meta,indent=2)+'\n')
    return 0
