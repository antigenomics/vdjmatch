"""Rebuild assay examples from raw AIRR tables through the VDJmatch CLI.

The LLW-versus-LLL sorted-population contrast is kept separate from the mirror's
Pgen-matched repertoire-control tasks. No classifier or threshold is fitted.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys


def prepare(benchmark, sorted_input=None):
    import polars as pl

    parts = []
    for number in (1, 2):
        raw = pl.read_csv(benchmark / f"vdjmatch/airr/dataset_{number}_TRB.tsv.gz", separator="\t")
        parts.append(raw.with_columns(pl.lit('TRB').alias('locus'),pl.lit(number).alias("dataset"),pl.lit('benchmark-mirror').alias('input_origin'),
            pl.col("epitope").alias("source_epitope"),
            pl.col("binder").alias("positive"),
            (pl.col("epitope") + (" sorted positive / negative" if number == 1 else
                " / matched repertoire controls")).alias("task")))
    if sorted_input is None:
        sorted_cells = parts[1].filter((pl.col("binder") == 1) &
            pl.col("source_epitope").is_in(["LLWNGPMAV", "LLLGIGILV"]))
        task='LLW / LLL sorted populations (mirror V calls)'
    else:
        from vdjmatch.io.airr import normalize_annotation_input
        sorted_cells,_=normalize_annotation_input(sorted_input)
        # Keep the legacy amino-acid field out of AIRR's nucleotide cdr3 namespace.
        if 'cdr3' in sorted_cells.columns:
            sorted_cells=sorted_cells.rename({'cdr3': 'source_cdr3'})
        label=next((c for c in ('source_epitope','epitope','antigen.epitope') if c in sorted_cells.columns),None)
        if label is None:
            raise ValueError('sorted input requires the original source population label')
        sorted_cells=sorted_cells.with_columns(pl.col(label).alias('source_epitope'))
        if sorted_cells['source_epitope'].null_count() or not sorted_cells['source_epitope'].is_in(['LLWNGPMAV','LLLGIGILV']).all():
            raise ValueError('sorted input must contain only the declared LLW and LLL populations; no rows are silently removed')
        sorted_cells=sorted_cells.rename({'query_id':'source_observation_id'}).with_columns(
            pl.lit('original-sorted-source').alias('input_origin'),pl.lit(2).alias('dataset'))
        defaults=dict(species='human',locus='TRB',mhc_a='HLA-A*02',mhc_b='B2M',mhc_class='MHCI',binder=1)
        sorted_cells=sorted_cells.with_columns(pl.lit(value).alias(c) for c,value in defaults.items() if c not in sorted_cells.columns)
        task='LLW / LLL sorted populations (original V calls)'
    contrast = sorted_cells.with_columns(
        pl.lit(task).alias("task"),
        (pl.col("source_epitope") == "LLWNGPMAV").cast(pl.Int64).alias("positive"),
        pl.lit("LLWNGPMAV").alias("epitope"))
    return pl.concat([*parts, contrast],how='diagonal_relaxed').with_row_index("query_id")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", type=Path, required=True, help="raw airr_benchmark mirror")
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument('--germline-background',type=Path,required=True,help='raw V/J/junction prior population for retained original scoring')
    parser.add_argument('--sorted-input',type=Path,help='original LLW/LLL sorted table; preserves its V/J calls instead of substituting the mirror')
    parser.add_argument("--output", type=Path, required=True, help="new output directory")
    parser.add_argument("--threads", type=int, default=4, help="CPU worker budget")
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error("--threads must be positive")
    # CLI owns its own worker budget; post-analysis uses one numerical helper pool.
    for key in ("POLARS_MAX_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[key] = "1"
    os.environ["MPLBACKEND"] = "Agg"
    import numpy as np
    import polars as pl
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, precision_recall_curve

    raw = prepare(args.benchmark,args.sorted_input)
    args.output.mkdir(parents=True, exist_ok=False)
    raw.write_csv(args.output / "queries.tsv", separator="\t")
    base = [sys.executable, "-m", "vdjmatch.cli", "evidence", str(args.output / "queries.tsv"),
        "--vdjdb", str(args.reference), "--locus", "TRB", "--targets-from-sample",
        "--reference-population", "curated-benchmark", "--exclude-exact",
        "--threads", str(args.threads),'--junction-ends','trim3','--radius-factor','22.5',
        '--position-weighting','uniform','--neighbour-weighting','kernel','--kernel-scorer','v-bounded','--background','vdjdb-other']
    commands={name:base+['--stem-depth',str(depth),'--distance-matrix',matrix,'--output-prefix',str(args.output/name)]
        for name,depth,matrix in [('gapblock',0,'clipped'),('stem4',4,'blosum62')]}
    commands['original']=[sys.executable,'-m','vdjmatch.cli','historical-density',str(args.output/'queries.tsv'),
        '--vdjdb',str(args.reference),'--locus','TRB','--targets-from-sample','--reference-population','curated-benchmark',
        '--germline-background',str(args.germline_background),'--likelihood-score','--threads',str(args.threads),
        '--output-prefix',str(args.output/'original')]
    import json
    from vdjmatch.db.cache import sha256
    recipes=[]
    from profile_cli import profile
    for name,command in commands.items():
        resources=profile(command,args.output/f'{name}.log',timeout=300)
        recipes.append(dict(method=name,command=command,**resources))
    pl.DataFrame([{k:v for k,v in recipe.items() if k!='command'} for recipe in recipes]).write_csv(
        args.output/'resources.tsv',separator='\t')
    keys = ["query_id", "epitope", "mhc_a", "mhc_b", "mhc_class"]
    scored=raw
    for name in commands:
        suffix='scores' if name=='original' else 'enrichment'
        column='likelihood_score' if name=='original' else 'kernel_score'
        evidence=pl.read_csv(args.output/f'{name}.{suffix}.tsv',separator='\t')
        if evidence.height!=raw.height or raw.select(keys).join(evidence.select(keys),on=keys,
                how='anti',nulls_equal=True).height:
            raise ValueError(f'{name} did not return every source query/pMHC key')
        scored=scored.join(evidence.select(*keys,pl.col(column).alias(name+'_score')),on=keys,
            how='left',validate='1:1',nulls_equal=True)
    assert scored.height==raw.height
    scored.write_csv(args.output / "observations.tsv", separator="\t")
    summaries, curves, counts, cases = [], {}, [], []
    for (task,), part in scored.group_by("task", maintain_order=True):
        y = part["positive"].to_numpy()
        positives, controls = int(y.sum()), int((y == 0).sum())
        for name in commands:
            score=part[name+'_score'].to_numpy();finite=np.isfinite(score)
            score=np.where(finite,score,np.nextafter(score[finite].min(),-np.inf) if finite.any() else -1.)
            fpr,tpr,thresholds=roc_curve(y,score,drop_intermediate=False)
            precision,recall,_=precision_recall_curve(y,score)
            summaries.append(dict(task=task,method=name,observations=len(y),positives=positives,controls=controls,
                missing_scores=int((~finite).sum()),ROC_AUC=roc_auc_score(y,score),average_precision=average_precision_score(y,score)))
            loss=np.empty(len(y))
            for label in (0,1):
                mask=y==label;other=np.sort(score[~mask])
                left=np.searchsorted(other,score[mask],side='left')
                right=np.searchsorted(other,score[mask],side='right')
                loss[mask]=((len(other)-right if label else left)+(right-left)*.5)/len(other)
            assert np.isclose(loss.mean(),1-roc_auc_score(y,score))
            cases.append(part.select('query_id','task','positive','source_epitope','junction_aa','v_call','j_call')
                .with_columns(pl.lit(name).alias('method'),pl.Series('ranking_score',score),
                    pl.Series('opposite_class_rank_loss',loss)))
            curves[task,name]=(fpr,tpr,recall,precision)
            for threshold,tp,fp in zip(thresholds,np.rint(tpr*positives),np.rint(fpr*controls)):
                counts.append(dict(task=task,method=name,threshold=float(threshold),TP=int(tp),FP=int(fp),FN=positives-int(tp),TN=controls-int(fp)))
    summary = pl.DataFrame(summaries)
    summary.write_csv(args.output / "accuracy.tsv", separator="\t")
    pl.DataFrame(counts).write_csv(args.output / "operating-counts.tsv", separator="\t")
    pl.concat(cases).write_csv(args.output/'rank-cases.tsv',separator='\t')
    sources=[args.benchmark/f'vdjmatch/airr/dataset_{i}_TRB.tsv.gz' for i in (1,2)]+[args.reference,args.germline_background]+([args.sorted_input] if args.sorted_input else [])
    (args.output/'manifest.json').write_text(json.dumps(dict(sources={str(p):sha256(p) for p in sources},recipes=recipes,
        observations=raw.height,sorted_origin='original' if args.sorted_input else 'mirror',
        labels='original binder/source population labels retained; positive defines each explicit contrast',
        scope='retained original log evidence, unstemmed gap-block and stem4; all higher-is-better within target',
        calibration='ranking comparison only; no mixed-score or selected-call P-value'),indent=2)+'\n')
    plt.rcParams.update({"font.size":8, "axes.spines.top":False, "axes.spines.right":False,
        "pdf.fonttype":42, "svg.fonttype":"none"})
    tasks=['NLVPMVATV sorted positive / negative',next(task for task,name in curves if task.startswith('LLW / LLL'))]
    fig,axes=plt.subplots(2,2,figsize=(7.09,5.6),layout='constrained')
    for number,task in enumerate(tasks):
        for name,label,colour in [('original','Original','#0072B2'),('gapblock','Gap-block','#D55E00'),('stem4','Stem4','#009E73')]:
            fpr,tpr,recall,precision=curves[task,name]
            row=summary.filter((pl.col('task')==task)&(pl.col('method')==name)).row(0,named=True)
            axes[number,0].plot(fpr,tpr,color=colour,lw=1.3,label=f'{label}; AUC {row["ROC_AUC"]:.3f}')
            axes[number,1].step(recall,precision,where='post',color=colour,lw=1.3,label=f'{label}; AP {row["average_precision"]:.3f}')
        axes[number,0].plot([0,1],[0,1],color='0.6',ls=':',lw=.8)
        axes[number,1].axhline(row['positives']/row['observations'],color='0.6',ls=':',lw=.8)
        label='NLV sorted positive / negative' if number==0 else 'LLW / LLL sorted populations'
        for column,ax in enumerate(axes[number]):
            ax.set(xlim=(0,1),ylim=(0,1.02),xlabel='False-positive rate' if column==0 else 'Recall',
                ylabel='True-positive rate' if column==0 else 'Precision',title=f'{label}\n{row["positives"]} positive / {row["controls"]} comparison rows')
            ax.legend(frameon=False,fontsize=7)
            ax.text(-.15,1.07,'abcd'[number*2+column],transform=ax.transAxes,fontweight='bold',fontsize=10)
    for extension in ("png", "svg", "pdf"):
        fig.savefig(args.output / f"sorted-populations.{extension}", dpi=300, bbox_inches="tight")
    print(summary.write_csv(separator="\t"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
