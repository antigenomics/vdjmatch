"""Score all locally available vaccination repertoires through one fresh CLI batch.

Preserve technical replicates, complete TRB read denominators and the declared
A*02 restriction. Rank by LLW reference-neighbour evidence rather than abundance.
Raw clonotypes, filenames and donor mappings stay in the private output directory.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys


def prepare(directory):
    import polars as pl
    from vdjmatch.io.columns import validate_counts

    metadata=pl.read_csv(directory/'metadata.txt',separator='\t',infer_schema_length=0)
    required={'file_name','donor','day','replica'}
    if not required.issubset(metadata.columns) or metadata['file_name'].n_unique()!=metadata.height:
        raise ValueError('metadata requires unique file_name and donor/day/replica fields')
    paths=sorted(directory.glob('*.airr.tsv.gz'))
    if not paths:
        raise ValueError('no local AIRR repertoires')
    unknown={p.name for p in paths}-set(metadata['file_name'])
    if unknown:
        raise ValueError('local AIRR tables missing from metadata')
    selected=metadata.filter(pl.col('file_name').is_in([p.name for p in paths]))
    if selected.select(pl.any_horizontal(pl.col(c).is_null() | (pl.col(c)=='') for c in required).any()).item():
        raise ValueError('selected metadata contains missing identities or time points')
    selected=selected.with_columns(pl.col('day').cast(pl.Int64))
    donors={name:f'Donor {i+1}' for i,name in enumerate(sorted(selected['donor'].unique()))}
    parts,denominators=[] ,[]
    keys=['v_call','junction_aa','j_call']
    for observation in selected.sort('donor','day','replica').iter_rows(named=True):
        path=directory/observation['file_name']
        names=pl.read_csv(path,separator='\t',n_rows=0).columns
        v='v_call' if 'v_call' in names else 'v_gene'
        j='j_call' if 'j_call' in names else 'j_gene'
        raw=pl.read_csv(path,separator='\t',columns=['locus','junction_aa',v,j,'duplicate_count'])
        validate_counts(raw)
        raw=raw.with_columns(pl.col('duplicate_count').cast(pl.Int64),
            pl.col(v).alias('v_call'),pl.col(j).alias('j_call'))
        trb=pl.col('locus')=='TRB'
        valid=pl.col('junction_aa').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)
        canonical=pl.col('junction_aa').str.contains(r'^C[ACDEFGHIKLMNPQRSTVWY]+[FW]$').fill_null(False)
        called=pl.col('v_call').is_not_null() & (pl.col('v_call').str.strip_chars()!='')
        eligible=raw.filter(trb & valid & called)
        clones=eligible.group_by(keys).agg(pl.col('duplicate_count').sum()).sort(keys)
        label=dict(donor=donors[observation['donor']],day=int(observation['day']),replicate=observation['replica'])
        denominators.append(dict(**label,raw_rows=raw.height,eligible_rows=eligible.height,
            non_trb_rows=raw.filter(~trb.fill_null(False)).height,
            invalid_junction_rows=raw.filter(trb & ~valid).height,
            noncanonical_anchor_rows=raw.filter(trb & valid & called & ~canonical).height,
            missing_v_rows=raw.filter(trb & valid & ~called).height,
            eligible_clonotypes=clones.height,eligible_read_counts=int(eligible['duplicate_count'].sum())))
        parts.append(clones.with_columns(pl.lit(value).alias(key) for key,value in label.items()))
    observations=pl.concat(parts)
    queries=observations.select(keys).unique().sort(keys).with_row_index('query_id').with_columns(
        pl.col('junction_aa').str.contains(r'^C[ACDEFGHIKLMNPQRSTVWY]+[FW]$').alias('canonical_anchor_residues'),
        pl.lit('TRB').alias('locus'),pl.lit('human').alias('species'),
        pl.lit('LLWNGPMAV').alias('epitope'),pl.lit('HLA-A*02').alias('mhc_a'),
        pl.lit('B2M').alias('mhc_b'),pl.lit('MHCI').alias('mhc_class'))
    observations=observations.join(queries.select('query_id',*keys),on=keys,how='left',
        validate='m:1',nulls_equal=True)
    assert observations['query_id'].null_count()==0
    return observations,queries,pl.DataFrame(denominators),selected,donors


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset',type=Path,required=True,help='local AIRR vaccination mirror, including metadata.txt')
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True,help='new private output directory')
    parser.add_argument('--threads',type=int,default=4,help='native CPU worker budget')
    parser.add_argument('--memory-gb',type=float,default=24,help='sampled CLI process-tree RSS ceiling in decimal GB, at most 32')
    args=parser.parse_args(argv)
    if args.threads<1 or not 0<args.memory_gb<=32:
        parser.error('threads must be positive and memory-gb must be in (0,32]')
    for key in ('POLARS_MAX_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
        os.environ[key]='1'
    os.environ['MPLBACKEND']='Agg'
    import numpy as np
    import polars as pl
    import matplotlib.pyplot as plt
    from vdjmatch.db.cache import sha256
    from profile_cli import profile

    args.output.mkdir(parents=True,exist_ok=False)
    observations,queries,coverage,metadata,donors=prepare(args.dataset)
    queries.write_csv(args.output/'queries.tsv',separator='\t')
    coverage.write_csv(args.output/'coverage.tsv',separator='\t')
    source_paths=[args.dataset/name for name in metadata['file_name']]+[args.dataset/'metadata.txt',args.reference]
    manifest=dict(sources={str(p):sha256(p) for p in source_paths},donor_mapping=donors,
        local_repertoires=metadata.height,local_donors=len(donors),unique_query_keys=queries.height,
        observation_unit='amino-acid junction/V/J clonotype per technical replicate; reads are not cells',
        restriction='author-declared HLA-A*02 for all donors; compatible reference union',
        sequence_eligibility='producer-declared full junction, standard amino acids, TRB and supplied V; noncanonical anchor residues retained and flagged; functionality unknown',
        calibration='relative target evidence and radius support only; no binder classification')
    command=[sys.executable,'-m','vdjmatch.cli','evidence',str(args.output/'queries.tsv'),
        '--vdjdb',str(args.reference),'--locus','TRB','--targets-from-sample',
        '--selection-only','--reference-population','all-valid','--threads',str(args.threads),
        '--verbose',
        '--output-prefix',str(args.output/'llw')]
    manifest['command']=command
    (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    usage=profile(command,args.output/'llw.log',timeout=1800,max_rss_bytes=int(args.memory_gb*1e9))
    selection=pl.read_csv(args.output/'llw.selection.tsv',separator='\t')
    status=pl.read_csv(args.output/'llw.queries.tsv',separator='\t')
    if status.height!=queries.height or selection.height!=queries.height:
        raise ValueError('CLI did not retain every query key')
    evidence=selection.select('query_id','kernel_score',(pl.col('exact_peptides')>0).alias('exact_supported')).join(
        status.select('query_id','status','available'),on='query_id',validate='1:1')
    evidence=evidence.join(queries.select('query_id','canonical_anchor_residues'),on='query_id',validate='1:1')
    view=observations.join(evidence,on='query_id',how='left',validate='m:1')
    view.write_csv(args.output/'observations.tsv',separator='\t')
    summaries=view.group_by('donor','day','replicate').agg(pl.len().alias('clonotypes'),
        pl.col('duplicate_count').sum().alias('read_counts'),
        pl.col('available').sum().alias('searchable_clonotypes'),
        (pl.col('status')=='matched').sum().alias('radius_supported_clonotypes'),
        pl.col('exact_supported').fill_null(False).sum().alias('exact_supported_clonotypes'),
        pl.col('duplicate_count').filter(pl.col('status')=='matched').sum().alias('radius_supported_read_counts'))
    summaries.sort('donor','day','replicate').write_csv(args.output/'support.tsv',separator='\t')
    pl.DataFrame([dict(workflow='All local LLW repertoires',query_keys=queries.height,
        native_threads=args.threads,**usage)]).write_csv(args.output/'resources.tsv',separator='\t')
    # Complete score-tied groups enter together; curves select by evidence, never abundance.
    curves=[]
    for labels,part in view.group_by('donor','day','replicate'):
        ranked=part.with_columns(pl.col('kernel_score').fill_nan(None)).group_by('kernel_score').agg(
            pl.len().alias('clonotypes'),pl.col('duplicate_count').sum().alias('reads')).sort(
                'kernel_score',descending=True,nulls_last=True)
        ranked=ranked.with_columns(pl.col('clonotypes').cum_sum().alias('cumulative_clonotypes'),
            (pl.col('reads').cum_sum()/part['duplicate_count'].sum()).alias('cumulative_read_fraction'),
            (pl.col('clonotypes').cum_sum()/part.height).alias('cumulative_clonotype_fraction'),
            *(pl.lit(value).alias(key) for key,value in zip(('donor','day','replicate'),labels)))
        curves.append(ranked)
    curves=pl.concat(curves).sort('donor','day','replicate','cumulative_clonotypes')
    curves.write_csv(args.output/'evidence-coverage.tsv',separator='\t')
    plt.rcParams.update({'font.size':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,2,figsize=(7.09,3.0),layout='constrained')
    palette=['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9']
    for i,donor in enumerate(sorted(donors.values())):
        for j,(labels,part) in enumerate(curves.filter(pl.col('donor')==donor).group_by('day','replicate',maintain_order=True)):
            axes[0].step(part['cumulative_clonotype_fraction'],part['cumulative_read_fraction'],
                where='post',color=palette[i%len(palette)],ls=['-','--',':','-.'][j%4],lw=1,label=donor if j==0 else None)
    axes[0].set(xlabel='Clonotype fraction, ranked by LLW evidence',ylabel='Cumulative TRB read-count fraction',xlim=(0,1),ylim=(0,1))
    axes[0].legend(frameon=False,fontsize=7)
    for i,donor in enumerate(sorted(donors.values())):
        part=summaries.filter(pl.col('donor')==donor)
        axes[1].scatter(np.full(part.height,i),part['radius_supported_clonotypes']/part['clonotypes'],
            color=palette[i%len(palette)],s=22)
    axes[1].set(xticks=range(len(donors)),xticklabels=[str(i+1) for i in range(len(donors))],
        xlabel='Donor',ylabel='Fraction with an LLW radius neighbour',ylim=(0,None))
    for ax,letter in zip(axes,'ab'):
        ax.text(-.17,1.04,letter,transform=ax.transAxes,fontweight='bold',fontsize=10)
    for ext in ('png','svg','pdf'):
        fig.savefig(args.output/f'llw-neighbour-evidence.{ext}',dpi=300,bbox_inches='tight')
    print(summaries.sort('donor','day','replicate').write_csv(separator='\t'))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
