"""Render recorded repertoire enrichment against every mock; no inference is refitted."""
from pathlib import Path
from textwrap import fill
import json


def register(sub):
    p = sub.add_parser('repertoire-report', help=__doc__)
    p.add_argument('counts', type=Path, help='Final real/mock count TSV (replicate -1 is observed)')
    p.add_argument('summary', type=Path, help='Final summary TSV with recorded empirical P and BH q values')
    p.add_argument('--metric', default='curve_mean')
    p.add_argument('--output-prefix', type=Path, required=True)
    p.set_defaults(func=main)


def render(counts, summary, metric, prefix):
    """Plot the recorded null distribution for every sample/target and export its evidence."""
    import numpy as np
    import polars as pl
    import matplotlib.pyplot as plt
    keys = ['sample_id', 'group_id']
    value, q = 'matched_' + metric, 'original_q_' + metric
    pvalue = 'original_p_' + metric
    for frame, required in [(counts, keys + ['replicate', value]), (summary, keys + [pvalue, q])]:
        if not set(required).issubset(frame.columns):
            raise ValueError('missing required report columns: ' + ', '.join(sorted(set(required)-set(frame.columns))))
    if summary.select(keys).is_duplicated().any():
        raise ValueError('summary requires unique sample/target identities')
    if counts.is_empty() or summary.is_empty():
        raise ValueError('nonempty count and summary tables required')
    rows = []
    groups = counts.partition_by(keys, as_dict=True, maintain_order=True)
    prefix = Path(prefix); prefix.parent.mkdir(parents=True, exist_ok=True)
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(str(prefix) + '.pdf') as pdf:
        for start in range(0, len(groups), 12):
            page_size = min(12, len(groups)-start)
            nrows = (page_size+3)//4
            fig, axes = plt.subplots(nrows, 4, squeeze=False, figsize=(12, 2.5*nrows+.6), layout='constrained')
            for ax, (key, part) in zip(axes.flat, list(groups.items())[start:start+12]):
                observed = part.filter(pl.col('replicate') == -1)
                mock = part.filter(pl.col('replicate') >= 0).sort('replicate')
                selected = summary.filter((pl.col('sample_id') == key[0]) & (pl.col('group_id') == key[1]))
                if observed.height != 1 or mock.height != 200 or mock['replicate'].n_unique() != 200 or selected.height != 1:
                    plt.close(fig)
                    raise ValueError('each sample/target requires one observed row, 200 distinct mocks and one summary')
                x = mock[value].to_numpy(); real = observed[value][0]; adjusted = selected[q][0]; nominal = selected[pvalue][0]
                if not np.isfinite(x).all() or not np.isfinite(real):
                    plt.close(fig); raise ValueError('finite endpoint values required')
                ax.hist(x, bins='auto', color='#999999', alpha=.75)
                ax.axvline(real, color='#0072B2', linewidth=2, label='Observed')
                sample_label = selected['sample_label'][0] if 'sample_label' in selected.columns else key[0]
                label = f'{sample_label} | P={nominal:.4g}' if nominal is not None else f'{sample_label} | P unavailable'
                if adjusted is not None: label += f' | q={adjusted:.4g}'
                target = str(key[1])
                if target.startswith('{"'): target = ' | '.join(map(str, json.loads(target).values()))
                ax.set_title(label + '\n' + fill(target, 35), fontsize=8)
                ax.set(xlabel='Mean discoverable clonotypes' if metric == 'curve_mean' else metric.replace('_', ' '), ylabel='Mock panels')
                ax.spines[['top','right']].set_visible(False)
                rows.append(dict(sample_id=key[0], group_id=key[1], observed=real, mock_panels=200,
                                 mock_mean=float(x.mean()), real_mock_ratio=float(real/x.mean()) if x.mean() else None, recorded_p=nominal, recorded_q=adjusted))
            for ax in list(axes.flat)[len(list(groups.items())[start:start+12]):]:
                ax.axis('off')
            fig.suptitle('Observed repertoire endpoint versus 200 matched mock panels', fontsize=13)
            axes.flat[0].legend(frameon=False, fontsize=8)
            pdf.savefig(fig)
            fig.savefig(str(prefix) + f'.page-{start//12+1}.png', dpi=180)
            plt.close(fig)
    evidence = pl.DataFrame(rows)
    evidence.write_csv(str(prefix) + '.evidence.tsv', separator='\t')
    return evidence


def main(args):
    import polars as pl
    import matplotlib
    matplotlib.use('Agg')
    render(pl.read_csv(args.counts, separator='\t'), pl.read_csv(args.summary, separator='\t'), args.metric, args.output_prefix)
