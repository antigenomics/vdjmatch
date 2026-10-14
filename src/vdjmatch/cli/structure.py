"""Reference-only graph, junction-length and V-usage composition."""
from pathlib import Path


def register(subparsers):
    p = subparsers.add_parser('reference-structure', help=__doc__)
    p.add_argument('--vdjdb', required=True)
    p.add_argument('--species', default='HomoSapiens')
    p.add_argument('--locus', choices=['TRA','TRB'], required=True)
    p.add_argument('--epitope', action='append')
    p.add_argument('--targets-from-sample', help='compatible pMHC unions declared in an AIRR table')
    p.add_argument('--reference-population', choices=['curated-benchmark','all-valid'], default='curated-benchmark')
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--output-prefix', required=True)
    p.set_defaults(func=main)


def composition(records, junctions):
    """Full and exact-junction-punctured length/V summaries; no query labels."""
    import math
    import numpy as np
    import polars as pl

    seqs = junctions['cdr3'].to_list()
    lengths = np.array([len(s) for s in seqs], dtype=np.float64)
    keys = records.select('cdr3','_allele').unique().with_columns(
        pl.col('_allele').str.split('*').list.first().alias('_gene'))
    counts = keys.group_by('_gene').len().sort('_gene')
    sizes = counts['len'].to_numpy().astype(float)
    total = keys.height
    flog = lambda x: x*math.log(x) if x else 0.
    v_log = math.fsum(flog(x) for x in sizes)
    length_counts = dict(zip(*np.unique(lengths, return_counts=True)))
    length_log = math.fsum(flog(x) for x in length_counts.values())
    sum_length = lengths.sum(); sum_square = (lengths*lengths).sum()
    records_per_seq = dict(records.group_by('cdr3').len().iter_rows())
    removed_keys = dict(keys.group_by('cdr3').len().iter_rows())
    # Each sequence removes its full contribution from every V category.
    changes = keys.group_by('cdr3','_gene').len().rename({'len':'_removed'}).join(
        counts,on='_gene',validate='m:1').sort('cdr3','_gene')
    v_changes = {}
    for sequence, part in changes.group_by('cdr3',maintain_order=True):
        before = part['len'].to_list(); after = (part['len']-part['_removed']).to_list()
        v_changes[sequence[0]] = (math.fsum(flog(new)-flog(old) for old,new in zip(before,after)),
                                  sum(x==0 for x in after))

    def row(removed=None):
        n = len(seqs)-(removed is not None)
        length = len(removed) if removed is not None else 0
        k = removed_keys[removed] if removed is not None else 0
        remain = total-k
        delta, lost_genes = v_changes[removed] if removed is not None else (0.,0)
        v_entropy = max(0., math.log(remain)-(v_log+delta)/remain) if remain else None
        f_length = length_log
        if removed is not None:
            count = length_counts[float(length)]
            f_length += flog(count-1)-flog(count)
        mean = (sum_length-length)/n if n else None
        variance = max(0.,(sum_square-length*length)/n-mean*mean) if n else None
        return dict(clonotype_keys=remain, v_genes=len(sizes)-lost_genes,
            v_entropy=v_entropy, effective_v_genes=math.exp(v_entropy) if v_entropy is not None else None,
            length_mean=mean, length_sd=math.sqrt(variance) if variance is not None else None,
            length_entropy=max(0.,math.log(n)-f_length/n) if n else None,
            reference_records=records.height-(records_per_seq[removed] if removed is not None else 0),
            excluded_reference_records=records_per_seq[removed] if removed is not None else 0,
            excluded_clonotype_keys=k)
    return row(), [row(sequence) for sequence in seqs]


def main(a):
    import json
    import time
    import polars as pl
    from importlib.metadata import version
    import seqtree._core as native_core
    from .. import db
    from ..aggregate.candidates import PMHC
    from ..aggregate.structure import graph_profiles
    from ..cluster import overlap
    from ..db.cache import sha256
    from ..evalue.control import _organism
    from ..match.geometry import load_v_loops, resolve_v_alleles
    from ..runner.annotation import _validate_output_paths

    start = time.monotonic()
    if a.threads < 1:
        raise ValueError('positive threads required')
    protected = [a.vdjdb]+([a.targets_from_sample] if a.targets_from_sample else [])
    prefix = _validate_output_paths(a.output_prefix, protected)
    organism = _organism(a.species)
    species = {'human':'HomoSapiens','mouse':'MusMusculus'}[organism]
    raw = db.load(a.vdjdb, species=species, gene=a.locus, epitope=a.epitope)
    if a.reference_population == 'all-valid':
        raw = raw.with_columns(pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False).alias('reference_valid'))
    model = load_v_loops(organism)
    raw = raw.with_columns(resolve_v_alleles(raw['v'],model,locus=a.locus).alias('_allele'))
    available = pl.col('reference_valid') & pl.col('_allele').is_not_null() & (pl.col('cdr3').str.len_chars() >= 8)
    if a.targets_from_sample:
        from ..io.airr import _read_table
        from .historical import assayed_targets
        targets = _read_table(a.targets_from_sample)
        if a.epitope:
            targets = targets.filter(pl.col('epitope').is_in(a.epitope))
        if 'query_id' not in targets.columns:
            targets = targets.with_row_index('query_id')
        groups = [(task, selected) for task, _, selected in
            assayed_targets(targets.select('query_id'),targets,raw,a.species)]
    else:
        groups = [(dict(zip(PMHC,key)), part) for key,part in raw.group_by(PMHC,maintain_order=True)]
    # One batched native search for the whole reference population, then subset its
    # edges by declared target membership; no per-epitope native search loop.
    selected_nodes = [source.filter(available).select('cdr3') for _,source in groups]
    nodes = (pl.concat(selected_nodes).unique(maintain_order=True) if selected_nodes else raw.head(0).select('cdr3')).with_row_index('_node')
    pairs = overlap(nodes['cdr3'].to_list(),scope='1,0,0,1',threads=a.threads) if nodes.height else None
    edges = pairs.filter(pl.col('a_idx')<pl.col('b_idx')).select('a_idx','b_idx') if pairs is not None else None
    memberships = []
    for group_id, (_, source) in enumerate(groups):
        memberships.append(source.filter(available).select('cdr3').unique(maintain_order=True)
            .with_row_index('_local').join(nodes,on='cdr3',validate='1:1')
            .with_columns(pl.lit(group_id).alias('_group')))
    membership = pl.concat(memberships) if memberships else None
    grouped_edges = {}
    if membership is not None and edges is not None and edges.height:
        left = membership.select(pl.col('_node').cast(pl.Int64).alias('a_idx'),'_group',pl.col('_local').alias('a'))
        right = membership.select(pl.col('_node').cast(pl.Int64).alias('b_idx'),'_group',pl.col('_local').alias('b'))
        local = edges.join(left,on='a_idx',how='inner',validate='m:m').join(
            right,on=['b_idx','_group'],how='inner',validate='m:1')
        grouped_edges = {key[0]: frame.select('a','b') for key, frame in local.group_by('_group')}
    summaries, punctures = [], []
    for group_id, (task, source) in enumerate(groups):
        retained = source.filter(available)
        junctions = memberships[group_id]
        local_edges = grouped_edges[group_id].iter_rows() if group_id in grouped_edges else []
        full, removed = graph_profiles(junctions.height,local_edges)
        v_full, v_removed = composition(retained,junctions)
        context = dict(species=species,locus=a.locus,**task,source_observations=source.height,
                       unavailable_records=source.height-retained.height)
        summaries.append(context | full | v_full | dict(structure_status='available' if full['junctions'] else 'empty_reference'))
        for sequence, graph, usage in zip(junctions['cdr3'],removed,v_removed):
            punctures.append(context | dict(excluded_junction=sequence, structure_status='available' if graph['junctions'] else 'empty_after_puncture') | graph | usage)
    prefix.parent.mkdir(parents=True,exist_ok=True)
    summary_path = Path(str(prefix)+'.structures.tsv')
    puncture_path = Path(str(prefix)+'.punctures.tsv')
    if not summaries:
        raise ValueError('no declared reference targets')
    pl.DataFrame(summaries,infer_schema_length=None).write_csv(summary_path,separator='\t')
    if punctures:
        pl.DataFrame(punctures,infer_schema_length=None).write_csv(puncture_path,separator='\t')
    else:
        pl.DataFrame({'excluded_junction':pl.Series([],dtype=pl.String)}).write_csv(puncture_path,separator='\t')
    Path(str(prefix)+'.manifest.json').write_text(json.dumps(dict(reference=db.provenance(a.vdjdb),
        reference_population=a.reference_population, species=species,locus=a.locus,threads=a.threads,
        epitope=a.epitope, targets_sha256=sha256(Path(a.targets_from_sample)) if a.targets_from_sample else None,
        software={name:version(name) for name in ('vdjmatch','seqtree','polars','numpy')},
        native_sha256=sha256(Path(native_core.__file__)),
        v_model_sha256=sha256(Path(load_v_loops.__code__.co_filename).parents[1]/'resources/vloops'/f'{organism}_v_loops.tsv'),
        loader_sha256=sha256(Path(db.load.__code__.co_filename)),
        search_source_sha256=sha256(Path(overlap.__code__.co_filename)),
        native_batches=1 if nodes.height else 0, searched_unique_junctions=nodes.height,
        global_hamming_edges=edges.height if edges is not None else 0, targets=len(summaries),
        puncture_rows=len(punctures), graph_scope='full-junction Hamming<=1; no indels',
        graph_node_unit='unique full junction', v_usage_unit='unique full junction/resolved V allele keys; gene-level categories',
        entropy_units='natural logarithms', length_units='amino-acid residues, conserved anchors included',
        exclusion='puncture removes the complete full junction, all V annotations and incident edges; articulation fragments recomputed exactly',
        interpretation='reference-only morphology; no fitted weights, recognition probability or P/E-value',
        source_sha256=sha256(Path(__file__)), reducer_sha256=sha256(Path(graph_profiles.__code__.co_filename)),
        tables={str(p.name):sha256(p) for p in (summary_path,puncture_path)}, seconds=time.monotonic()-start),indent=2)+'\n')
    return 0
