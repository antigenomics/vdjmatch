"""masked-motif V-conditioned MF×IKF evidence; not a calibrated likelihood or P-value."""
import polars as pl

from ..aggregate.candidates import PMHC
from ..evalue.epitope_rules import _tokens

_STREAM_TOKEN_ROWS = 1_048_576
_STREAM_OUTPUT_ROWS = 65536


def masked_motif_scores(queries, reference, *, exclude_exact=False, target_reference=None, _consume=None):
    """Return supported V+masked-4-mer scores per query/full pMHC.

    Inputs use full junction ``cdr3`` and ``v``; queries have unique query_id,
    reference has PMHC columns. Call separately for each species and locus.
    Follow masked_motif: IMGT classes, all four masks, distinct junction presence,
    x_length,V,target>=2, x_V,target>=3, x_target>=5, and sum(max(MF_V×IKF))
    per motif divided by query length. All model support gates and global motif
    frequencies remain fixed. Benchmark exclusion masks exact-junction target/V
    contributions only in the hit numerator, without refitting the model.
    Missing/unsupported motifs emit no row, never a prior-only specificity.
    Optional target_reference defines compatible-restriction unions labelled with
    each query's declared PMHC; global motif frequencies still use reference.
    Reference records lacking a nonblank V are ineligible for this V-conditioned
    model; exclude them once before building global and target counts. An empty
    eligible population emits no scores. Queries require complete V, and invalid
    junctions are rejected even on V-ineligible reference records.
    Optional internal ``_consume`` receives bounded supported-score frames and
    returns None. It prepares the same fixed model once, then bounds query token
    and expanded-hit blocks independently; no query-by-target grid is allocated.
    """
    if queries['query_id'].null_count() or queries['query_id'].n_unique()!=queries.height:
        raise ValueError('unique non-null query IDs required')
    for frame in (queries,reference,*([target_reference] if target_reference is not None else [])):
        if frame.filter(~pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)).height:
            raise ValueError('valid complete junction required')
    present_v=pl.col('v').cast(pl.String).str.strip_chars().ne('').fill_null(False)
    if queries.filter(~present_v).height:
        raise ValueError('complete junction and V required')
    reference=reference.with_columns(pl.col('v').cast(pl.String)).filter(present_v)
    if target_reference is not None:
        target_reference=target_reference.with_columns(pl.col('v').cast(pl.String)).filter(present_v)
    ref=reference.select(*PMHC,'cdr3','v').with_columns(pl.col('v').str.replace(r'\*.*$',''))
    ref=ref.unique(maintain_order=True).rename({'cdr3':'seq0','v':'v0'}).with_row_index('_id')
    ref=ref.with_columns(pl.col('seq0').str.len_chars().alias('_length'))
    q=queries.select('query_id','cdr3','v',*(PMHC if target_reference is not None else [])).rename({'cdr3':'seq0','v':'v0'}).with_row_index('_id')
    q=q.with_columns(pl.col('v0').str.replace(r'\*.*$',''))
    def tokens(frame):
        return _tokens(frame,[0],'imgt_wild4',full_junction=True,wildcard_positions=(0,1,2,3))
    rt=tokens(ref).join(ref.drop('v0'),on='_id',validate='m:1')
    target_rt=rt
    if target_reference is not None:
        target=target_reference.select(*PMHC,'cdr3','v').with_columns(pl.col('v').str.replace(r'\*.*$',''))
        target=target.unique(maintain_order=True).rename({'cdr3':'seq0','v':'v0'}).with_row_index('_id')
        target=target.with_columns(pl.col('seq0').str.len_chars().alias('_length'))
        target_rt=tokens(target).join(target.drop('v0'),on='_id',validate='m:1')
    strata=[*PMHC,'v0','_length','t0']
    model=target_rt.group_by(strata).len().rename({'len':'x_lva'}).with_columns(
        pl.col('x_lva').sum().over([*PMHC,'v0','t0']).alias('x_va'),
        pl.col('x_lva').sum().over([*PMHC,'t0']).alias('x_a'),
        pl.lit(ref.height).alias('X')).join(rt.group_by('t0').len().rename({'len':'x'}),on='t0',validate='m:1')
    # Model support gates are fixed before matching any query.
    model=model.filter((pl.col('x_lva')>=2)&(pl.col('x_va')>=3)&(pl.col('x_a')>=5))
    keys=['v0','t0',*(PMHC if target_reference is not None else [])]
    removed=None
    if exclude_exact:
        removed=target_rt.group_by('seq0',*PMHC,'v0','t0').len().rename({'len':'_removed'})

    def reduce(query_tokens):
        hits=query_tokens.join(model,on=keys,how='inner').with_columns(pl.col('x_va').alias('hit_x_va'))
        if removed is not None:
            hits=hits.join(removed,on=['seq0',*PMHC,'v0','t0'],how='left',validate='m:1',nulls_equal=True).with_columns(
                (pl.col('hit_x_va').cast(pl.Int64)-pl.col('_removed').fill_null(0)).alias('hit_x_va')).drop('_removed')
        hits=hits.filter(pl.col('hit_x_va')>0).with_columns(
            (pl.col('hit_x_va')/pl.col('x')*(pl.col('X')/pl.col('x')).log()).alias('_mfikf'))
        return (hits.group_by('query_id',*PMHC,'t0').agg(pl.col('_mfikf').max(),pl.col('seq0').first())
            .group_by('query_id',*PMHC).agg((pl.col('_mfikf').sum()/pl.col('seq0').first().str.len_chars()).alias('motif_score'),
                pl.len().alias('supported_motifs')))

    if _consume is None:
        return reduce(tokens(q).join(q.drop('v0'),on='_id',validate='m:1'))
    if not callable(_consume):
        raise ValueError('motif batch consumer must be callable')
    fanout=model.group_by(keys).len().rename({'len':'_fanout'})
    longest=q['seq0'].str.len_chars().max() or 0
    token_bound=4*max(1,longest-3)
    if token_bound>_STREAM_TOKEN_ROWS:
        raise ValueError('one motif query exceeds the streamed token budget')
    for part in q.iter_slices(max(1,_STREAM_TOKEN_ROWS//token_bound)):
        qt=tokens(part).join(part.drop('v0'),on='_id',validate='m:1')
        expanded=qt.join(fanout,on=keys,how='left').group_by('_id').agg(pl.col('_fanout').sum())
        maximum=expanded['_fanout'].max() or 1
        if maximum>_STREAM_TOKEN_ROWS:
            raise ValueError('one motif query exceeds the streamed expanded-hit budget')
        for block in part.iter_slices(max(1,_STREAM_TOKEN_ROWS//maximum)):
            result=reduce(qt.filter(pl.col('_id').is_between(block['_id'].min(),block['_id'].max())))
            for output in result.iter_slices(_STREAM_OUTPUT_ROWS):
                _consume(output)
    return None


def masked_motif_presence_counts(frame, *, mode='imgt_wild4'):
    """Bounded counts of full/reduced 3-mer or masked-4-mer presence by v0.

    ``_id,seq0,v0`` identify unique input keys, full junctions and nonnegative
    contiguous stratum IDs. Repeated windows count once per input key. Return
    an int64 stratum/token histogram and its fixed vocabulary. ``mode`` uses
    the existing tokenizer's exact3/imgt3/wild4/imgt_wild4 definitions.
    Uses the existing tokenizer in blocks; no per-sequence Python/search calls.
    """
    from itertools import product
    import numpy as np
    from ..evalue.epitope_rules import IMGT_CLASSES, MOTIFS

    if mode not in MOTIFS:
        raise ValueError('unknown motif mode')
    if (not frame.schema['v0'].is_integer() or frame['_id'].null_count() or
            frame['_id'].n_unique() != frame.height or frame['v0'].null_count() or
            frame.filter((pl.col('v0') < 0) | ~pl.col('seq0').str.contains(
                r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False)).height):
        raise ValueError('motif counts require unique IDs, valid full junctions and nonnegative strata')
    letters = list(dict.fromkeys(IMGT_CLASSES.values())) if 'imgt' in mode else list('ACDEFGHIKLMNPQRSTVWY')
    k = 4 if 'wild' in mode else 3
    tokens = []
    for mask in (range(4) if k == 4 else [None]):
        for parts in product(letters, repeat=3):
            word = ''.join(parts)
            tokens.append(word[:mask] + 'X' + word[mask:] if mask is not None else word)
    vocabulary = pl.DataFrame({'t0': tokens}).with_row_index('_token_id')
    groups = int(frame['v0'].max()) + 1 if frame.height else 0
    if groups * len(tokens) * 8 > 1024**3:
        raise ValueError('motif histogram exceeds the 1 GiB count-buffer budget')
    counts = np.zeros((groups, len(tokens)), dtype=np.int64)
    longest = frame['seq0'].str.len_chars().max() or 0
    # Measured expanded-token memory is below 256 B/row; budget a full 256 MiB
    # block, including its tokenizer joins/unique/group reduction temporaries.
    rows = min(65536, max(1, (256 * 1024**2) // ((4 if k == 4 else 1) * max(1, longest - k + 1) * 256)))
    for part in frame.iter_slices(rows):
        frequencies = _tokens(part, [0], mode, full_junction=True,
            wildcard_positions=(0, 1, 2, 3)).group_by('v0', 't0').len().join(
                vocabulary, on='t0', validate='m:1')
        counts[frequencies['v0'].to_numpy(), frequencies['_token_id'].to_numpy()] += frequencies['len'].to_numpy()
    return counts, vocabulary
