"""Epitope-document IDF, distinct from control-junction rarity."""
import polars as pl

from .structure import _gene_expression
from ..evalue.epitope_rules import _tokens, MOTIFS


def epitope_motif_documents(rows, *, locus, mode, min_records=30):
    """Count eligible epitope documents per V/length/token presence stratum.

    Caller supplies valid rows of one species/locus. Eligibility counts source
    records before distinct junction/V keys. Full junctions and all four mask
    positions match the control-relative diagnostic tokenizer. Each epitope
    contributes once per token, regardless of receptor or record multiplicity.
    """
    if locus not in {'TRA', 'TRB'} or mode not in MOTIFS or min_records < 1:
        raise ValueError('valid locus, motif mode and positive min_records required')
    eligible = rows.group_by('epitope').len().filter(pl.col('len') >= min_records).select('epitope')
    keys = rows.join(eligible, on='epitope', how='semi').select('epitope',
        pl.col('cdr3').alias('seq0'), _gene_expression(locus, 'V').alias('v0')).unique().with_columns(
            pl.col('seq0').str.len_chars().alias('junction_length')).with_row_index('_id')
    # Bound expanded tokenizer/join scratch using the existing 256 B/token budget.
    longest = keys['seq0'].str.len_chars().max() or 0
    block = max(1, min(65536, (256*1024**2)//(256*max(1, longest)*(4 if 'wild' in mode else 1))))
    parts = []
    for part in keys.iter_slices(block):
        tokens = _tokens(part, [0], mode, full_junction=True, wildcard_positions=(0, 1, 2, 3))
        parts.append(tokens.join(part.select('_id', 'epitope', 'junction_length'),
            on='_id', validate='m:1').select('epitope', 'v0', 'junction_length', 't0').unique())
    tokens = pl.concat(parts) if parts else pl.DataFrame(schema={
        'epitope': pl.String, 'v0': pl.String, 'junction_length': pl.UInt32, 't0': pl.String})
    output = []
    for scope, strata in [('exact_length', ['v0', 'junction_length']), ('all_lengths', ['v0'])]:
        documents = keys.group_by(strata).agg(pl.col('epitope').n_unique().alias('epitope_documents'))
        frequency = tokens.group_by(*strata, 't0').agg(
            pl.col('epitope').n_unique().alias('epitope_document_frequency'))
        table = frequency.join(documents, on=strata, nulls_equal=True, validate='m:1').with_columns(
            pl.lit(scope).alias('length_scope'))
        if scope == 'all_lengths':
            table = table.with_columns(pl.lit(None, dtype=pl.UInt32).alias('junction_length'))
        output.append(table.rename({'v0': 'v_gene', 't0': 'motif'}))
    return pl.concat(output, how='diagonal_relaxed').with_columns(
        (1+((pl.col('epitope_documents')+1)/(pl.col('epitope_document_frequency')+1)).log()).alias('epitope_idf')
    ).sort('length_scope', 'v_gene', 'junction_length', 'motif', nulls_last=True)
