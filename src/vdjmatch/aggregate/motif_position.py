"""Descriptive positional junction contrasts conditional on both genes and length."""
import polars as pl

from .structure import _gene_expression

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'


def position_contrast(reference, control, *, locus, v_gene, j_gene, junction_length, motif=None):
    """Return residue and optional exact-motif fractions for unique junction/V/J keys.

    Alleles are pooled using the existing gene projection. Positions are 1-based
    full-junction coordinates including conserved anchors; motif position zero
    denotes presence anywhere, counting each key once. Control event weights
    are deliberately ignored. Missing control strata yield null differences.
    """
    if locus not in {'TRA', 'TRB'} or junction_length < 1:
        raise ValueError('TRA/TRB and a positive junction length required')
    genes = pl.DataFrame({'v': [v_gene], 'j': [j_gene]}).select(
        _gene_expression(locus, 'V').alias('v'), _gene_expression(locus, 'J').alias('j'))
    v_gene, j_gene = genes.row(0)
    if v_gene is None or j_gene is None:
        raise ValueError('unambiguous matching V and J genes required')
    if motif is not None and (not motif or len(motif) > junction_length or
                              any(aa not in AMINO_ACIDS for aa in motif)):
        raise ValueError('motif must be canonical amino acids and fit within the junction')

    def selected(frame):
        return frame.select('cdr3', _gene_expression(locus, 'V').alias('v'),
            _gene_expression(locus, 'J').alias('j')).filter(
                (pl.col('v') == v_gene) & (pl.col('j') == j_gene) &
                (pl.col('cdr3').str.len_chars() == junction_length) &
                pl.col('cdr3').str.contains(r'^[ACDEFGHIKLMNPQRSTVWY]+$').fill_null(False))

    observed = selected(reference)
    background = selected(control.rename({'junction_aa': 'cdr3'})).unique()
    records = observed.height
    observed = observed.unique()
    if not observed.height:
        raise ValueError('no reference keys in the requested V/J/length stratum')
    n, m = observed.height, background.height

    def finish(table):
        return table.with_columns(pl.lit(v_gene).alias('v_gene'), pl.lit(j_gene).alias('j_gene'),
            pl.lit(junction_length).alias('junction_length'), pl.lit(records).alias('reference_records'),
            pl.lit(n).alias('reference_keys'), pl.lit(m).alias('control_keys'),
            (pl.col('reference_count') / n).alias('reference_fraction'),
            (pl.col('control_count') / m if m else pl.lit(None, dtype=pl.Float64)).alias('control_fraction'),
            pl.lit('available' if m else 'unsupported_control_stratum').alias('control_status')
        ).with_columns((pl.col('reference_fraction') - pl.col('control_fraction')).alias('difference'))

    def counts(frame):
        return frame.select(pl.col('cdr3').str.split('').alias('amino_acid'),
            pl.int_ranges(1, junction_length + 1).alias('position')).explode(
                'position', 'amino_acid').group_by('position', 'amino_acid').len()

    grid = pl.DataFrame({'position': range(1, junction_length + 1)}).join(
        pl.DataFrame({'amino_acid': list(AMINO_ACIDS)}), how='cross')
    positions = finish(grid.join(counts(observed).rename({'len': 'reference_count'}),
        on=['position', 'amino_acid'], how='left').join(
            counts(background).rename({'len': 'control_count'}), on=['position', 'amino_acid'],
            how='left').with_columns(pl.col('reference_count', 'control_count').fill_null(0).cast(pl.Int64)))
    motifs = None
    if motif is not None:
        expressions = [pl.col('cdr3').str.contains(motif, literal=True).sum().alias('0')]
        expressions += [(pl.col('cdr3').str.slice(i, len(motif)) == motif).sum().alias(str(i+1))
            for i in range(junction_length-len(motif)+1)]
        motifs = finish(pl.DataFrame({'position': range(junction_length-len(motif)+2),
            'motif': motif, 'reference_count': observed.select(expressions).row(0),
            'control_count': background.select(expressions).row(0)}).with_columns(
                pl.col('reference_count', 'control_count').cast(pl.Int64)))
    return positions.sort('position', 'amino_acid'), motifs
