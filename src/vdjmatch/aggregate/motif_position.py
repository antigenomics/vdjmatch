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


def generative_position_contrast(positions, motifs, *, species, locus, threads=1):
    """Append exact OLGA conditional probabilities using one native masked-DP batch.

    Sum joint Pgen over every model allele pair denoted by the selected genes,
    then divide by the all-wildcard V/J/length mass. No sequence sampling,
    fitting, pseudocounts or empirical-control information enters this model.
    Motif probabilities apply to each start separately; the overlapping
    anywhere union is deliberately unsupported, rather than summing starts.
    """
    from importlib import import_module
    from pathlib import Path
    import numpy as np
    from vdjtools.model.bundled import load_bundled
    from vdjtools.model.native import pack, pgen_aa_degenerate_batch
    from ..db.cache import sha256
    from ..evalue.control import _organism
    from ..evalue.mock_reference import _pgen_model_identity

    if threads < 1:
        raise ValueError('positive native thread budget required')
    organism = _organism(species)
    model = load_bundled(locus, 'olga', organism=organism, collapse=False)
    alleles_for = import_module('vdjtools.overlap.alice')._alleles_for
    _, vi, ji = pack(model)
    v_gene, j_gene, length = positions.select('v_gene', 'j_gene', 'junction_length').unique().row(0)
    vs, js = alleles_for(vi, v_gene), alleles_for(ji, j_gene)
    patterns = [['X']*length]
    for position, aa in positions.select('position', 'amino_acid').iter_rows():
        allowed = ['X']*length
        allowed[position-1] = aa
        patterns.append(allowed)
    starts = motifs.filter(pl.col('position') > 0) if motifs is not None else None
    if starts is not None:
        for position, motif in starts.select('position', 'motif').iter_rows():
            allowed = ['X']*length
            allowed[position-1:position-1+len(motif)] = list(motif)
            patterns.append(allowed)
    queries = [(pattern, va, ja) for va in vs for ja in js for pattern in patterns]
    masses = np.asarray(pgen_aa_degenerate_batch(model, [q[0] for q in queries],
        v=[q[1] for q in queries], j=[q[2] for q in queries], threads=threads), dtype=float
        ).reshape(len(vs)*len(js), len(patterns)).sum(axis=0)
    if not np.isfinite(masses).all() or (masses < 0).any():
        raise ValueError('generative model returned invalid probability masses')
    denominator = float(masses[0])

    def append(frame, values):
        return frame.with_columns(pl.Series('generative_pgen', values, dtype=pl.Float64),
            pl.lit(denominator).alias('generative_stratum_pgen')
        ).with_columns(pl.when(pl.col('generative_pgen').is_null()).then(pl.lit('unsupported_anywhere_union'))
            .when(pl.col('generative_stratum_pgen') == 0).then(pl.lit('unsupported_generative_stratum'))
            .otherwise(pl.lit('available')).alias('generative_status'),
            (pl.col('generative_pgen')/denominator if denominator else pl.lit(None, dtype=pl.Float64)
            ).alias('generative_fraction')
        ).with_columns((pl.col('reference_fraction')-pl.col('generative_fraction')).alias('generative_difference'))

    positions = append(positions, masses[1:1+positions.height])
    if motifs is not None:
        motifs = append(motifs, [None, *masses[1+positions.height:]])
    identity = _pgen_model_identity(organism, locus)
    identity.update(collapse_alleles=False, v_alleles=vs, j_alleles=js,
        native_batches=1, native_queries=len(queries), distinct_patterns=len(patterns),
        stratum_pgen=denominator, threads=threads,
        source_sha256={fn.__name__: sha256(Path(fn.__code__.co_filename))
            for fn in (load_bundled, pack, pgen_aa_degenerate_batch, alleles_for, generative_position_contrast)})
    return positions, motifs, identity
