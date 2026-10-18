"""Selection in gene-consensus aligned germline-loop residue projections."""
import polars as pl

from .candidates import PMHC
from .structure import compare_gene_usage, _selection_information
from ..evalue.epitope_rules import IMGT_CLASSES

_AXES = ['loop', 'alignment_column', 'alphabet']
_CONTEXT = ['species', 'locus', *PMHC, 'segment']


def germline_projection_map(model, *, locus):
    """Project V genes only where every bundled allele agrees.

    Columns are one-based model alignment columns, not structural residue numbers.
    Full and eleven-class alphabets are separate projections. A dot is a gap in
    an otherwise populated loop; an entirely empty loop template is unavailable.
    Omitted alleles never select ``*01``. Caller retains original gene/allele calls.
    """
    if locus not in {'TRA', 'TRB'}:
        raise ValueError('germline projection requires TRA or TRB')
    genes = {}
    for allele, row in model.items():
        if row['locus'] == locus:
            genes.setdefault(allele.split('*')[0], []).append(row)
    rows = []
    for loop in ('cdr1', 'cdr2', 'cdr25'):
        widths = {len(row[loop]) for alleles in genes.values() for row in alleles}
        if len(widths) != 1:
            raise ValueError('germline templates require one aligned width per loop')
        for gene, alleles in sorted(genes.items()):
            missing = any(not (set(row[loop]) & set(IMGT_CLASSES)) for row in alleles)
            for column in range(next(iter(widths))):
                raw = {row[loop][column] for row in alleles}
                if raw - (set(IMGT_CLASSES) | {'.', '*', 'X'}):
                    raise ValueError('germline templates require standard amino acids or aligned gaps')
                for alphabet in ('full', 'imgt'):
                    values = raw if alphabet == 'full' else {IMGT_CLASSES.get(r, r) for r in raw}
                    status = ('missing_loop_template' if missing else 'non_amino_acid_template' if raw & {'*', 'X'}
                        else 'allele_dependent' if len(values) != 1 else 'resolved')
                    rows.append(dict(v_gene=gene, loop=loop, alignment_column=column+1, alphabet=alphabet,
                        residue=next(iter(values)) if status == 'resolved' else None,
                        modeled_residues=''.join(sorted(values)), modeled_alleles=len(alleles), projection_status=status))
    return pl.DataFrame(rows, schema_overrides={'residue': pl.String})


def iter_germline_information(usage, control, mapping):
    """Yield point factors and KL for one loop-column/alphabet projection at a time.

    Consume final V-gene composition counts. Preserve distinct junction/V-gene
    key mass when pooling genes; do not deduplicate keys after residue projection.
    Both pooled and exact-length information reuse the gene-count likelihood
    reducer. KL describes the gene-consensus projection, not inferred allele
    frequencies. Residue enrichment lower bounds also include every compatible
    or unresolved control gene in the denominator. No smoothing or position sum.
    """
    if mapping.select('v_gene', *_AXES).is_duplicated().any() or mapping.is_empty():
        raise ValueError('unique nonempty germline projection map required')
    reference = usage.filter(pl.col('segment') == 'V').select(*_CONTEXT, 'length_scope', 'junction_length',
        'gene', 'reference_records', 'junction_gene_keys')
    background = control.filter(pl.col('segment') == 'V').select('segment', 'length_scope', 'junction_length',
        'gene', 'control_events', 'control_junction_gene_keys')
    strata = [*_CONTEXT, 'length_scope', 'junction_length']
    null_strata = ['segment', 'length_scope', 'junction_length']
    for axis, part in mapping.group_by(_AXES, maintain_order=True):
        lookup = part.select(pl.col('v_gene').alias('_source_gene'), pl.col('residue').alias('gene'))
        def project(frame, groups, counts, total_names):
            result = frame.rename({'gene': '_source_gene'}).join(lookup, on='_source_gene', how='left', validate='m:1')
            result = result.group_by(*groups, 'gene').agg(*[pl.col(c).sum() for c in counts])
            return result.with_columns(*[pl.col(c).sum().over(groups).alias(t) for c,t in zip(counts,total_names)])
        observed = project(reference, strata, ['reference_records', 'junction_gene_keys'],
            ['total_reference_records', 'total_junction_gene_keys']).with_columns(
                (pl.col('reference_records')/pl.col('total_reference_records')).alias('record_fraction'),
                (pl.col('junction_gene_keys')/pl.col('total_junction_gene_keys')).alias('key_fraction'))
        null = project(background, null_strata, ['control_events', 'control_junction_gene_keys'],
            ['total_control_events', 'total_control_junction_gene_keys'])
        alphabet = ''.join(sorted((set(IMGT_CLASSES) if axis[2] == 'full' else set(IMGT_CLASSES.values())) | {'.'}))
        possible = background.rename({'gene': 'v_gene'}).join(part.select('v_gene', 'modeled_residues', 'projection_status'),
            on='v_gene', how='left', validate='m:1').with_columns(
                pl.when(pl.col('projection_status').is_in(['resolved', 'allele_dependent']).fill_null(False))
                .then(pl.col('modeled_residues')).otherwise(pl.lit(alphabet)).str.split('').alias('gene'))
        possible = possible.explode('gene').group_by(*null_strata, 'gene').agg(
            pl.col('control_junction_gene_keys').sum().alias('control_possible_v_keys'))
        compared = compare_gene_usage(observed, null).with_columns(*[pl.lit(v).alias(k) for k,v in zip(_AXES,axis)])
        points, info = _selection_information(compared, [*_CONTEXT, *_AXES])
        points = points.join(possible, on=[*null_strata, 'gene'], how='left', nulls_equal=True,
            validate='m:1', maintain_order='left').with_columns(pl.when(pl.col('gene').is_not_null() &
                pl.col('total_control_junction_gene_keys').is_not_null()).then(
                    pl.col('control_possible_v_keys').fill_null(0)).otherwise(None).alias('control_possible_v_keys')).with_columns(
                (pl.col('control_possible_v_keys')/pl.col('total_control_junction_gene_keys')).alias('control_possible_fraction'))
        points = points.with_columns(pl.when(pl.col('gene').is_not_null() & (pl.col('control_possible_fraction') > 0))
            .then((pl.col('key_fraction')/pl.col('control_possible_fraction')).log(2)).otherwise(None).alias('residue_bits_lower'),
            pl.when(pl.col('gene').is_null()).then(pl.lit('unresolved_reference_residue'))
            .when(pl.col('total_control_junction_gene_keys').is_null()).then(pl.lit('no_control_length'))
            .when(pl.col('control_possible_v_keys')==0).then(pl.lit('residue_not_observed_in_control'))
            .otherwise(pl.lit('available')).alias('lower_bound_status'))
        renames = {'gene':'residue', 'junction_gene_keys':'reference_v_keys', 'total_junction_gene_keys':'reference_stratum_v_keys',
            'control_junction_gene_keys':'control_v_keys', 'total_control_junction_gene_keys':'control_stratum_v_keys',
            'key_fraction':'reference_fraction', 'control_key_fraction':'control_fraction', 'gene_bits':'consensus_bits', 'gene_bits_imputed':'consensus_bits_imputed',
            'control_status':'consensus_control_status'}
        points = points.rename(renames).with_columns(pl.col('consensus_control_status').replace(
            {'unresolved_reference_gene':'unresolved_reference_residue', 'gene_not_observed':'residue_not_observed'}))
        columns = [*_CONTEXT, *_AXES, 'length_scope', 'junction_length', 'residue', 'reference_records',
            'reference_v_keys', 'reference_stratum_v_keys', 'control_v_keys', 'control_stratum_v_keys',
            'reference_fraction', 'control_fraction', 'consensus_bits', 'control_possible_v_keys', 'control_possible_fraction',
            'residue_bits_lower', 'length_bits', 'joint_bits', 'consensus_control_status', 'lower_bound_status', 'consensus_bits_imputed', 'length_bits_imputed']
        yield points.select(columns).sort(*_CONTEXT, 'length_scope', 'junction_length', 'residue', nulls_last=True), info.rename({
            'gene_information_bits':'consensus_information_bits', 'reference_keys':'reference_v_keys', 'control_keys':'control_v_keys',
            'unresolved_comparison_keys':'unavailable_comparison_v_keys', 'unavailable_length_keys':'unavailable_length_v_keys'})
