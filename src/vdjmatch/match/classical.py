"""Frozen manuscript alignment model, without fitting or rounded matrices.

The bundled segment scores are already self-normalized. Junction alignments
use raw double VDJAM similarities and a separate linear NW selection charge
``max(abs(VDJAM))`` per gap residue. The fitted gap coefficient is applied only
after that alignment has been selected.
"""
from __future__ import annotations

from hashlib import sha256
from importlib import resources
from io import BytesIO

import numpy as np
import polars as pl

_ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'


def gene_key(column: str) -> pl.Expr:
    """Literal gene lookup key; ambiguous calls remain unknown, not first-called."""
    value = pl.col(column).cast(pl.String).str.strip_chars().str.to_uppercase()
    return pl.when(value.is_not_null() & ~value.str.contains(','))\
        .then(value.str.split('*').list.first()).otherwise(None)


def frozen_model(species: str, locus: str) -> dict:
    """Load one bundled frozen model and its compact coefficient-weighted terms.

    Returns the full double similarity grid, V/J gene names and square offset
    matrices, ``[intercept, cdr3 coefficient, gap coefficient]`` and asset hashes.
    No caller files, fitted parameters, annotation inference or computed cache
    are used. Missing model genes receive ID -1 in the caller and remain unscored.
    """
    from ..evalue.control import _organism

    organism = _organism(species)
    name = {'human': 'HomoSapiens', 'mouse': 'MusMusculus'}[organism]
    if locus not in ('TRA', 'TRB'):
        raise ValueError('frozen manuscript model requires TRA or TRB')
    assets = {key: resources.files('vdjmatch.resources').joinpath(key).read_bytes()
              for key in ('vdjam.txt', 'score_coef.txt', 'segm_score.txt')}
    coef = pl.read_csv(BytesIO(assets['score_coef.txt']), separator='\t').filter(
        (pl.col('species') == name) & (pl.col('gene') == locus) & (pl.col('fun') == 'cloglog'))
    if coef.height != 1:
        raise ValueError('exactly one bundled cloglog coefficient row is required')
    coefficients = coef.row(0, named=True)
    terms = ['(Intercept)', 'cdr1.score', 'cdr2.score', 'cdr3.score', 'indels', 'v.score', 'j.score']
    if not np.isfinite([coefficients[key] for key in terms]).all():
        raise ValueError('bundled frozen coefficients must be finite')
    sim = pl.read_csv(BytesIO(assets['vdjam.txt']), separator='\t')
    if sim.height != 400 or sim.select('aa.1', 'aa.2').n_unique() != 400:
        raise ValueError('bundled VDJAM must contain all 400 amino-acid pairs')
    ids = {aa: i for i, aa in enumerate(_ALPHABET)}
    grid = np.full((20, 20), np.nan, dtype=np.float64)
    i = sim['aa.1'].replace_strict(ids).to_numpy()
    j = sim['aa.2'].replace_strict(ids).to_numpy()
    grid[i, j] = sim['score'].to_numpy()
    if not np.isfinite(grid).all():
        raise ValueError('bundled double VDJAM must be complete and finite')
    segment = pl.read_csv(BytesIO(assets['segm_score.txt']), separator='\t', infer_schema_length=None).filter(
        (pl.col('species') == name) & (pl.col('gene') == locus)).with_columns(
            gene_key('id.1').alias('_first'), gene_key('id.2').alias('_second'))
    names = []; offsets = []
    for kind in ('V', 'J'):
        part = segment.filter(pl.col('segment') == kind)
        available = set(part['_first']) | set(part['_second'])
        if None in available:
            raise ValueError('bundled segment gene keys cannot be missing')
        genes = sorted(available)
        if not genes or part.select('_first', '_second').n_unique() != part.height:
            raise ValueError('bundled segment gene-pair keys must be nonempty and unique')
        lookup = {gene: index for index, gene in enumerate(genes)}
        values = coefficients['v.score'] * part['segm.score'] if kind == 'V' else coefficients['j.score'] * part['segm.score']
        if kind == 'V':
            values = values + coefficients['cdr1.score'] * part['cdr1.score'] + coefficients['cdr2.score'] * part['cdr2.score']
        matrix = np.full((len(genes), len(genes)), np.nan, dtype=np.float64)
        matrix[part['_first'].replace_strict(lookup).to_numpy(),
               part['_second'].replace_strict(lookup).to_numpy()] = values.to_numpy()
        if not np.isfinite(matrix).all():
            raise ValueError('bundled segment offsets must form a complete finite matrix')
        names.append(genes); offsets.append(matrix)
    return dict(species=name, locus=locus, similarity=grid, similarity_symbols=_ALPHABET,
        selection_gap=float(np.max(np.abs(grid))), group_names=names, group_offsets=offsets,
        coefficients=np.asarray([coefficients['(Intercept)'], coefficients['cdr3.score'], coefficients['indels']], dtype=np.float64),
        source_coefficients=coefficients,
        **{'sha256_' + key.replace('.', '_'): sha256(value).hexdigest() for key, value in assets.items()})


def fitted_probability(eta):
    """Stable frozen cloglog inverse: ``-expm1(-exp(eta))`` (not a posterior)."""
    values = np.asarray(eta, dtype=np.float64)
    if not np.isfinite(values).all():
        raise ValueError('finite frozen linear predictors are required')
    with np.errstate(over='ignore', under='ignore'):
        return -np.expm1(-np.exp(values))


def alignment_probability(model, raw_similarity, gap_residues, query_self, reference_self,
                          v_offset=0., j_offset=0.):
    """Apply the frozen model to an already selected alignment's numeric terms.

    ``raw_similarity`` sums aligned residue similarities without NW gap charges;
    full-junction self sums retain conserved anchors. Segment offsets already
    include their own frozen coefficients and self normalization.
    """
    b, c3, gamma = model['coefficients']
    gaps = np.asarray(gap_residues)
    if np.any(gaps < 0) or np.any(gaps != np.floor(gaps)):
        raise ValueError('gap residues must be nonnegative integers')
    cdr3 = np.asarray(raw_similarity) - np.maximum(query_self, reference_self)
    return fitted_probability(b + np.asarray(v_offset) + np.asarray(j_offset) + c3*cdr3 + gamma*gaps)
