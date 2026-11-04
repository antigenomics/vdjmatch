"""Synthetic single-chain and paired label lookups; no download or scientific evaluation."""
import polars as pl

from vdjmatch import Annotator

single = Annotator.from_labels([('CASSQYF', 'class-a'), ('CATSRYF', 'class-b')], locus='TRB')
calls = single.annotate(['CASSQYF', 'CWWWWWWWF'])
assert calls['vdjmatch_epitope'].to_list() == ['class-a', None]
print(calls)

paired = Annotator.from_labels([('CAVQYF', 'CASSQYF', 'class-a'),
                               ('CAARYF', 'CATSRYF', 'class-b')])
queries = pl.DataFrame({'junction_aa_alpha': ['CAVQYF', 'CAARYF'],
                        'junction_aa_beta': ['CASSQYF', 'CATSRYF']})
calls = paired.annotate_paired(queries)
assert calls['vdjmatch_epitope'].to_list() == ['class-a', 'class-b']
print(calls)
