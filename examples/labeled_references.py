"""Synthetic single-chain and paired label lookups; no download or scientific evaluation."""
from vdjmatch import Annotator

import polars as pl

single = Annotator.from_labels([('CASSQYF', 'class-a'), ('CATSRYF', 'class-b')], locus='TRB')
calls = single.annotate(['CASSQYF', 'CWWWWWWWF'])
assert calls['vdjmatch_epitope'].to_list() == ['class-a', None]
print(calls.write_csv(separator="\t"), end="")

paired = Annotator.from_labels([('CAVQYF', 'CASSQYF', 'class-a'),
                               ('CAARYF', 'CATSRYF', 'class-b')])
queries = pl.DataFrame({'junction_aa_alpha': ['CAVQYF', 'CAARYF'],
                        'junction_aa_beta': ['CASSQYF', 'CATSRYF']})
calls = paired.annotate_paired(queries)
assert calls['vdjmatch_epitope'].to_list() == ['class-a', 'class-b']
print(calls.write_csv(separator="\t"), end="")
