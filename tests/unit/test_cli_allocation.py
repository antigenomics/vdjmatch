"""CLI budgets resolve allocation before dispatch, without running searches."""
import argparse
import os
from unittest.mock import patch

import pytest

from vdjmatch.cli.__main__ import main


@pytest.mark.parametrize('arguments', [
    ['annotate', 'sample.tsv'],
    ['annotate', 'sample.tsv', '--control', 'control.tsv'],
    ['neighbours', 'sample.tsv', '--vdjdb', 'ref.tsv', '--locus', 'TRB', '--output-prefix', 'out'],
    ['evidence', 'sample.tsv', '--vdjdb', 'ref.tsv', '--locus', 'TRB', '--output-prefix', 'out'],
    ['search', 'sample.tsv', '--reference', 'ref.tsv', '--sequence-col', 'junction_aa', '--output-prefix', 'out'],
    ['match', 'sample.tsv'], ['index'], ['update'],
])
def test_cli_allocation_budget(arguments):
    original = argparse.ArgumentParser.parse_args
    observed = []

    def parse(parser, *args, **kwargs):
        options = original(parser, *args, **kwargs)
        options.func = lambda options: observed.append((options.threads, os.environ['POLARS_MAX_THREADS'])) or 0
        return options

    with patch('vdjtools.cores.available_cores', return_value=3), patch.object(
            argparse.ArgumentParser, 'parse_args', parse), patch.dict(os.environ, POLARS_MAX_THREADS='1'):
        assert main(arguments) == 0
        assert observed == [(3, '1')]
        assert main(arguments + ['--threads', '0']) == 0
        with pytest.raises(SystemExit) as error:
            main(arguments + ['--threads', '4'])
        assert error.value.code == 2
        assert os.environ['POLARS_MAX_THREADS'] == '1'


def test_match_config_threads_and_initialized_pool(tmp_path):
    config = tmp_path / 'config.json'
    config.write_text('{"threads": 2}')
    original = argparse.ArgumentParser.parse_args
    observed = []

    def parse(parser, *args, **kwargs):
        options = original(parser, *args, **kwargs)
        options.func = lambda options: observed.append(options.threads) or 0
        return options

    import polars as pl
    with patch('vdjtools.cores.available_cores', return_value=3), patch.object(
            argparse.ArgumentParser, 'parse_args', parse), patch.dict(os.environ, POLARS_MAX_THREADS='1'):
        assert main(['match', '--config', str(config), 'sample.tsv']) == 0
        assert observed == [2]
        with patch.object(pl, 'thread_pool_size', return_value=3), pytest.raises(SystemExit) as error:
            main(['match', '--threads', '2', 'sample.tsv'])
        assert error.value.code == 2
        assert observed == [2]
        assert os.environ['POLARS_MAX_THREADS'] == '1'
