"""Scientific defaults use meaningful replication; small test recipes stay explicit."""
import argparse
import inspect


def test_scientific_resampling_defaults_and_profile_override():
    from vdjmatch.aggregate.profile_noise import profile_noise
    from vdjmatch.aggregate.reference_coverage import abundance_coverage
    from vdjmatch.evalue.combined import calibrate_combined
    from vdjmatch.cli.reference_profiles import register
    from vdjmatch.cli.diagnostics import register as diagnostics
    assert inspect.signature(profile_noise).parameters['shuffles'].default == 999
    assert inspect.signature(abundance_coverage).parameters['nboot'].default == 200
    assert inspect.signature(calibrate_combined).parameters['replicates'].default == 200
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest='cmd')
    register(sub)
    diagnostics(sub)
    profile = parser.parse_args(['reference-profiles', '--vdjdb', 'raw', '--output-prefix', 'out'])
    assert profile.shuffles == 999
    profile = parser.parse_args(['reference-profiles', '--vdjdb', 'raw', '--output-prefix', 'out', '--shuffles', '200'])
    assert profile.shuffles == 200
    args = parser.parse_args(['reference-diagnostics', '--vdjdb', 'raw', '--locus', 'TRB', '--output-prefix', 'out'])
    assert args.mock_replicates == args.nboot == 200
