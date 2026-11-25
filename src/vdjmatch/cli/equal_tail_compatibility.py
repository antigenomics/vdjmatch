"""Recognize the one verified execution-only transition from the fixed 128 MiB score slab."""

# A change in any other source, scientific setting or reference still invalidates the recipe.
_BUDGET_REVISION = {
    'vdjmatch.evalue.equal_tail': (
        '213413167d56c1804e4001c7cb819e763e868871f570ba906107a19f959dd0fa',
        'ace92db47f7701b8eebd1dd65a342ce6136b0d0776613e373c2d75f107057b6c'),
    'vdjmatch.evalue.equal_tail_selection': (
        '306fe97299799ab438c92b49a1a465be4faaaf39752be156b96b5bc2c063afd1',
        '5a88f9f5057996fe1b3301c7f10ee0a895a91c8a7d2bd1f6e990323737f7dc14'),
    'vdjmatch.cli.equal_tail': (
        '400a4b4f3c546fd835c7e60d3832001e4ef737c53145bbc96d3fd6724e82d3a7',
        '2818f9dd6f7b373eb31a0f98ed8a7d20b2f1b10d987ab975c15defe6c7c8880f'),
    'calibration': (
        'fdf591b46bfa5be1f748e2dd03c1e89e0445dd75ac940d52ad549a49cf4a6cc6',
        'd215533729045bca7bb50af22aa235d1d7369b73f710187769e89d0545f8e7fc'),
}


def compatible_binding(frozen, current):
    if frozen == current:
        return True
    if not isinstance(frozen, dict) or not isinstance(frozen.get('source_sha256'), dict):
        return False
    pins = frozen['source_sha256'].copy()
    for name, (old, new) in _BUDGET_REVISION.items():
        if pins.get(name) != old or current['source_sha256'].get(name) != new:
            return False
        pins[name] = new
    return (frozen | {'source_sha256': pins}) == current
