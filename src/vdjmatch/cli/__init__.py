"""Shared command-line argument conversion."""


def cpu_budget(value):
    """Resolve automatic CPU counts before standalone subcommand dispatch."""
    from argparse import ArgumentTypeError
    from vdjtools.cores import kernel_threads

    try:
        return kernel_threads(int(value))
    except (TypeError, ValueError) as exc:
        raise ArgumentTypeError(f"threads: {exc}") from exc
