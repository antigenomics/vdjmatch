"""Multi-sample annotation runner."""
__all__ = ["annotate_sample"]


def __getattr__(name):
    if name == "annotate_sample":
        from .multisample import annotate_sample
        return annotate_sample
    raise AttributeError(name)
