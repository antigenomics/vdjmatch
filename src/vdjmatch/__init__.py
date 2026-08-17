"""vdjmatch — fast, control-calibrated annotation of TCR antigen specificity.

Built on the ``seqtree`` fuzzy-search core: annotates AIRR clonotypes against VDJdb,
reporting control-calibrated E-values and enriched antigen labels. See ROADMAP in the
seqtree repo and the project README for the design.

``vdjmatch.precursor`` (T-cell precursor frequency for an epitope) needs the optional
``vdjmatch[precursor]`` extra and is deliberately **not** imported here — importing vdjmatch must
never pull in the recombination model.
"""
from importlib.metadata import PackageNotFoundError, version as _version

try:                                          # single source of truth: the installed distribution
    __version__ = _version("vdjmatch")
except PackageNotFoundError:                  # pragma: no cover - source tree without an install
    __version__ = "0+unknown"

from .api import Annotator, annotate  # noqa: E402  high-level annotation API

__all__ = ["Annotator", "annotate", "__version__"]
