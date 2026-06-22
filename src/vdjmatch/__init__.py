"""vdjmatch — fast, control-calibrated annotation of TCR antigen specificity.

Built on the ``seqtree`` fuzzy-search core: annotates AIRR clonotypes against VDJdb,
reporting control-calibrated E-values and enriched antigen labels. See ROADMAP in the
seqtree repo and the project README for the design.
"""

__version__ = "0.0.1"

from .api import Annotator, annotate  # noqa: E402  high-level annotation API

__all__ = ["Annotator", "annotate", "__version__"]
