"""Sphinx configuration for the vdjmatch documentation."""
from importlib.metadata import version as package_version

project = "vdjmatch"
release = package_version(project)
rst_prolog = f".. |release| replace:: {release}\n"
author = "ISALGO laboratory"
copyright = "2026, ISALGO laboratory"

extensions = [
    "sphinx.ext.mathjax",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
    "nbsphinx",
]

# autodoc imports vdjmatch with its heavy / optional deps mocked.
autodoc_mock_imports = ["seqtree", "polars", "huggingface_hub", "numpy", "scipy", "Bio", "ezodf", "psutil"]
autodoc_typehints = "description"
autodoc_member_order = "bysource"

intersphinx_mapping = {"python": ("https://docs.python.org/3", None)}

templates_path = ["_templates"]
exclude_patterns = ["_build", "**.ipynb_checkpoints"]

html_theme = "pydata_sphinx_theme"
html_title = f"vdjmatch {release}"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "github_url": "https://github.com/antigenomics/vdjmatch",
    "navigation_with_keys": True,
    "show_prev_next": False,
    "show_nav_level": 2,
    "navigation_depth": 3,
    "collapse_navigation": False,
    "header_links_before_dropdown": 4,
    "show_toc_level": 2,
}
html_sidebars = {"**": ["site-nav"], "index": []}
nbsphinx_execute = "never"
