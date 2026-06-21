"""Sphinx configuration for the vdjmatch documentation."""

project = "vdjmatch"
author = "ISALGO laboratory"
copyright = "2026, ISALGO laboratory"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "nbsphinx",
]

# autodoc (when API pages are added) imports vdjmatch with its heavy deps mocked.
autodoc_mock_imports = ["seqtree", "polars"]
autodoc_typehints = "description"

intersphinx_mapping = {"python": ("https://docs.python.org/3", None)}

templates_path = ["_templates"]
exclude_patterns = ["_build", "**.ipynb_checkpoints"]

html_theme = "pydata_sphinx_theme"
html_title = "vdjmatch"
html_theme_options = {
    "github_url": "https://github.com/antigenomics/vdjmatch",
    "navigation_with_keys": True,
}
nbsphinx_execute = "never"
