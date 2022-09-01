import os
import sys
from datetime import datetime

import scmorph

sys.path.insert(0, os.path.abspath("../../src/scmorph/"))
nitpicky = True  # Warn about broken links.
needs_sphinx = "2.0"  # Nicer param docs

# -- Project information -----------------------------------------------------

project = "scmorph"
copyright = f"{datetime.now():%Y}, Jesko Wagner"
author = "Jesko Wagner"
version = scmorph.__version__

# -- General configuration ---------------------------------------------------

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx_autodoc_typehints",
    "nbsphinx",
    "scanpydoc.rtd_github_links",
    "scanpydoc.theme",
    "scanpydoc.definition_list_typed_field",
    "scanpydoc.autosummary_generate_imported",
]

source_suffix = [".rst", ".md"]

templates_path = ["_templates"]

exclude_patterns = [""]

# API documentation
autosummary_generate = True
autodoc_member_order = "bysource"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [("Params", "Parameters")]
typehints_defaults = "braces"
nitpick_ignore = [("py:class", "type")]

# -- Options for HTML output -------------------------------------------------

html_theme = "scanpydoc"
html_theme_options = dict(navigation_depth=4, logo_only=True, docsearch_index="scmorph")

html_context = dict(
    display_github=True,  # Integrate GitHub
    github_user="jeskowagner",  # Username
    github_repo="scmorph",  # Repo name
    github_version="master",  # Version
    conf_py_path="/docs/",  # Path in the checkout to the docs root
)
html_static_path = ["_static"]
html_show_sphinx = False
html_logo = "_static/img/scmorph_logo_white.svg"

html_static_path = ["_static"]
