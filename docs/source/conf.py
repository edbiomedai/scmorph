import sys
from datetime import datetime
from pathlib import Path

import scmorph

HERE = Path(__file__).parent

sys.path[:0] = [str(HERE.parent)]

for generated in HERE.glob("anndata.*.rst"):
    generated.unlink()

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
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",
    "nbsphinx",
    "matplotlib.sphinxext.plot_directive",
    "scanpydoc",
    "scanpydoc.rtd_github_links",
    "scanpydoc.theme",
    "scanpydoc.definition_list_typed_field",
    "scanpydoc.autosummary_generate_imported",
]

source_suffix = [".rst", ".md"]

templates_path = ["_templates"]

exclude_patterns = [""]

intersphinx_mapping = dict(
    anndata=("https://anndata.readthedocs.io/en/latest/", None),
    scanpy=("https://scanpy.readthedocs.io/en/latest/", None),
    pyod=("https://pyod.readthedocs.io/en/latest/", None),
)

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
nitpick_ignore = [("py:class", "matplotlib.colors.Colormap")]

# -- Options for HTML output -------------------------------------------------

html_theme = "scanpydoc"
html_theme_options = dict(navigation_depth=4, logo_only=True, docsearch_index="scmorph")

html_context = dict(
    display_github=True,  # Integrate GitHub
    github_user="edbiomedai",  # Username
    github_repo="scmorph",  # Repo name
    github_version="main",  # Version
    conf_py_path="/docs/",  # Path in the checkout to the docs root
)
html_static_path = ["_static"]
html_show_sphinx = False
html_logo = "_static/img/scmorph_logo_white.svg"

html_static_path = ["_static"]

# Options for plot output

plot_include_source = True
plot_formats = [("png", 90)]
plot_html_show_formats = False
plot_html_show_source_link = False
plot_working_directory = HERE.parent  # Project root
