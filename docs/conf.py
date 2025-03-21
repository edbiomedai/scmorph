# Configuration file for the Sphinx documentation builder.

# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------
import sys
from typing import TYPE_CHECKING
from pathlib import Path, PurePosixPath
from datetime import datetime
from importlib.metadata import metadata

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE / "extensions"))

if TYPE_CHECKING:
    pass


# -- Project information -----------------------------------------------------

# NOTE: If you installed your project in editable mode, this might be stale.
#       If this is the case, reinstall it to refresh the metadata
info = metadata("scmorph")
project_name = info["Name"]
author = info["Author"]
copyright = f"{datetime.now():%Y}, {author}"
version = info["Version"]
urls = dict(pu.split(", ") for pu in info.get_all("Project-URL"))
repository_url = urls["Source"]

# The full version, including alpha/beta/rc tags
release = info["Version"]

bibtex_bibfiles = ["references.bib"]
bibtex_reference_style = "author_year"

templates_path = ["_templates"]
nitpicky = True  # Warn about broken links
needs_sphinx = "4.0"

html_context = {
    "display_github": True,  # Integrate GitHub
    "github_user": "edbiomedai",
    "github_repo": project_name,
    "github_version": "main",
    "conf_py_path": "/docs/",
    "accent_color": "#B6C6E6",
}

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings.
# They can be extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    "myst_nb",
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinxcontrib.bibtex",
    "sphinx_autodoc_typehints",
    "sphinx_tabs.tabs",
    "sphinx.ext.mathjax",
    "IPython.sphinxext.ipython_console_highlighting",
    "sphinxext.opengraph",
    "git_ref",  # needs to be before scanpydoc.rtd_github_links
    "scanpydoc",  # needs to be before sphinx.ext.linkcode
    "sphinx.ext.linkcode",
    "matplotlib.sphinxext.plot_directive",
    *[p.stem for p in (HERE / "extensions").glob("*.py")],
]

autosummary_generate = True
autodoc_member_order = "groupwise"
default_role = "literal"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_custom_sections = [("Params", "Parameters")]
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
myst_heading_anchors = 6  # create anchors for h1-h6
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
    "html_admonition",
]
myst_url_schemes = ("http", "https", "mailto")
nb_output_stderr = "remove"
nb_execution_mode = "off"
nb_merge_streams = True
typehints_defaults = "braces"

source_suffix = {
    ".rst": "restructuredtext",
    ".ipynb": "myst-nb",
    ".myst": "myst-nb",
}

intersphinx_mapping = {
    "anndata": ("https://anndata.readthedocs.io/en/latest/", None),
    "scanpy": ("https://scanpy.readthedocs.io/en/latest/", None),
    "pyod": ("https://pyod.readthedocs.io/en/latest/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "seaborn": ("https://seaborn.pydata.org/", None),
}


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]

html_title = project_name

html_show_sphinx = False
html_logo = "_static/img/scmorph_logo_embeddedFont.svg"

html_theme_options = {
    "repository_url": repository_url,
    "use_repository_button": True,
    "path_to_docs": "docs/",
    "navigation_with_keys": True,
}

ogp_social_cards = {
    "image": "_static/img/scmorph_logo.png",
}


pygments_style = "default"

nitpick_ignore = [
    ("py:class", "matplotlib.*"),
    ("py:class", "seaborn.matrix.FacetGrid"),
    ("py:class", "scmorph.qc.images._Classifier"),
]

# Options for plot output
plot_include_source = True
plot_formats = [("png", 150)]
plot_html_show_formats = False
plot_html_show_source_link = False
plot_working_directory = HERE.parent  # Project root

# link config
rtd_links_prefix = PurePosixPath("src")
