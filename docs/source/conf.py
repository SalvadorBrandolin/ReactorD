# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import pathlib
import sys

CURRENT_PATH = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))
REACTORD_PATH = CURRENT_PATH.parent.parent

sys.path.insert(0, str(REACTORD_PATH))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "nbsphinx",
    # "sphinxcontrib.bibtex",
    "sphinx_copybutton",
]

templates_path = ["_templates"]
exclude_patterns = []

# =============================================================================
# EXTRA CONF
# =============================================================================

autodoc_member_order = "bysource"

# =============================================================================
# NUMPY DOC
# =============================================================================

numpydoc_class_members_toctree = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# General information about the project.
project = "ReactorD"

copyright = (
    "Brandolin, Salvador Eduardo - "
    "Santos, Maricel Del Valle - "
    "Parodi, Adrian - "
    "Rovezzi, Juan Pablo - "
    'Scilipoti, Jose Antonio - Copyright (c) 2022"'
)

author = (
    "Brandolin, Salvador Eduardo; "
    "Parodi, Adrian; "
    "Rovezzi, Juan Pablo; "
    "Santos, Maricel Del Valle; "
    "Scilipoti, Jose Antonio"
)

release = "0.0.1a1"

# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

import sphinx_rtd_theme

html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_favicon = "_static/favicon.ico"
