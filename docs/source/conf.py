# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

# explicitly assign the index page
master_doc = 'index'

# -- Project information -----------------------------------------------------

project = 'Paganini'
copyright = '2019, Maciej Bendkowski and Sergey Dovgal'
author = 'Maciej Bendkowski and Sergey Dovgal'

# The full version, including alpha/beta/rc tags
release = '1.1.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.extlinks'
]

# Do not attempt importing packages over the internet
autodoc_mock_imports = ['sympy', 'numpy', 'cvxpy']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# The LaTeX engine to build the docs.
latex_engine = 'xelatex'

# If True, the PDF build from the LaTeX files created by Sphinx will use xindy
# rather than makeindex.
latex_use_xindy = False

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Should we show "Created using Sphinx" in the HTML footer?
html_show_sphinx = False

mathjax_config = {
    'CommonHTML': {'linebreaks': {'automatic': True}},
    'HTML-CSS': {'linebreaks': {'automatic': True}},
    'SVG': {'linebreaks': {'automatic': True}},
}
