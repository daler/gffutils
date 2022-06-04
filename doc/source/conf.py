# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../..'))
import gffutils
import sphinx_rtd_theme

# -- Project information -----------------------------------------------------

project = 'gffutils'
copyright = '2013-2022, Ryan Dale'
author = 'Ryan Dale'

# The full version, including alpha/beta/rc tags
release = '0.11'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.doctest',
    'sphinx_rtd_theme',
    'autoapi.extension',
    'sphinx.ext.autodoc',
]

doctest_global_setup = """
import gffutils
"""
default_role = "literal"
version = gffutils.__version__
release = gffutils.__version__
autoapi_dirs = ['../../gffutils']
autoapi_generate_api_docs = False
pygments_style = 'sphinx'
issue_tracker_url = 'https://github.com/daler/gffutils/issues/{issue}'
templates_path = ['_templates']
exclude_patterns = []
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
