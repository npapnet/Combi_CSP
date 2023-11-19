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
sys.path.insert(0, os.path.abspath('../CombiCSP'))


def get_version(rel_path):
    """Get the version string from a file.
    
    Assuming the version line is in the form: __version__ = '0.1.0'
    strips out the version and remove leading and trailing whitespace and quotes
               
    Args:
        rel_path (str): The relative path to the file.

    Raises:
        RuntimeError: If the version string is not found.

    Returns:
        str: The version string.
    """
    with open(rel_path, 'r', encoding='utf-8') as fp:
        for line in fp:
            if line.startswith('__version__'):
                delim = '"' if '"' in line else "'"
                return line.split(delim)[1]
    raise RuntimeError("Unable to find version string.")

__version__ = get_version('../../CombiCSP/__init__.py')


# -- Project information -----------------------------------------------------

project = 'CombiCSP'
copyright = '2023, N. Papadakis, G. Arnaoutakis'
author = 'N. Papadakis, G. Arnaoutakis'

# The full version, including alpha/beta/rc tags
release = 'v1.5.0a0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
    'sphinxcontrib.mermaid',
]
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
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']