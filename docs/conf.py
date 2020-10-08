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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

from sphinx_astropy.conf import *
#from sphinx_astropy.conf.v1 import *

# -- Project information -----------------------------------------------------

project = 'Spec_pipeline'
copyright = '2020, Roberto Assef'
author = 'Roberto Assef'

# The full version, including alpha/beta/rc tags
release = '0.9.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
# extensions = [
#     'sphinx.ext.autodoc',
#     'sphinx.ext.autosummary',
#     # 'sphinx.ext.doctest',
#     # 'sphinx.ext.intersphinx',
#     # 'sphinx.ext.todo',
#     'numpydoc',
#     # 'sphinx.ext.ifconfig',
#     # 'sphinx.ext.viewcode',
#     # 'sphinx.ext.intersphinx',
#     # 'sphinx.ext.todo',
#     # 'sphinx.ext.coverage',
#     # 'sphinx.ext.pngmath',
#     # 'sphinx.ext.inheritance_diagram',
#     # 'astropy.sphinx.ext.numpydoc',
#     # 'astropy.sphinx.ext.astropyautosummary',
#     # 'astropy.sphinx.ext.automodsumm',
#     # 'astropy.sphinx.ext.automodapi',
#     # 'astropy.sphinx.ext.tocdepthfix',
#     # 'astropy.sphinx.ext.doctest',
#     # 'astropy.sphinx.ext.changelog_links',
#     # 'astropy.sphinx.ext.viewcode',  # Use patched version of viewcode
#     # 'astropy.sphinx.ext.smart_resolver'
# ]

# autoclass_content = 'both'
#
# autodoc_default_options = {
#     #'members': True,
#     #'member-order': 'bysource',
#     #'special-members': '__init__',
#     #'undoc-members': True,
#     #'exclude-members': '__weakref__'
# }

#autosectionlabel_prefix_document = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'alabaster'
# html_theme = 'bootstrap-astropy'
#
# html_theme_options = {
#     'logotext1': 'packagename',  # white,  semi-bold
#     'logotext2': '',  # orange, light
#     'logotext3': ':docs',   # white,  light
#     'astropy_project_menubar': True
#     }

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']
