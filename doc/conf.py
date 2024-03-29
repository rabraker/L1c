# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

import sys

project = "L1c"
copyright = "2019, Arnold Braker"
author = "Arnold Braker"
html_context = {"header_txt": "libl1C User Guide"}
# The short X.Y version
version = ""
# The full version, including alpha/beta/rc tags
release = "0.0"


# -- General configuration ---------------------------------------------------


extensions = [
    "breathe",
    "sphinx.ext.autodoc",
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "sphinxcontrib.fulltoc",
    "sphinxcontrib.bibtex",
]

breathe_projects = {"L1c": "/home/arnold/matlab/l1c/doc/xml/"}
breathe_default_project = "L1c"

mathjax_config = {
    "extensions": ["tex2jax.js"],
    "jax": ["input/TeX", "output/HTML-CSS"],
}


# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = None


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'default'
html_theme = "basic"
# html_theme = 'alabaster'
# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {
#     'description': "Compressed Sensing in C"
# }

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
html_sidebars = {
    "**": [
        "localtoc.html",
    ]
}

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = "L1cdoc"


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    # 'papersize': 'letterpaper',
    # The font size ('10pt', '11pt' or '12pt').
    # 'pointsize': '10pt',
    # 'preamble': '',
    # Latex figure (float) alignment
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, "L1c.tex", "L1c Documentation", "Arnold Braker", "manual"),
]


# -- Options for manual page output ------------------------------------------
# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(master_doc, "l1c", "L1c Documentation", [author], 1)]


# -- Options for Texinfo output ----------------------------------------------
# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        "L1c",
        "L1c Documentation",
        author,
        "L1c",
        "One line description of project.",
        "Miscellaneous",
    ),
]


# -- Extension configuration -------------------------------------------------

# -- Options for todo extension ----------------------------------------------

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True
