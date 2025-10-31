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

# --- Make Sphinx import the local repo code first ---
HERE = os.path.dirname(__file__)                         # docs/source
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))  # repo root
sys.path.insert(0, REPO_ROOT)

print("Sphinx autodoc will import from REPO_ROOT:", REPO_ROOT)

#sys.path.insert(0, os.path.abspath('source/'))
sys.path.insert(0, os.path.abspath('../'))


# -- Project information -----------------------------------------------------

project = 'CompactObject'
copyright = '2024, Chun Huang, João Cartaxo, Tuhin Malik, Wenli Yuan,Shashwat Sourav, Xuezhi Liu, Tianzhe Zhou, John Groger, Zhiheng Wang, Nicole Osborn, Nathan Whitsett'
author = 'Chun Huang, João Cartaxo, Tuhin Malik, Wenli Yuan, Shashwat Sourav, Xuezhi Liu, Tianzhe Zhou, John Groger, Zhiheng Wang, Nicole Osborn, Nathan Whitsettt'

# The full version, including alpha/beta/rc tags
version = '1.9.9'
release = '1.9.9'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
                'sphinx.ext.napoleon',
                'sphinx.ext.intersphinx',
                'sphinx.ext.viewcode',
                'sphinx.ext.coverage',
                'sphinx.ext.mathjax',
                'sphinx.ext.githubpages',
                'sphinx.ext.autosummary',
                'nbsphinx'
]

intersphinx_mapping = {'sphinx': ('http://www.sphinx-doc.org/en/master', None),
           'numpy': ('https://docs.scipy.org/doc/numpy', None),
           'UltraNest': ('https://johannesbuchner.github.io/UltraNest/index.html', None),
           'scipy': ('https://docs.scipy.org/doc/', None),
           'CodeAstro': ('https://semaphorep.github.io/codeastro/', None)
           }

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
source_suffix = ['.rst','.md']
master_doc = 'index'
language = 'en'
pygments_style = 'sphinx'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_sidebars = {
    '**': [
        'relations.html',  # needs 'show_related': True theme option to display
        'searchbox.html',
    ]
}

# If true, links to the reST sources are added to the pages.
html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
html_show_copyright = True

# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'CompactObjectdoc'

latex_documents = [
    (master_doc, 'CompactObject.tex', 'CompactObject Documentation',
     'Chun Huang', 'manual'),
]

# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'CompactObject', 'CompactObject Documentation',
     [author], 1)
]

# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'CompactObject', 'CompactObject Documentation',
     author, 'CompactObject', 'One line description of project.',
     'Miscellaneous'),
]


