# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SIRENA'
copyright = '2024, B. Cobo, M.T. Ceballos'
author = 'B. Cobo, M.T. Ceballos'
release = '10.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.duration',]
#extensions = ['sphinx.ext.duration','sphinxcontrib.bibtex']

templates_path = ['_templates']
exclude_patterns = []
#bibtex_bibfiles = ['references.bib']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_theme = 'classic'
html_static_path = ['_static']
html_logo = 'images/SIRENA_black.png'
