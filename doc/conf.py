# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Radio calibration toolkit'
copyright = '2023, Tomas Fodran'
author = 'Tomas Fodran'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration


extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinx.ext.autodoc',
    'sphinx.ext.githubpages',
    'myst_parser',
  #  "sphinx.ext.napoleon",
    'nbsphinx',
    'sphinx_rtd_size'
]

# source_suffix = '.rst'
source_suffix = ['.rst', '.md']

templates_path = ['_templates', ]
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'html']

language = 'en'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = 'alabaster'
html_static_path = ['_static']

import sphinx_rtd_theme

html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

html_logo = '_static/logo.png'
html_favicon = '_static/favicon.ico'

html_theme_options = {
#    'logo_only': True,
    'display_version': True,
#    'navigation_depth': 4,
    'sticky_navigation': True,
    'style_external_links': True,
    'prev_next_buttons_location': 'bottom',
    'style_nav_header_background': '#2980B9'
}

html_css_files = [
    'custom.css',
]

html_js_files = [
    'custom.js',
]

# -- Options for todo extension ----------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/todo.html#configuration

todo_include_todos = True


##############

cmd_line_template = "sphinx-apidoc --module-first -f -o {outputdir} {moduledir}"

import os
import sys
sys.path.insert(0, os.path.abspath('../../radiocalibrationtoolkit/'))
from radiocalibrationtoolkit import *
#nbsphinx_execute = 'never'
add_module_names = False

sphinx_rtd_size_width = "70%"

