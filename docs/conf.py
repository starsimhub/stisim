# -*- coding: utf-8 -*-
#
# Sphinx Demo docs build configuration file, created by
# sphinx-quickstart on Fri Aug 26 16:52:16 2016.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# docs root, use os.path.abspath to make it absolute, like shown here.
#
import os
import subprocess
import sys
#import sphinx_rtd_theme
import configparser
from datetime import datetime

if sys.platform in ["linux", "darwin"]:
    subprocess.check_output(["make", "generate-api"], cwd=os.path.dirname(os.path.abspath(__file__)))
else:
    subprocess.check_output(["make.bat", "generate-api"], cwd=os.path.dirname(os.path.abspath(__file__)))

# Rename "stisim" to "API reference"
filename = 'api/modules.rst' # This must match the Makefile
with open(filename) as f: # Read existing file
    lines = f.readlines()
lines[0] = "API reference\n" # Blast away the existing heading and replace with this
lines[1] = "=============\n" # Ensure the heading is the right length
with open(filename, "w") as f: # Write new file
    f.writelines(lines)


# -- General configuration ------------------------------------------------

# If your docs needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.mathjax', # use latex to create math equations
    'sphinx.ext.githubpages', # host on gh-pages
    'sphinx.ext.autodoc', # create api reference
    'sphinx.ext.todo', 
    'plantweb.directive', # create diagrams from text files
    'sphinxcontrib.programoutput', # show program output
    'sphinx.ext.intersphinx', # link to other sphinx docsets
    'sphinxext.remoteliteralinclude', # include remote code examples
    'sphinx.ext.viewcode' # link to view source code
]


plantuml = 'plantweb'

autodoc_default_options = {
    'member-order': 'bysource',
    'members': None,
    'exclude-members': '__all__'
}

autodoc_mock_imports = []

napoleon_google_docstring = True
# napoleon_numpy_docstring = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
source_suffix = ['.rst']

# The encoding of source files.
#
# source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'stisim'
copyright = f'1999 - {datetime.today().year}, Gates Foundation. All rights reserved.'
author = u'Institute for Disease Modeling'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
#current_path = os.path.dirname(__file__)
#version_path = os.path.join(current_path, '..', '.bumpversion.cfg')
#config = configparser.ConfigParser()
#config.read(version_path)

#version = config['bumpversion']['current_version']
# The full version, including alpha/beta/rc tags.
# release = u'1.0'

# The language for content autogenerated by Sphinx. Refer to docs
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#
# today = ''
#
# Else, today_fmt is used as the format for a strftime call.
#
# today_fmt = '%B %d, %Y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The reST default role (used for this markup: `text`) to use for all
# documents.
#
# default_role = None

# If true, '()' will be appended to :func: etc. cross-reference text.
#
# add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
#
# add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
#
# show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# A list of ignored prefixes for module index sorting.
# modindex_common_prefix = []

# If true, keep warnings as "system message" paragraphs in the built documents.
# keep_warnings = False

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

# RST epilog is added to the end of every topic. Useful for replace
# directives to use across the docset.

# rst_epilog = "\n.. include:: /variables.txt"

# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the docs for
# a list of builtin themes.
#
html_theme = 'pydata_sphinx_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# docs.
# #
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 3,
    "show_prev_next": True,
    "icon_links": [
        {"name": "IDM docs", "url": "https://docs.idmod.org", "icon": "fas fa-home"},
        {
            "name": "GitHub",
            "url": "https://github.com/starsimhub/stisim",
            "icon": "fab fa-github-square",
        },
    ],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "secondary_sidebar_items": ["navbar-side"],
    "header_links_before_dropdown": 5,
    "footer_start": ["copyright", "footer_start"],
    "footer_end": ["theme-version", "footer_end"],
}
html_sidebars = {
    "**": ["sidebar-nav-bs", "page-toc"],
}
#html_logo = "images/idm-logo-transparent.png"
html_favicon = "images/favicon.ico"
html_static_path = ['_static']
html_baseurl = "https://docs.idmod.org/projects/stisim/en/latest"
html_context = {
    'rtd_url': 'https://docs.idmod.org/projects/stisim/en/latest',
    "versions_dropdown": {
        "latest": "devel (latest)",
        "stable": "current (stable)",
    },
    "default_mode": "light",
}

# Add customizations
def setup(app):
    app.add_css_file("theme_overrides.css")

# Add any paths that contain custom themes here, relative to this directory.
#html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# The name for this set of Sphinx documents.
# "<project> v<release> docs" by default.
#
# html_title = u'Sphinx Demo v1.0'

# A shorter title for the navigation bar.  Default is the same as html_title.
#
# html_short_title = None

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
#
#html_logo = "images/IDM_white.png"

# The name of an image file (relative to this directory) to use as a favicon of
# the docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#
html_favicon = "images/favicon.ico"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".

html_static_path = ['_static']

html_css_files = ['theme_overrides.css']

#html_js_files = ['show_block_by_os.js']

# Add any extra paths that contain custom files (such as robots.txt or
# .htaccess) here, relative to this directory. These files are copied
# directly to the root of the docs.
#
# html_extra_path = []

# If not None, a 'Last updated on:' timestamp is inserted at every page
# bottom, using the given strftime format.
# The empty string is equivalent to '%b %d, %Y'.
#
# html_last_updated_fmt = None

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
#
# html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
#
# html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
#
# html_additional_pages = {}

# If false, no module index is generated.
#
# html_domain_indices = True

# If false, no index is generated.
#
# html_use_index = True

# If true, the index is split into individual pages for each letter.
#
# html_split_index = False

# If true, links to the reST sources are added to the pages.
#
html_show_sourcelink = False

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
#
html_show_sphinx = False

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
#
# html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
#
html_use_opensearch = 'www.idmod.org/docs/'


# This is the file name suffix for HTML files (e.g. ".xhtml").
# html_file_suffix = None

# Language to be used for generating the HTML full-text search index.
# Sphinx supports the following languages:
#   'da', 'de', 'en', 'es', 'fi', 'fr', 'hu', 'it', 'ja'
#   'nl', 'no', 'pt', 'ro', 'ru', 'sv', 'tr', 'zh'
#
# html_search_language = 'en'

# A dictionary with options for the search language support, empty by default.
# 'ja' uses this config value.
# 'zh' user can custom change `jieba` dictionary path.
#
# html_search_options = {'type': 'default'}

# The name of a javascript file (relative to the configuration directory) that
# implements a search results scorer. If empty, the default will be used.
#
# html_search_scorer = 'scorer.js'

# Output file base name for HTML help builder.
htmlhelp_basename = 'stisim'

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'stisim-docs.tex', u'stisim',
     u'Institute for Disease Modeling', 'manual'),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#
# latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#
# latex_use_parts = False

# If true, show page references after internal links.
#
# latex_show_pagerefs = False

# If true, show URL addresses after external links.
#
# latex_show_urls = False

# Documents to append as an appendix to all manuals.
#
# latex_appendices = []

# It false, will not define \strong, \code,     itleref, \crossref ... but only
# \sphinxstrong, ..., \sphinxtitleref, ... To help avoid clash with user added
# packages.
#
# latex_keep_old_macro_names = True

# If false, no module index is generated.
#
# latex_domain_indices = True


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'stisim-docs', u'stisim-hiv',
     [author], 1)
]

# If true, show URL addresses after external links.
#
man_show_urls = True

# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'stisim-docs', u'stisim',
     author, 'Institute for Disease Modeling', 'How to use STIsim for simulating sexually transmitted infections.',
     'Miscellaneous'),
]

# Documents to append as an appendix to all manuals.
#
# texinfo_appendices = []

# If false, no module index is generated.
#
# texinfo_domain_indices = True

# How to display URL addresses: 'footnote', 'no', or 'inline'.
#
# texinfo_show_urls = 'footnote'

# If true, do not generate a @detailmenu in the "Top" node's menu.
#
# texinfo_no_detailmenu = False


# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'python': ('https://docs.python.org/3', None),
                       'starsim': ('https://docs.idmod.org/projects/starsim/en/latest/', None),
                       'fpsim': ('https://docs.idmod.org/projects/fpsim/en/latest/', None),
                       'hpvsim': ('https://docs.idmod.org/projects/hpvsim/en/latest/', None),
                       }



