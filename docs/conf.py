# -*- coding: utf-8 -*-
#
# rf documentation build configuration file, created by
# sphinx-quickstart on Mon Nov  4 11:31:39 2013.
#
# This file is execfile()d with the current directory set to its containing dir.


import sys, os, traceback

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute.
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
root = os.path.abspath('../')
sys.path.insert(0, root)


class Mock(object):
    def __init__(self, *args, **kwargs):
        pass

    def __call__(self, *args, **kwargs):  #@UnusedVariable
        return Mock()

    @classmethod
    def __getattr__(cls, name):
        if name in ('__file__', '__path__'):
            return '/dev/null'
        elif name[0] == name[0].upper():
            mockType = type(name, (), {})
            mockType.__module__ = __name__
            return mockType
        else:
            return Mock()

# Mock all modules in rf which raise an import error
for i in range(20):
    try:
        import rf
        import rf.batch
        import rf.imaging
    except ImportError:
        exc_type, exc_value, tb = sys.exc_info()
        codeline = traceback.extract_tb(tb)[-1][-1]
        missing_module = codeline.split()[1]
        if 'rf' in missing_module:
            raise
        # mock missing module and all parent modules
        for c in range(missing_module.count('.'), -1, -1):
            m = missing_module.rsplit('.', c)[0]
            if sys.modules.get(m, None) is None:
                print('Mocking module %s' % m)
                sys.modules[m] = Mock()
    else:
        break


def get_all_param_type_declarations():
    from glob import glob
    from re import search
    ignobj = []
    text = ''
    for fname in glob(os.path.join('..', 'rf', '*.py')):
        with open(fname) as f:
            text = text + f.read()
    for line in text.splitlines():
        match = search(':param\s+(\S+)\s+\S+:', line)
        if match is None:
            match = search(':type.*:\s*(.*?)\s*$', line)
        if match is not None:
            ignobj.append(match.group(1))
    return ignobj


# Show warnings for unreferenced targets
nitpicky = True
ignobj = get_all_param_type_declarations()
igncls = ('object', 'exceptions.Exception', 'json.decoder.JSONDecoder')
nitpick_ignore = ([('py:obj', n) for n in ignobj] +
                  [('py:class', n) for n in igncls])

# Add any Sphinx extension module names here
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.doctest',
              'sphinx.ext.intersphinx',
              'sphinx.ext.autosummary',
              'sphinx.ext.viewcode',
              'alabaster']

autodoc_default_flags = ['members', 'undoc-members', 'show-inheritance']
default_role = 'py:obj'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build']

# The suffix of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'rf'
copyright = u'2013-2016, Tom Eulenfeld'

# The full version, including alpha/beta/rc tags.
release = rf.__version__
# The short X.Y version.
version = release.rsplit(".", 1)[0]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# -- Options for HTML output ---------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#html_theme = 'nature'  #nature, sphinxdoc

import alabaster

html_theme_path = [alabaster.get_path()]
html_theme = 'alabaster'
html_sidebars = {
    '**': [
        'about.html',
        'mynavigation.html',
        'searchbox.html',
    ]
}

html_theme_options = {
    'logo': 'logo_rf.svg',
    'github_user': 'trichter',
    'github_repo': 'rf',
    'description': 'Receiver function calculation in seismology',
    'show_powered_by': False,
    'page_width': '1240px',
    #'sidebar_width': '220px',
    'extra_nav_links': {'Project on Github': 'https://github.com/trichter/rf'}

}

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
#html_logo = None

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#html_favicon = ''

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
#html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
#html_show_copyright = True

# Output file base name for HTML help builder.
htmlhelp_basename = 'rfdoc'


# -- Options for LaTeX output --------------------------------------------------

latex_elements = {} # stuff for the LaTeX preamble

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
latex_documents = [
  ('index', 'rf.tex', u'rf Documentation',
   u'Tom Eulenfeld', 'manual'),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False


# -- Options for manual page output --------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    ('index', 'rf', u'rf Documentation',
     [u'Tom Eulenfeld'], 1)
]


# -- Options for Texinfo output ------------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
  ('index', 'rf', u'rf Documentation',
   u'Tom Eulenfeld', 'rf', 'Receiver function calculation in seismology',
   'Miscellaneous'),
]

# -- Options for Epub output ---------------------------------------------------

# Bibliographic Dublin Core info.
epub_title = u'rf'
epub_author = u'Tom Eulenfeld'
epub_publisher = u'Tom Eulenfeld'
epub_copyright = u'2013-2016, Tom Eulenfeld'


# Configuration for intersphinx
intersphinx_mapping = {'obspy': ('http://docs.obspy.org/', None),
                       #'python': ('https://docs.python.org/2.7/', None)
                       }
