# -*- coding: utf-8 -*-
#
# rf documentation build configuration file, created by
# sphinx-quickstart on Mon Nov  4 11:31:39 2013.
#
# This file is execfile()d with the current directory set to its containing dir.
# build with sphinx-build -aEn . _build

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
             ]


autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
}

default_role = 'py:obj'
templates_path = ['_templates']
exclude_patterns = ['_build']

project = 'rf'
copyright = '2013-2024, Tom Eulenfeld'

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

html_theme = 'furo'
html_logo = '_static/logo_rf.svg'
html_show_sphinx = True
html_static_path = ['_static']
html_title = f'Receiver function calculation in seismology <br>(v{version} docs)'
html_theme_options = {
    'top_of_page_buttons': [],
}

# Output file base name for HTML help builder.
htmlhelp_basename = 'rfdoc'

# Configuration for intersphinx
intersphinx_mapping = {'obspy': ('https://docs.obspy.org/', None),
                       #'python': ('https://docs.python.org/2.7/', None)
                       }
