#!/usr/bin/env python
import os.path
import re

from setuptools import find_packages, setup


def find_version(*paths):
    fname = os.path.join(os.path.dirname(__file__), *paths)
    with open(fname) as fp:
        code = fp.read()
    match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", code, re.M)
    if match:
        return match.group(1)
    raise RuntimeError("Unable to find version string.")


version = find_version('rf', '__init__.py')

with open('README.rst') as f:
    README = f.read()
if 'dev' not in version:  # get image for correct version from travis-ci
    README = README.replace('branch=master', 'branch=v%s' % version)
readme = README.split('\n')
DESCRIPTION = readme[2]
LONG_DESCRIPTION = '\n'.join(readme[5:])

ENTRY_POINTS = {
    'console_scripts': ['rf-runtests = rf.tests:run',
                        'rf = rf.batch:run_cli']}

REQUIRES = ['decorator', 'matplotlib', 'numpy', 'scipy',
            'setuptools', 'obspy>=1.0',
            'cartopy', 'geographiclib', 'shapely', 'toeplitz', 'tqdm']
# optional: joblib, obspyh5
# documentation: sphinx, alabaster, obspy

setup(name='rf',
      version=version,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      url='https://github.com/trichter/rf',
      author='Tom Eulenfeld',
      author_email='tom.eulenfeld@gmail.de',
      license='MIT',
      packages=find_packages(),
      package_dir={'rf': 'rf'},
      install_requires=REQUIRES,
      entry_points=ENTRY_POINTS,
      include_package_data=True,
      zip_safe=False
      )
