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


VERSION = find_version('rf', '__init__.py')
DESCRIPTION = 'Receiver function calculation in seismology'
STATUS = """|buildstatus|

.. |buildstatus| image:: https://api.travis-ci.org/trichter/rf.png?
    branch=master
   :target: https://travis-ci.org/trichter/rf"""

ENTRY_POINTS = {
    'console_scripts': ['rf-runtests = rf.tests:run',
                        'rf = rf.batch:run_cli']}

REQUIRES = ['decorator', 'matplotlib', 'numpy', 'scipy',
            'setuptools', 'obspy>=1.0',
            'cartopy', 'geographiclib', 'shapely', 'toeplitz', 'tqdm']

# optional: joblib, obspyh5
# documentation: sphinx, alabaster, obspy

CLASSIFIERS = [
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.7',
    'Topic :: Scientific/Engineering :: Physics'
    ]

if 'dev' not in VERSION:  # get image for correct version from travis-ci
    STATUS = STATUS.replace('branch=master', 'branch=v%s' % VERSION)

setup(name='rf',
      version=VERSION,
      description=DESCRIPTION,
      long_description=STATUS,
      url='https://github.com/trichter/rf',
      author='Tom Eulenfeld',
      author_email='tom.eulenfeld@gmail.de',
      license='MIT',
      packages=find_packages(),
      package_dir={'rf': 'rf'},
      install_requires=REQUIRES,
      entry_points=ENTRY_POINTS,
      include_package_data=True,
      zip_safe=False,
      classifiers=CLASSIFIERS
      )
