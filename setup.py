#!/usr/bin/env python
from setuptools import find_packages, setup
import os.path

with open(os.path.join('rf', '_version.py')) as f:
    VERSION = f.read().split('=')[1].strip().strip("'")
with open('README.rst') as f:
    README = f.read()
if not 'dev' in VERSION:  # get image for correct version from travis-ci
    README = README.replace('branch=master', 'branch=v%s' % VERSION)
DESCRIPTION = README.split('\n')[2]
LONG_DESCRIPTION = '\n'.join(README.split('\n')[5:])

ENTRY_POINTS = {
    'console_scripts': ['rf-runtests = rf.tests.suite:main',
                        'rf = rf.batch:main']}

setup(name='rf',
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      url='https://github.com/trichter/rf',
      author='Tom Eulenfeld',
      author_email='tom.eulenfeld@gmail.de',
      license='MIT',
      packages=find_packages(),
      package_dir={'rf': 'rf'},
      install_requires=['obspy>=0.10', 'toeplitz', 'geographiclib'],
      entry_points=ENTRY_POINTS,
      include_package_data=True,
      zip_safe=False
      )
