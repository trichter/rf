#!/usr/bin/env python
from setuptools import find_packages, setup
import os.path

with open(os.path.join('rf', '_version.py')) as f:
    VERSION = f.read().split('=')[1].strip().strip("'")
with open('README.rst') as f:
    README = f.read()
if not 'dev' in VERSION:  # get image for correct version from travis-ci
    README = README.replace('branch=master', 'branch=v%s' % VERSION)
readme = README.split('\n')
DESCRIPTION = readme[2]
LONG_DESCRIPTION = '\n'.join(readme[5:7] + readme[9:10] + readme[12:])

ENTRY_POINTS = {
    'console_scripts': ['rf-runtests = rf.tests:run',
                        'rf = rf.batch:run_cli']}

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
