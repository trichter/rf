#!/usr/bin/env python
from setuptools import find_packages, setup
import os.path

with open(os.path.join('rf', '_version.py')) as f:
    VERSION = f.read().split('=')[1].strip().strip("'")
with open('README.rst') as f:
    README = f.read()
if not 'dev' in VERSION: # get image for correct version from travis-ci
    README = README.replace('branch=master', 'branch=v%s' % VERSION)


setup(name='rf',
      version=VERSION,
      description='Receiver function calculation in seismology',
      long_description=README,
      url='https://github.com/trichter/rf',
      author='Tom Richter',
      author_email='richter@gfz-potsdam.de',
      license='MIT',
      packages=find_packages(),
      package_dir={'rf': 'rf'},
      requires=['obspy', 'toeplitz', 'geographiclib'],
      entry_points={'console_scripts':
                    ['rf-runtests = rf.tests.suite:main',
                     'rf = rf.batch:main']},
      include_package_data=True,
      zip_safe=False
      )
