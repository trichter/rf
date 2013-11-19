#!/usr/bin/env python
"""
rf
==
Receiver function calculation in seismology
-------------------------------------------

For details see the documentation_.

.. _documentation: http://rf.readthedocs.org/
"""

from setuptools import find_packages, setup
import os.path


with open(os.path.join('rf', '_version.py')) as f:
    VERSION = f.read().split('=')[1].strip().strip("'")

setup(name='rf',
      version=VERSION,
      description='Receiver function calculation in seismology',
      long_description=__doc__,
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
