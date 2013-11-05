#!/usr/bin/env python
"""
rf
==
Receiver function calculation in seismology
-------------------------------------------

For details see the documentation_.

.. _documentation: http://rf.readthedocs.org/
"""

from setuptools import find_packages
from numpy.distutils.core import Extension, setup
import os.path

ext = Extension(name='rf._xy',
                sources=['rf/src/psmout.f', 'rf/src/pspier.f',
                         'rf/src/sppier.f'])

with open(os.path.join('rf', '_version.py')) as f:
    VERSION = f.read().split('=')[1].strip().strip("'")
                         
setup(name='rf',
      version=VERSION,
      description='Receiver function calculation in seismology',
      long_description=__doc__,
      url='https://github.com/trichter/rf',
      download_url='https://github.com/trichter/rf/archive/master.zip',
      author='Tom Richter',
      author_email='richter@gfz-potsdam.de',
      license='MIT',
      packages=find_packages(),
      requires=['obspy', 'toeplitz'],
      ext_modules=[ext],
      entry_points={'console_scripts':
                    ['rf-runtests = rf.tests.suite:main']},
      package_data={'rf': ['data/*.dat']},
      zip_safe=False
      )
