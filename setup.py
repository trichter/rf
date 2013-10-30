#!/usr/bin/env python
"""
RF
==
Receiver function calculation in seismology
-------------------------------------------
"""

from setuptools import find_packages
from numpy.distutils.core import Extension, setup

ext = Extension(name='rf._xy',
                sources=['rf/src/psmout.f', 'rf/src/pspier.f',
                         'rf/src/sppier.f'])

setup(name='rf',
      version='0.0.1',
      description='Receiver function calculation in seismology',
      long_description=__doc__,
      url='https://github.com/trichter/rf',
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
