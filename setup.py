# Copyright 2013-2016 Tom Eulenfeld, MIT license
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
LONG_DESCRIPTION = (
    'Please look at the project site for tutorials and information.')

ENTRY_POINTS = {
    'console_scripts': ['rf-runtests = rf.tests:run',
                        'rf = rf.batch:run_cli']}

REQUIRES = ['decorator', 'matplotlib>=2', 'numpy', 'scipy',
            'setuptools', 'obspy>=1.0.3',
            'cartopy', 'geographiclib', 'shapely', 'toeplitz', 'tqdm']

EXTRAS_REQUIRE = {
    'doc': ['sphinx', 'alabaster'],  # and decorator, obspy
    'h5': ['obspyh5>=0.3']}

CLASSIFIERS = [
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Topic :: Scientific/Engineering :: Physics'
    ]


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
      install_requires=REQUIRES,
      extras_require=EXTRAS_REQUIRE,
      entry_points=ENTRY_POINTS,
      include_package_data=True,
      zip_safe=False,
      classifiers=CLASSIFIERS
      )
