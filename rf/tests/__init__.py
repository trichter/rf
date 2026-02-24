"""
Tests for the rf package.
"""

import importlib.resources as imp_resources
import sys
import unittest

import matplotlib
matplotlib.use('Agg')


def run():
    loader = unittest.TestLoader()
    test_dir = imp_resources.files('rf') / 'tests'
    suite = loader.discover(test_dir)
    runner = unittest.runner.TextTestRunner()
    ret = not runner.run(suite).wasSuccessful()
    sys.exit(ret)
