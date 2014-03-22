"""
Submodule fetching all test suites and returning a merged test suite.
"""
import unittest
import test_simple_model
import test_deconvolve
import test_rfstream
import test_batch


def suite():
    return unittest.TestSuite([test_simple_model.suite(),
                               test_deconvolve.suite(),
                               test_rfstream.suite(),
                               test_batch.suite()])


def main():
    unittest.TextTestRunner().run(suite())

if __name__ == '__main__':
    main()
