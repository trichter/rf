"""
Submodule fetching all test suites and returning a merged test suite.
"""
import unittest
import test_deconvolve
import test_rfstream


def suite():
    return unittest.TestSuite([test_deconvolve.suite(),
                               test_rfstream.suite()])


def main():
    unittest.TextTestRunner().run(suite())

if __name__ == '__main__':
    main()
