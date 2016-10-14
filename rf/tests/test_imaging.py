# Copyright 2013-2016 Tom Eulenfeld, MIT license
"""
Tests for imaging module.
"""

import unittest
import warnings

from rf.util import minimal_example_rf

try:
    import cartopy
except ImportError:
    cartopy = None


class ImagingTestCase(unittest.TestCase):

    def test_plot_rf(self):
        stream = minimal_example_rf()
        stream.select(component='Q').plot_rf()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            stream.plot_rf()

    @unittest.skipIf(cartopy is None, 'cartopy not installed')
    def test_plot_ppoints(self):
        from rf.imaging import plot_ppoints
        stream = minimal_example_rf()
        plot_ppoints(stream.ppoints(50), inventory=stream)


def suite():
    return unittest.makeSuite(ImagingTestCase, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
