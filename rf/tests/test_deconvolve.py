"""
Tests for deconvolve module.
"""
from numpy.random import random, seed
import numpy as np
import scipy.linalg
import unittest
import rf


class DeconvolveTestCase(unittest.TestCase):

    def setUp(self):
        # set specific seed value such that random numbers are reproducible
        seed(42)
        self.Z = random(412) - 0.5
        self.N = random(412) - 0.5
        self.E = random(412) - 0.5

    def test_toeplitz_real_sym(self):
        src = random(50) - 0.5
        rsp = random(50) - 0.5
        toep = scipy.linalg.toeplitz(src)
        x = np.dot(scipy.linalg.inv(toep), rsp)  # compare to scipy.linalg
        x2 = rf.deconvolve._toeplitz_real_sym(src, rsp)
        np.testing.assert_array_almost_equal(x, x2, decimal=3)

    def test_deconvolution(self):
        #TODO: test_deconvolution
        pass


def suite():
    return unittest.makeSuite(DeconvolveTestCase, 'test')

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
