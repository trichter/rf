"""
Tests for deconvolve module.
"""
from numpy.random import random, seed
import numpy as np
import scipy.linalg
from scipy.signal import get_window, convolve
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
        ms = rf.read_rf()
        ms.decimate(10)
        for i in range(len(ms)):
            ms[i].stats.channel = ms[i].stats.channel[:2] + 'LQT'[i]
        t = np.linspace(0, 30, len(ms[0]))
        hann1 = get_window('hann', 10)
        hann2 = get_window('hann', 50)
        ms[0].data[:] = 0
        ms[0].data[40:50] = hann1
        ms[0].data[50:60] = -hann1

        ms[1].data[:] = 0
        ms[1].data[100:150] = hann2
        ms[1].data[240:290] = hann2
        ms_orig = ms.copy()
        data3 = convolve(ms[1].data, ms[0].data, 'full')[50:350]/np.sum(np.abs(ms[0].data))
        ms[1].data = data3
        ms[2].data = -ms[1].data
        for tr in ms:
            tr.stats.sampling_rate = 1
            tr.stats.onset = tr.stats.starttime + 40
        ms.deconvolve(method='time')
        #ms.deconvolve(method='freq', tshift=45.)
#        print(ms)

#        import pylab as plt
#        plt.plot(t,  ms[0].data)
#        plt.plot(t, ms_orig[0].data)
#        plt.plot(t, data3)
#        plt.plot(t, ms[1].data)
#        plt.plot(t, ms_orig[1].data)

#        plt.plot(ms[0].data)
#        plt.plot(ms_orig[0].data)
#        plt.plot(ms[1].data)
#        plt.plot(ms_orig[1].data)

#        ms.plot()
#        plt.show()


#        print(ms)


        #TODO: test_deconvolution


def suite():
    return unittest.makeSuite(DeconvolveTestCase, 'test')

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
