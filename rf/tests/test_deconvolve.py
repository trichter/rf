"""
Tests for deconvolve module.
"""
import unittest

from numpy.random import random, seed
import numpy as np
import scipy.linalg
from scipy.signal import get_window, convolve
import rf.deconvolve
from rf.rfstream import read_rf, rfstats


def test_deconvolve_Lpeak(testcase, stream, *args, **kwargs):
    s1 = stream.copy()
    s1.deconvolve(*args, **kwargs)
    Ltrace = s1.select(component='L')[0]
    st = Ltrace.stats
    onsetL = Ltrace.data.argmax() * st.delta
    onset = st.onset - st.starttime
    msg = ('L component maxium at %.2fs, but should be at %.2fs. '
           'deconvolve args: %s, %s') % (onsetL, onset, args, kwargs)
    testcase.assertLess(abs(onsetL-onset), 0.1, msg=msg)


def test_deconvolve_Qpeak(testcase, stream, *args, **kwargs):
    s1 = stream.copy()
    s1.deconvolve(*args, winsrc='S', source_components='Q', **kwargs)
    Qtrace = s1.select(component='Q')[0]
    st = Qtrace.stats
    onsetQ = Qtrace.data.argmax() * st.delta
    onset = st.onset - st.starttime
    msg = ('Q component maxium at %.2fs, but should be at %.2fs. '
           'deconvolve args: %s, %s') % (onsetQ, onset, args, kwargs)
    testcase.assertLess(abs(onsetQ-onset), 0.2, msg=msg)


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

    def test_deconvolution_of_stream_Lpeak_position(self):
        stream = read_rf()[:3]
        rfstats(stream=stream)
        stream.filter('bandpass', freqmin=0.4, freqmax=1)
        stream.trim2(5, 95, reftime='starttime')
        stream.rotate('ZNE->LQT')
        # check that maximum in L component is at 0s (at P onset)
        test_deconvolve_Lpeak(self, stream, 'time')
        test_deconvolve_Lpeak(self, stream, 'freq')
        test_deconvolve_Lpeak(self, stream, 'time', tshift=15)
        test_deconvolve_Lpeak(self, stream, 'freq', tshift=15)
        stream.trim2(5, 70, reftime='starttime')
        test_deconvolve_Lpeak(self, stream, 'time')
        test_deconvolve_Lpeak(self, stream, 'freq')
        test_deconvolve_Lpeak(self, stream, 'time', tshift=0)
        test_deconvolve_Lpeak(self, stream, 'freq', tshift=0)

    def test_deconvolution_of_stream_Qpeak_position(self):
        # S receiver deconvolution
        from pkg_resources import resource_filename
        fname = resource_filename('rf', 'example/minimal_example_S.tar.gz')
        stream = read_rf(fname)[:3]
        rfstats(stream=stream, phase='S')
        stream.filter('bandpass', freqmin=0.2, freqmax=0.5)
        stream.trim2(10, 120, reftime='starttime')
        stream.rotate('ZNE->LQT')
        # check that maximum in Q component is at 0s (at S onset)
        test_deconvolve_Qpeak(self, stream, 'time')
        test_deconvolve_Qpeak(self, stream, 'freq')
        test_deconvolve_Qpeak(self, stream, 'time', tshift=15)
        test_deconvolve_Qpeak(self, stream, 'freq', tshift=15)

    def test_deconvolution_of_convolution(self):
        # TODO: test_deconvolution_of_convolution
        from obspy import read
        from rf import RFStream
        ms = RFStream(read())
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
        # ms.deconvolve(method='freq', tshift=45.)
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


def suite():
    return unittest.makeSuite(DeconvolveTestCase, 'test')

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
