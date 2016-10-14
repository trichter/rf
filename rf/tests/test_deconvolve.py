# Copyright 2013-2016 Tom Eulenfeld, MIT license
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
    msg = 'maximum of L component is %.2f, but should 1'
    max_ = Ltrace.data.max()
    testcase.assertGreater(max_, 0.999, msg=msg % max_)


def test_deconvolve_Qpeak(testcase, stream, *args, **kwargs):
    s1 = stream.copy()
    winsrc = kwargs.pop('winsrc', 'S')
    s1.deconvolve(*args, winsrc=winsrc, source_components='Q', **kwargs)
    Qtrace = s1.select(component='Q')[0]
    st = Qtrace.stats
    onsetQ = Qtrace.data.argmax() * st.delta
    onset = st.onset - st.starttime
    msg = ('Q component maxium at %.2fs, but should be at %.2fs. '
           'deconvolve args: %s, %s') % (onsetQ, onset, args, kwargs)
    testcase.assertLess(abs(onsetQ-onset), 0.01, msg=msg)
    msg = 'maximum of Q component is %.2f, but should 1'
    max_ = Qtrace.data.max()
    testcase.assertGreater(max_, 0.999, msg=msg % max_)


class DeconvolveTestCase(unittest.TestCase):

    def test_toeplitz_real_sym(self):
        # set specific seed value such that random numbers are reproducible
        seed(0)
        src = random(50) - 0.5
        rsp = random(50) - 0.5
        toep = scipy.linalg.toeplitz(src)
        x = np.dot(scipy.linalg.inv(toep), rsp)  # compare to scipy.linalg
        x2 = rf.deconvolve._toeplitz_real_sym(src, rsp)
        np.testing.assert_array_almost_equal(x, x2, decimal=3)

    def test_deconvolution_of_stream_Lpeak_position(self):
        stream = read_rf()[:3]
        rfstats(stream)
        stream.filter('bandpass', freqmin=0.4, freqmax=1)
        stream.trim2(5, 95, reftime='starttime')
        stream.rotate('ZNE->LQT')
        # check that maximum in L component is at 0s (at P onset)
        test_deconvolve_Lpeak(self, stream, 'time')
        test_deconvolve_Lpeak(self, stream, 'freq')
        test_deconvolve_Lpeak(self, stream, 'time', winsrc=(-20, 40, 5))
        test_deconvolve_Lpeak(self, stream, 'freq', winsrc=(-20, 40, 5))
        stream.trim2(5, 70, reftime='starttime')
        test_deconvolve_Lpeak(self, stream, 'time')
        test_deconvolve_Lpeak(self, stream, 'freq')
        test_deconvolve_Lpeak(self, stream, 'time', winsrc=(-20, 40, 5))
        test_deconvolve_Lpeak(self, stream, 'freq', winsrc=(-20, 40, 5))
        test_deconvolve_Lpeak(self, stream, 'time', winsrc=(-5, 18, 5))
        test_deconvolve_Lpeak(self, stream, 'freq', winsrc=(-5, 18, 5))

    def test_deconvolution_of_stream_Qpeak_position(self):
        # S receiver deconvolution
        from pkg_resources import resource_filename
        fname = resource_filename('rf', 'example/minimal_example_S.tar.gz')
        stream = read_rf(fname)[:3]
        rfstats(stream, phase='S')
        stream.filter('bandpass', freqmin=0.2, freqmax=0.5)
        stream.trim2(10, 120, reftime='starttime')
        stream.rotate('ZNE->LQT')
        # check that maximum in Q component is at 0s (at S onset)
        test_deconvolve_Qpeak(self, stream, 'time')
        test_deconvolve_Qpeak(self, stream, 'freq')
        test_deconvolve_Qpeak(self, stream, 'time', winsrc=(-5, 18, 5))
        test_deconvolve_Qpeak(self, stream, 'freq', winsrc=(-5, 18, 5))
        test_deconvolve_Qpeak(self, stream, 'time', winsrc=(-20, 40, 5))
        test_deconvolve_Qpeak(self, stream, 'freq', winsrc=(-20, 40, 5))

    def test_deconvolution_of_convolution(self):
        from rf.rfstream import RFStream, RFTrace
        data = np.zeros(400)
        data_src = np.zeros(400)
        hann1 = get_window('hann', 10)
        hann2 = get_window('hann', 50)
        data_src[40:50] = hann1
        data_src[50:60] = -hann1
        data[100:150] = hann2
        data[240:290] = 0.5 * hann2
        data_rsp = convolve(data_src, data, 'full')[50:450]/3.
        stream1 = RFStream([RFTrace(data=data_src), RFTrace(data=data_rsp)])
        for i, tr in enumerate(stream1):
            tr.stats.channel = tr.stats.channel[:2] + 'LQT'[i]
            tr.stats.onset = tr.stats.starttime + 40
        stream2 = stream1.copy()
        stream1.deconvolve(spiking=10)
        stream2.deconvolve(method='freq', waterlevel=0.1)

#        import matplotlib.pyplot as plt
#        plt.subplot(121)
#        plt.plot(data, label='desired')
#        plt.plot(data_src, label='source')
#        plt.plot(data_rsp, label='convolution')
#        plt.plot(stream1[0].data, label='deconv src')
#        plt.plot(stream1[1].data, label='deconv')
#        plt.legend()
#        plt.subplot(122)
#        plt.plot(data)
#        plt.plot(data_src)
#        plt.plot(data_rsp)
#        plt.plot(stream2[0].data)
#        plt.plot(stream2[1].data)
#        plt.show()

        # (shift from middle of source (50) to onset (40)
        peakpos = np.argmax(data) - 10
        self.assertEqual(peakpos, np.argmax(stream1[1].data))
        self.assertEqual(peakpos, np.argmax(stream2[1].data))


def suite():
    return unittest.makeSuite(DeconvolveTestCase, 'test')

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
