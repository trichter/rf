# Copyright 2013-2016 Tom Eulenfeld, MIT license
"""
Tests for harmonics module.
"""
import unittest

from rf.rfstream import read_rf, rfstats
from rf.harmonics import harmonics
import warnings

class HarmonicsTestCase(unittest.TestCase):
    def test_harmonics_RT(self):
        stream = read_rf()
        rfstats(stream)
        stream.filter('bandpass',freqmin=0.33,freqmax=1)
        stream.trim2(-5,15)
        stream.rf(deconvolve='multitaper',rotate='NE->RT',\
                    gauss=0.5,K=3,tband=4,T=10,olap=0.75,normalize=0)
        with warnings.catch_warnings():  # warns for sparse data; generally bad but ok for test
            warnings.simplefilter("ignore")
            harm = harmonics(stream,components='RT',scalars=(1,1),method='time')
        trtest = harm.select(location='mod',channel='0')[0]
        self.assertEqual(trtest.data.argmax(),100)
        self.assertEqual(trtest.data.max().round(4),1.1839)

    def test_harmonics_R(self):
        stream = read_rf()
        rfstats(stream)
        stream.filter('bandpass',freqmin=0.33,freqmax=1)
        stream.trim2(-5,15)
        stream.rf(deconvolve='multitaper',rotate='NE->RT',\
                    gauss=0.5,K=3,tband=4,T=10,olap=0.75,normalize=0)
        with warnings.catch_warnings():  # warns for sparse data; generally bad but ok for test
            warnings.simplefilter("ignore")
            harm = harmonics(stream,components='R',method='time')
        trtest = harm.select(location='mod',channel='0')[0]
        self.assertEqual(trtest.data.argmax(),100)
        self.assertEqual(trtest.data.max().round(4),1.7725)

def suite():
    return unittest.makeSuite(HarmonicsTestCase, 'test')

if __name__ == '__main__':
    unittest.main(defaultTest='suite')

