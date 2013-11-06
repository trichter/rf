"""
Tests for simple_model module.
"""
import unittest
from rf.simple_model import moveout, ppoint
from obspy import read
from rf import RFStream


class SimpleModelTestCase(unittest.TestCase):

    def setUp(self):
        self.stream = RFStream(read())
        self.stream._write_test_header()

    def test_ppoint(self):
        xy_plat_P = 50.29670878  # 200km sppier
        xy_plon_P = -97.24776461  #200km sppier
        xy_plat_S = 50.29862202  # 200km pspier
        xy_plon_S = -98.96394212  #200km pspier

        st = self.stream[0].stats
        ppoint(self.stream, 200, phase='P')
        self.assertLess(abs((st.plat - xy_plat_P) /
                            (st.plat - st.station_latitude)), 1)  # very big
        self.assertLess(abs((st.plon - xy_plon_P) /
                            (st.plon - st.station_longitude)), 0.05)

#        print 'station lon ', st.station_longitude
#        print 'xy lon    P ', xy_plon_P
#        print 'model lon P ', st.plon
#        print
#        print 'station lat ', st.station_latitude
#        print 'xy lat    P ', xy_plat_P
#        print 'model lat P ', st.plat
#        print

        ppoint(self.stream, 200, phase='S')
        self.assertLess(abs((st.plat - xy_plat_S) /
                            (st.plat - st.station_latitude)), 1)  # very big
        self.assertLess(abs((st.plon - xy_plon_S) /
                            (st.plon - st.station_longitude)), 0.05)

#        print 'station lon ', st.station_longitude
#        print 'xy lon    S ', xy_plon_S
#        print 'model lon S ', st.plon
#        print
#        print 'station lat ', st.station_latitude
#        print 'xy lat    S ', xy_plat_S
#        print 'model lat S ', st.plat

    def test_moveout(self):
#        st1 = self.stream.copy()
#        moveout(self.stream)
#        st2 = self.stream
#        st = st1[:1] + st2[:1]
#        st = st1 + st2
#        st.plot(automerge=False)
        moveout(self.stream)



def suite():
    return unittest.makeSuite(SimpleModelTestCase, 'test')

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
