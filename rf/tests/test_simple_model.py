"""
Tests for simple_model module.
"""
import unittest
from rf.simple_model import load_model
from obspy import read
from rf import RFStream


class SimpleModelTestCase(unittest.TestCase):

    def setUp(self):
        self.stream = RFStream(read())
        self.stream._write_test_header()
        self.model = load_model()

    def test_ppoint(self):
        xy_plat_P = 50.29670878  # 200km sppier
        xy_plon_P = -97.24776461  # 200km sppier
        xy_plat_S = 50.29862202  # 200km pspier
        xy_plon_S = -98.96394212  # 200km pspier

        st = self.stream[0].stats
        self.model.ppoint(st, 200, phase='P')
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

        self.model.ppoint(st, 200, phase='S')
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
        self.model.moveout(self.stream)
        i = 20
        for tr in self.stream:
            tr.stats.distance = i
            i = i + 50
        self.model.moveout(self.stream)


def suite():
    return unittest.makeSuite(SimpleModelTestCase, 'test')

if __name__ == '__main__':
    unittest.main(defaultTest='suite')
