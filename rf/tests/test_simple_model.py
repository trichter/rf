"""
Tests for simple_model module.
"""
import numpy as np
from obspy import read
from obspy.core.util.geodetics import degrees2kilometers
from obspy.taup import TauPyModel
from rf import RFStream
from rf.simple_model import load_model
import unittest


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

    def test_ppointvsobspytaup_S2P(self):
        slowness = 12.33
        evdep = 12.4
        evdist = 67.7
        pp1 = self.model.ppoint_distance(200, slowness, phase='P')
        model = TauPyModel(model='iasp91')
        arrivals = model.get_ray_paths(evdep, evdist, ('S250p',))
        arrival = arrivals[0]
        index = np.searchsorted(arrival.path['depth'][::-1], 200)
        pdist = arrival.path['dist']
        pp2 = degrees2kilometers((pdist[-1] - pdist[-index-1]) * 180 / np.pi)
        self.assertLess(abs(pp1-pp2)/pp2, 0.2)

    def test_ppointvsobspytaup_P2S(self):
        slowness = 6.28
        evdep = 12.4
        evdist = 67.7
        depth = 200
        pp1 = self.model.ppoint_distance(depth, slowness)
        model = TauPyModel(model='iasp91')
        arrivals = model.get_ray_paths(evdep, evdist, ('P250s',))
        arrival = arrivals[0]
        index = np.searchsorted(arrival.path['depth'][::-1], depth)
        pdist = arrival.path['dist']
        pp2 = degrees2kilometers((pdist[-1] - pdist[-index-1]) * 180 / np.pi)
        self.assertLess(abs(pp1-pp2)/pp2, 0.1)

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
