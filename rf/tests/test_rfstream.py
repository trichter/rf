"""
Tests for rfstream module.
"""
import os.path
import tempfile
import unittest

from obspy import read, readEvents
from obspy.core import AttribDict
from rf import RFStream, rfstats
from rf.rfstream import obj2stats
from rf.rfstream import HEADERS, FORMATHEADERS, STATION_GETTER, EVENT_GETTER


class RFStreamTestCase(unittest.TestCase):

    def setUp(self):
        self.event = readEvents()[0]
        self.station = AttribDict({'latitude': 41.818 - 66.7,
                                   'longitude': 79.689,
                                   'elevation': 365.4})
        self.temp = os.path.join(tempfile.gettempdir(), 'RF_TEST')

    def test_io_header(self):
        def test_io_format(format):
            stream1 = stream.copy()
            fname = self.temp + '_IO_FORMAT.' + format.upper()
            if format == 'sh':
                format = 'q'
                fname = self.temp + '.QHD'
            stream1.write(fname, format.upper())
            stream2 = RFStream(stream=read(fname))
            st1 = stream1[0].stats
            st2 = stream2[0].stats
            for head in HEADERS:
                self.assertAlmostEqual(st1[head], st2[head], 4, msg=head)
        stream = RFStream(stream=read())[:1]
        stream._write_test_header()
        for format in FORMATHEADERS:
            test_io_format(format)

    def test_obj2stats(self):
        stats = obj2stats(event=self.event, station=self.station)
        for head, _ in STATION_GETTER + EVENT_GETTER:
            self.assertIn(head, stats)

    def test_rfstats(self):
        stats = rfstats(station=self.station, event=self.event)
        for head in HEADERS:
            self.assertIn(head, stats)
        # event is exactly north from station and around 66.7 degrees away
        self.assertTrue(abs(stats.distance - 66.7) < 1.)
        self.assertTrue(abs(stats.back_azimuth % 360.) < 0.1)
        self.assertTrue(abs(stats.slowness - 6.4) < 0.1)

    def test_simple(self):
        stream = RFStream(stream=read())
        stream._write_test_header()
        stream.rf()
        stream.moveout()
        stream.ppoint(50)

    def test_rf(self):
        pass

    def test_moveout(self):
        pass

    def test_piercing_points(self):
        pass


def suite():
    return unittest.makeSuite(RFStreamTestCase, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
