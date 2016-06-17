"""
Tests for profile module.
"""
import unittest

import numpy as np
from rf.profile import get_profile_boxes
from rf.tests.test_rfstream import test_io_header
from rf.util import minimal_example_rf


class ProfileTestCase(unittest.TestCase):

    def test_profile(self):
        dx = np.linspace(-50, 50, 41)
        boxes = get_profile_boxes((-21, -69.5), 85, dx, width=300)
        stream = minimal_example_rf()
        stream.extend(stream[:3])  # to actually stack something
        profile = stream.select(component='Q').profile(boxes)
        self.assertEqual(len(profile), 3)
        str_ = 'profile | -5.0s - 22.0s | -6.25km slow: 6.40 (Ps moveout)'
        self.assertEqual(str(profile[0]), str_)
        test_io_header(self, profile[:1])


def suite():
    return unittest.makeSuite(ProfileTestCase, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
