# Copyright 2013-2016 Tom Eulenfeld, MIT license
"""
Tests for batch module.
"""
from glob import glob
import unittest
import os
from pkg_resources import load_entry_point
import sys
import warnings

import matplotlib
matplotlib.use('Agg')

from rf.batch import init_data, run_cli as script
from rf.tests.util import quiet, tempdir
try:
    import obspyh5
except ImportError:
    obspyh5 = None


def substitute(old, new):
    fname = 'conf.json'
    with open(fname) as f:
        text = f.read()
    text = text.replace(old, new)
    with open(fname, 'w') as f:
        f.write(text)


def test_format(testcase, format):
    join = os.path.join

    # QHD
    with tempdir():
        script(['create', '-t'])
        substitute('#"format": "Q"', '"format": "%s"' % format)
        script(['data', 'data'])
        script(['calc', 'moveout', 'data', 'mout1'])
        script(['data', 'calc', 'datarf'])
        script(['moveout', 'datarf', 'mout2'])
        script(['stack', 'mout1', 'stack'])
        script(['profile', 'mout1', 'profile'])
        if format in ('Q', 'SAC'):
            patterns = [join('data', '*', '*'), join('mout1', '*', '*'),
                        join('mout2', '*', '*'), join('stack', '*'),
                        join('profile', '*') if format == 'SAC' else
                        'profile*.*']
        else:
            patterns = ['data.h5', 'mout1.h5', 'mout2.h5', 'stack.h5',
                        'profile.h5']
        nums = [len(glob(p)) for p in patterns]
        nums2 = {'Q': [14, 14, 14, 2, 2],
                 'SAC': [21, 21, 21, 3, 6],
                 'H5': [1, 1, 1, 1, 1]}
        testcase.assertEqual(nums, nums2[format])
        if format in ('Q', 'H5'):
            script(['convert', 'mout1', 'mout_SAC', 'SAC'])
        if format in ('H5', 'SAC'):
            script(['convert', 'mout1', 'mout_Q', 'Q'])
        if obspyh5 and format in ('Q', 'SAC'):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                script(['convert', 'mout1', 'mout_H5', 'H5'])
        script(['plot', 'mout1', 'plot'])
        script(['plot-profile', 'profile', 'plot_profile'])
        testcase.assertEqual(len(glob(join('plot', '*.pdf'))), 3)
        testcase.assertEqual(len(glob(join('plot_profile', '*.pdf'))), 3)


class BatchTestCase(unittest.TestCase):
    """
    For now batch tests are not run on windows.
    See https://travis-ci.org/trichter/rf/jobs/589375488 for failures.
    The batch module has a low priority.
    """
    def setUp(self):
        # turn off progressbar
        import rf.batch
        rf.batch.tqdm = lambda: None

    @unittest.skipIf(sys.platform.startswith("win"), "fails on Windows")
    def test_entry_point(self):
        ep_script = load_entry_point('rf', 'console_scripts', 'rf')
        try:
            with quiet():
                ep_script(['-h'])
        except SystemExit:
            pass

    @unittest.skipIf(sys.platform.startswith("win"), "fails on Windows")
    def test_batch_command_interface_Q(self):
        test_format(self, 'Q')

    @unittest.skipIf(sys.platform.startswith("win"), "fails on Windows")
    def test_batch_command_interface_SAC(self):
        test_format(self, 'SAC')

    @unittest.skipIf(obspyh5 is None, 'obspyh5 not installed')
    @unittest.skipIf(sys.platform.startswith("win"), "fails on Windows")
    def test_batch_command_interface_H5(self):
        test_format(self, 'H5')

    def test_plugin_option(self):
        f = init_data('plugin', plugin='rf.tests.test_batch : gw_test')
        self.assertEqual(f(nework=4, station=2), 42)


def gw_test(**kwargs):
    return 42


def suite():
    return unittest.makeSuite(BatchTestCase, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
