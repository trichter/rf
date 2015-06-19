"""
Tests for batch module.
"""
import unittest
import os
from pkg_resources import load_entry_point
import shutil
from subprocess import check_call

import warnings

from rf.batch import _no_pbar, init_data, run_cli as script
from rf.tests.util import quiet, tempdir
try:
    import obspyh5
except ImportError:
    obspyh5 = None

_no_pbar()


def _call(command):
    print command
    check_call(command, shell=False)


class BatchTestCase(unittest.TestCase):

    def test_entry_point(self):
        ep_script = load_entry_point('rf', 'console_scripts', 'rf')
        try:
            with quiet():
                ep_script(['-h'])
        except SystemExit:
            pass

    def test_batch_command_interface(self):
        travis = os.environ.get('TRAVIS') == 'true'
        def substitute(old, new):
            fname = 'conf.json'
            with open(fname) as f:
                text = f.read()
            text = text.replace(old, new)
            with open(fname, 'w') as f:
                f.write(text)

        # QHD
        with quiet(), tempdir():
            script(['create'])
            script(['calc'])
            script(['moveout'])
            script(['convert', 'Prf_Ps', 'SAC'])
            if obspyh5:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    script(['convert', 'Prf_Ps', 'H5'])
            script(['stack', 'Prf_Ps'])
            if not travis:
                script(['plot', 'Prf_Ps'])

        # SAC
        with quiet(), tempdir():
            script(['create'])
            substitute('"format": "Q"', '"format": "SAC"')
            script(['calc'])
            script(['moveout'])
            script(['convert', 'Prf_Ps', 'Q'])
            if obspyh5:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    script(['convert', 'Prf_Ps', 'H5'])
            script(['stack', 'Prf_Ps'])
            if not travis:
                script(['plot', 'Prf'])
        # H5
        if obspyh5:
            with quiet(), tempdir(), warnings.catch_warnings():
                warnings.simplefilter("ignore")
                script(['create'])
                substitute('"format": "Q"', '"format": "H5"')
                script(['calc'])
                script(['moveout'])
                script(['convert', 'Prf_Ps', 'Q'])
                script(['convert', 'Prf_Ps', 'SAC'])
                script(['stack', 'Prf_Ps'])
                if not travis:
                    script(['plot', 'Prf'])

    def test_plugin_option(self):
        f = init_data('plugin', plugin='rf.tests.test_batch : gw_test')
        self.assertEqual(f(nework=4, station=2), 42)

def gw_test(**kwargs):
    return 42

def suite():
    return unittest.makeSuite(BatchTestCase, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
