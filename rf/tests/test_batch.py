"""
Tests for batch module.
"""
import unittest
import os
import shutil
import tempfile
from subprocess import check_call
from pkg_resources import load_entry_point
from rf.batch import _no_pbar, main as script
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
            ep_script(['-h'])
        except SystemExit:
            pass

    def test_batch_command_interface(self):
        travis = os.environ.get('TRAVIS') == 'true'
        def substitute(old, new):
            fname = os.path.join(temp_path, 'conf.py')
            with open(fname) as f:
                text = f.read()
            text = text.replace(old, new)
            with open(fname, 'w') as f:
                f.write(text)

        # QHD
        temp_path = os.path.join(tempfile.gettempdir(), 'RF_test')
        if os.path.exists(temp_path):
            shutil.rmtree(temp_path)
        script(['init', '-t', temp_path])
        os.chdir(temp_path)
        script(['calc', 'P'])
        script(['moveout', 'Prf', 'Ps'])
        script(['convert', 'Prf_Ps', 'SAC'])
        if obspyh5:
            script(['convert', 'Prf_Ps', 'H5'])
        script(['stack', 'Prf_Ps'])
        if not travis:
            script(['plot', 'Prf_Ps'])

        # SAC
        if os.path.exists(temp_path):
            shutil.rmtree(temp_path)
        script(['init', '-t', temp_path])
        substitute("format = 'Q'", "format = 'SAC'")
        os.chdir(temp_path)
        script(['calc', 'P'])
        script(['moveout', 'Prf', 'Ps'])
        script(['convert', 'Prf_Ps', 'Q'])
        if obspyh5:
            script(['convert', 'Prf_Ps', 'H5'])
        script(['stack', 'Prf_Ps'])
        if not travis:
            script(['plot', 'Prf'])

        # H5
        if obspyh5:
            if os.path.exists(temp_path):
                shutil.rmtree(temp_path)
            script(['init', '-t', temp_path])
            substitute("format = 'Q'", "format = 'H5'")
            os.chdir(temp_path)
            script(['calc', 'P'])
            script(['moveout', 'Prf', 'Ps'])
            script(['convert', 'Prf_Ps', 'Q'])
            script(['convert', 'Prf_Ps', 'SAC'])
            script(['stack', 'Prf_Ps'])
            if not travis:
                script(['plot', 'Prf'])


def suite():
    return unittest.makeSuite(BatchTestCase, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
