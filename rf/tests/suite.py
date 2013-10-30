import unittest
#from test_deconvolve import suite as suite1
#from test_rfstream import suite as suite2
import test_deconvolve
import test_rfstream


def suite():
    return unittest.TestSuite([test_deconvolve.suite(),
                               test_rfstream.suite()])

def main():
    unittest.TextTestRunner().run(suite())

if __name__ == '__main__':
    main()