"""Test that verify some edge cases."""

import mpyfit
import unittest
import numpy


class SingleParameter(unittest.TestCase):

    @staticmethod
    def func(p):
        return p[0]

    def test_fit(self):
        p = [10]
        p, result = mpyfit.fit(self.func, p)
        self.assertLess(result['niter'], 10)
        self.assertEqual(result['status'][0], 4)
        self.assertAlmostEqual(p[0], 0)

    # Python 2/3 compatibility
    def assertRaisesRegex23(self, error, regex):
        try:
            return self.assertRaisesRegex(error, regex)  # Python 3
        except AttributeError:
            return self.assertRaisesRegexp(error, regex)  # Python 2

    def test_fixed(self):
        p = [10]
        parinfo = [{'fixed': True}]
        with self.assertRaisesRegex23(RuntimeError, "mpfit function error -19"):
            p, result = mpyfit.fit(self.func, p, parinfo=parinfo)

    def test_limits(self):
        p = [10]
        parinfo = [{'limits': (10, 10)}]
        with self.assertRaisesRegex23(RuntimeError, "mpfit function error -22"):
            p, result = mpyfit.fit(self.func, p, parinfo=parinfo)


if __name__ == '__main__':
    unittest.main()
