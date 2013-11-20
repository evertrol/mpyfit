"""A number of hard-to-fit functions taken from Wikipedia:
http://en.wikipedia.org/wiki/Test_functions_for_optimization

The number of iterations for convergence may differ, so if these tests
fail, carefully check if this is simply because the test fit has not
converged yet. (A.k.a. tests need to be improved.)

"""

import mpyfit
import unittest
import numpy


class SphereFunction(unittest.TestCase):

    @staticmethod
    def func(p):
        return p[0]**2 + p[1]**2

    def test_fit(self):
        p = [1, 1]
        p, result = mpyfit.fit(self.func, p)
        self.assertEqual(result['status'][0], 5)
        self.assertAlmostEqual(p[0], 0)
        self.assertAlmostEqual(p[1], 0)

    def test_fit10000iterations(self):
        p = [1, 1]
        p, result = mpyfit.fit(self.func, p, maxiter=50000)
        self.assertEqual(result['status'][0], 4)
        self.assertAlmostEqual(p[0], 0)
        self.assertAlmostEqual(p[1], 0)


class RosenbrockFunction(unittest.TestCase):

    @staticmethod
    def func(p):
        return (1 - p[0])**2 + 100 * (p[1] - p[0]*p[0])**2

    def test_fit(self):
        p = [0, 0]
        p, result = mpyfit.fit(self.func, p)
        # Does not even come close
        self.assertNotAlmostEqual(round(p[0], 1), 1.0)
        self.assertNotAlmostEqual(round(p[1], 1), 1.0)

    def test_fit50000iterations(self):
        p = [0, 0]
        p, result = mpyfit.fit(self.func, p, maxiter=50000)
        self.assertAlmostEqual(p[0], 1., places=4)
        self.assertAlmostEqual(p[1], 1., places=4)


class BealesFunction(unittest.TestCase):

    @staticmethod
    def func(p):
        return ((1.5 - p[0] + p[0]*p[1])**2 + 
                (2.25 - p[0] + p[0]*p[1]*p[1])**2 +
                (2.625 - p[0] + p[0]*p[1]*p[1]*p[1])**2)

    def test_fit(self):
        p = [0, 0]
        p, result = mpyfit.fit(self.func, p, maxiter=10000)
        self.assertEqual(result['status'][0], 2)
        self.assertAlmostEqual(p[0], 3.0, places=5)
        self.assertAlmostEqual(p[1], 0.5, places=5)



if __name__ == '__main__':
    unittest.main()
