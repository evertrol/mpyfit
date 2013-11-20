"""
These tests are the test fits that come with C mpfit
"""

import mpyfit
import unittest
import numpy


class LinearFunction( unittest.TestCase):

    @staticmethod
    def func(p, args):
        x, y, error = args
        return (y - p[0] - p[1]*x)/error

    def test_fit(self):
        x = numpy.array([-1.7237128E+00,1.8712276E+00,-9.6608055E-01,
                         -2.8394297E-01,1.3416969E+00,1.3757038E+00,
                         -1.3703436E+00,4.2581975E-02,-1.4970151E-01,
                         8.2065094E-01])
        y = numpy.array([1.9000429E-01,6.5807428E+00,1.4582725E+00,
                         2.7270851E+00,5.5969253E+00,5.6249280E+00,
                         0.787615,3.2599759E+00,2.9771762E+00,
                         4.5936475E+00])
        error = 0.07 * numpy.ones(x.shape)
        p = [1.0, 1.0]
        pactual = [3.20, 1.78]

        args = (x, y, error)
        p, result = mpyfit.fit(self.func, p, args)

        dof = result['nfunc'] - result['nfree']
        self.assertEqual(dof, 8)
        self.assertEqual(result['status'][0], 1)
        self.assertAlmostEqual(result['bestnorm'], 2.756285)
        self.assertTrue(abs(p[0] - pactual[0]) < result['parerrors'][0])
        self.assertTrue(abs(p[1] - pactual[1]) < result['parerrors'][1])


class QuadraticFunction(unittest.TestCase):

    @staticmethod
    def func(p, args):
        x, y, error = args
        return (y - p[0] - p[1]*x - p[2]*x*x)/error

    def setUp(self):
        self.p = numpy.ones((3,))
        self.pactual = [4.7, 0.0, 6.2]
        self.x = numpy.asarray([-1.7237128E+00,1.8712276E+00,-9.6608055E-01,
                                -2.8394297E-01,1.3416969E+00,1.3757038E+00,
                                -1.3703436E+00,4.2581975E-02,-1.4970151E-01,
                                8.2065094E-01])
        self.y = numpy.asarray([2.3095947E+01,2.6449392E+01,1.0204468E+01,
                                5.40507,1.5787588E+01,1.6520903E+01,
                                1.5971818E+01,4.7668524E+00,4.9337711E+00,
                                8.7348375E+00])
        self.error = 0.2 * numpy.ones(self.x.shape)

    def test_fit(self):
        args = (self.x, self.y, self.error)
        p, result = mpyfit.fit(self.func, self.p, args)

        dof = result['nfunc'] - result['nfree']
        self.assertEqual(result['status'][0], 1)
        self.assertEqual(dof, 7)
        self.assertTrue(abs(p[0] - self.pactual[0]) < result['parerrors'][0])
        self.assertFalse(abs(p[1] - self.pactual[1]) < result['parerrors'][1])
        self.assertTrue(abs(p[2] - self.pactual[2]) < result['parerrors'][2])

    def test_frozen(self):
        args = (self.x, self.y, self.error)
        # Fix the second term
        parinfo = [{}, {'fixed': True}, {}]
        self.p[1] = 0
        p, result = mpyfit.fit(self.func, self.p, args=args, parinfo=parinfo)

        dof = result['nfunc'] - result['nfree']
        self.assertEqual(result['status'][0], 1)
        self.assertEqual(dof, 8)
        self.assertTrue(abs(p[0] - self.pactual[0]) < result['parerrors'][0])
        self.assertEqual(p[1], 0.0)
        self.assertTrue(abs(p[2] - self.pactual[2]) < result['parerrors'][2])


class GaussFunction(unittest.TestCase):

    @staticmethod
    def func(p, args):
        x, y, error = args
        xc = x - p[2]
        return (y - p[1] * numpy.exp(-0.5*(x-p[2])*(x-p[2]) / (p[3]*p[3])) -
                p[0]) / error

    def setUp(self):
        self.x = numpy.array([-1.7237128E+00,1.8712276E+00,-9.6608055E-01,
                              -2.8394297E-01,1.3416969E+00,1.3757038E+00,
                              -1.3703436E+00,4.2581975E-02,-1.4970151E-01,
                              8.2065094E-01])
        self.y = numpy.array([-4.4494256E-02,8.7324673E-01,7.4443483E-01,
                              4.7631559E+00,1.7187297E-01,1.1639182E-01,
                              1.5646480E+00,5.2322268E+00,4.2543168E+00,
                              6.2792623E-01])
        self.p = numpy.ones((4,))
        self.p[0] = 0
        self.pactual = [0.0, 4.7, 0.0, 0.5]
        self.error = 0.5 * numpy.ones(self.x.shape)

    def test_fit(self):
        args = (self.x, self.y, self.error)
        p, result = mpyfit.fit(self.func, self.p, args)

        dof = result['nfunc'] - result['nfree']
        self.assertEqual(result['status'][0], 1)
        self.assertEqual(dof, 6)
        self.assertFalse(abs(p[0] - self.pactual[0]) < result['parerrors'][0])
        self.assertTrue(abs(p[1] - self.pactual[1]) < result['parerrors'][1])
        self.assertTrue(abs(p[2] - self.pactual[2]) < result['parerrors'][2])
        self.assertFalse(abs(p[3] - self.pactual[3]) < result['parerrors'][3])

    def test_freeze(self):
        self.p[2] = 0
        self.p[3] = 0.1
        parinfo = [{'fixed': True}, {}, {'fixed': True}, {}]

        args = (self.x, self.y, self.error)
        p, result = mpyfit.fit(self.func, self.p, args, parinfo=parinfo)

        dof = result['nfunc'] - result['nfree']
        self.assertEqual(result['status'][0], 1)
        self.assertEqual(result['nfunc'] - result['nfree'], 8)
        self.assertEqual(p[0], 0.0)
        self.assertTrue(abs(p[1] - self.pactual[1]) < 1.5*result['parerrors'][1])
        self.assertEqual(p[2], 0.0)
        self.assertTrue(abs(p[3] - self.pactual[3]) < result['parerrors'][3])

    def test_limit(self):
        self.p[2] = 0
        self.p[3] = 0.1
        parinfo = [{'fixed': True}, {}, {'fixed': True},
                   {'limits': (-0.2, 0.3)}]

        args = (self.x, self.y, self.error)
        p, result = mpyfit.fit(self.func, self.p, args, parinfo=parinfo)

        dof = result['nfunc'] - result['nfree']
        self.assertEqual(result['status'][0], 1)
        self.assertEqual(result['nfunc'] - result['nfree'], 8)
        self.assertEqual(p[0], 0.0)
        self.assertEqual(result['parerrors'][0], 0)
        self.assertAlmostEqual(p[1], 5.533173)
        self.assertAlmostEqual(result['parerrors'][1], 0.3395397)
        self.assertEqual(p[2], 0.0)
        self.assertEqual(result['parerrors'][2], 0)
        self.assertTrue(p[3] >= -0.2 and p[3] <= 0.3)


if __name__ == '__main__':
    unittest.main()
