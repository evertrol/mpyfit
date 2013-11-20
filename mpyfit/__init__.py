import numpy
from mpyfit.mpfit import mpfit

__version__ = "1.0"

EPSFCN = numpy.finfo(numpy.float64).eps


def fit(function, p, args=(), parinfo=None, ftol=1e-10, xtol=1e-10, gtol=1e-10,
        epsfcn=EPSFCN, stepfactor=100, covtol=1e-14, maxiter=200, maxfev=0):
    """Robust non-linear least squares curve fitting

    A wrapper around C mpfit.

    Arguments
    ---------

    - function: the function to minimize.

    - parameters: the parameters to optimize.

    Keyword arguments
    -----------------

    - args: extra arguments pass to the function, given as a tuple.

    - parinfo: list of dictionaries with extra parameter information.
      See below for more details.

    - ftol: relative chi-square convergence criterium.

    - xtol: relative parameter convergence criterium.

    - gtol: Orthogonality convergence criterium.

    - epsfcn: Finite derivative step size. Default: double precision
      limits (numpy.finfo(numpy.float64).eps

    - stepfactor: initial step bound.

    - covtol: range tolerance for covariance calculation.

    - maxiter: maximum number of iterations. If maxiter == 0, then
      basic error checking is done, and parameter errors/covariances
      are estimated based on input parameter values, but no fitting
      iterations are done.

    - maxfev: maximum number of function evaluations, or 0 for no
      limit. Default 0 (no limit)

    Returns
    -------

    A two tuple: the first element contains a copy of the fitted
    parameters (i.e., the parameters given as input keep their origin
    values); the second elements is a dictionary with the following
    items:

    - bestnorm: the final chi-squared value.

    - orignorm: the starting chi-squared value.

    - niter: the number of iterations.

    - nfev: the number of function evaluations.

    - status: the fitting status code. See below for details.

    - npar: the total number of parameters.

    - nfree: the total number of free (fitted) parameters.

    - npegged: the total number of pegged parameters

              (parameters that ran into their boundary constraint).

    - nfunc: number of residuals (equal to the number of data
              points).

    - residuals: final residuals.

    - parerrors: 1 sigma estimated parameter uncertainties.

    - covariances: n by n covariance matrix of the parameters.


    Parameter info dictionary
    -------------------------

    A list containing dictionaries can be passed as keyword argument
    parinfo. Each dictionary contains extra information about the
    corresponding parameter. The keys in the dictionary are all
    optional, and can be the following:

    - fixed : bool

      Whether to freeze the parameter (True), or fit the parameter
      (False). The frozen value is that of the initial value given to
      the `fit()` function.

    - limits : 2-tuple containg floats or None

      Set the lower or upper limit for a variable (first and second
      element of the tuple, respectively). Set a limit to None to
      indicate no limit. For example, `(-10, 10)` limits the parameter
      to the range `[-10, 10]`, while `(-10, None)` limits it to `>= -10`.
      The default is no limits (same as `(None, None)`).

    - step : float

      Stepsize when using the finite difference in computing the
      derivative.

    - relstep : float

      Relative stepsize when using the finite difference in
      computing the derivative.

    - side : int

      Sidedness of finite difference derivative:

      * 0 : one-sided derivative computed automatically
      * 1 : one-sided derivative (f(x+h) - f(x)  )/h
      * -1 : one-sided derivative (f(x)   - f(x-h))/h
      * 2 : two-sided derivative (f(x+h) - f(x-h))/(2*h)
      * 3 : user-computed analytical derivatives

    - deriv_debug : bool

      Derivative debug mode. If True, compute both analytical and
      numerical derivatives and print them to the console for
      comparison. NOTE: when debugging, do *not* set side = 3, but
      rather to the kind of numerical derivative you want to compare
      the user-analytical one to (0, 1, -1, or 2).

    - deriv_reltol : float

      Relative tolerance for derivative debug printout.

    - deriv_abstol : float

      Absolute tolerance for derivative debug printout.


    Example
    -------

    >>> from __future__ import print_function
    >>> import mpyfit
    >>> mpyfit.fit  #doctest: +ELLIPSIS
    <function fit at 0x...>
    >>> import numpy    

    >>> # Define the actual function
    >>> def func(x, p):
    ...     return p[0] + p[1] * numpy.sin(p[2]*x - p[3])
    ... 
    >>> # A simple minimization function:
    >>> def least(p, args):
    ...     x, y = args
    ...     return func(x, p) - y
    ... 
    >>> p = [1, 1.5, 0.2, 0.5]
    >>> x = numpy.linspace(-10, 10, 20)
    >>> y = func(x, p)
    >>> # Add some noise
    >>> y += numpy.random.normal(0, 0.01, y.shape)

    >>> pstart = [1, 1, 0.1, 1]
    >>> pfit, results = mpyfit.fit(least, pstart, (x, y))

    >>> print([round(p, 1) for p in pfit])  # doctest: +NORMALIZE_WHITESPACE
    [1.0, 1.5, 0.2, 0.5]

    """

    if parinfo is None:
        parinfo = [{}] * len(p)
    p, result = mpfit(
        function, numpy.asarray(p), args=args, parinfo=parinfo, ftol=ftol, 
        xtol=xtol, gtol=gtol, epsfcn=epsfcn, stepfactor=stepfactor, 
        covtol=covtol, maxiter=maxiter, maxfev=maxfev)
    status = {
        1: 'Convergence in chi-squared value', 
        2: 'Convergence in parameter value',
        3: 'Convergence in both chi-squared and parameter value',
        4: 'Convergence in orthogonality',
        5: 'Maximum number of iterations reached',
        6: 'ftol is too small: no further improvement',
        7: 'xtol is too small: no further improvement',
        8: 'gtol is too small: no further improvement',
        0: 'General input parameter error',
        -16: 'User function produced non-finite values',
        -17: 'No user function was supplied',
        -18: 'No user data points were supplied',
        -19: 'No free parameters',
        -20: 'Memory allocation error',
        -21: 'Initial values inconsistent with constraints',
        -22: 'Initial constraints inconsistent',
        -23: 'General input parameter error',
        -24: 'Not enough degrees of freedom',
    }
    result['status'] = (result['status'], status[result['status']])
    return p, result
