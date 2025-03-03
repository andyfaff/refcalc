import numpy
from numpy import (
    atleast_1d,
    eye,
    mgrid,
    argmin,
    zeros,
    shape,
    squeeze,
    vectorize,
    asarray,
    sqrt,
    Inf,
    asfarray,
    isinf,
)


class Result(dict):
    """Represents the optimization result.

    Attributes
    ----------
    x : ndarray
        The solution of the optimization.
    success : bool
        Whether or not the optimizer exited successfully.
    status : int
        Termination status of the optimizer. Its value depends on the
        underlying solver. Refer to `message` for details.
    message : str
        Description of the cause of the termination.
    fun, jac, hess : ndarray
        Values of objective function, Jacobian and Hessian (if available).
    nfev, njev, nhev : int
        Number of evaluations of the objective functions and of its
        Jacobian and Hessian.
    nit : int
        Number of iterations performed by the optimizer.
    maxcv : float
        The maximum constraint violation.

    Notes
    -----
    There may be additional attributes not listed above depending of the
    specific solver. Since this class is essentially a subclass of dict
    with attribute accessors, one can see which attributes are available
    using the `keys()` method.
    """

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __repr__(self):
        if self.keys():
            m = max(map(len, list(self.keys()))) + 1
            return "\n".join([k.rjust(m) + ": " + repr(v) for k, v in self.items()])
        else:
            return self.__class__.__name__ + "()"


def fminbound(func, x1, x2, args=(), xtol=1e-5, maxfun=500, full_output=0, disp=1):
    """Bounded minimization for scalar functions.

    Parameters
    ----------
    func : callable f(x,*args)
        Objective function to be minimized (must accept and return scalars).
    x1, x2 : float or array scalar
        The optimization bounds.
    args : tuple, optional
        Extra arguments passed to function.
    xtol : float, optional
        The convergence tolerance.
    maxfun : int, optional
        Maximum number of function evaluations allowed.
    full_output : bool, optional
        If True, return optional outputs.
    disp : int, optional
        If non-zero, print messages.
            0 : no message printing.
            1 : non-convergence notification messages only.
            2 : print a message on convergence too.
            3 : print iteration results.


    Returns
    -------
    xopt : ndarray
        Parameters (over given interval) which minimize the
        objective function.
    fval : number
        The function value at the minimum point.
    ierr : int
        An error flag (0 if converged, 1 if maximum number of
        function calls reached).
    numfunc : int
      The number of function calls made.

    See also
    --------
    minimize_scalar: Interface to minimization algorithms for scalar
        univariate functions. See the 'Bounded' `method` in particular.

    Notes
    -----
    Finds a local minimizer of the scalar function `func` in the
    interval x1 < xopt < x2 using Brent's method.  (See `brent`
    for auto-bracketing).

    """
    options = {"xtol": xtol, "maxiter": maxfun, "disp": disp}

    res = _minimize_scalar_bounded(func, (x1, x2), args, **options)
    if full_output:
        return res["x"], res["fun"], res["status"], res["nfev"]
    else:
        return res["x"]


def _minimize_scalar_bounded(
    func, bounds, args=(), xtol=1e-5, maxiter=500, disp=0, **unknown_options
):
    maxfun = maxiter
    # Test bounds are of correct form
    if len(bounds) != 2:
        raise ValueError("bounds must have two elements.")
    x1, x2 = bounds

    if not (is_array_scalar(x1) and is_array_scalar(x2)):
        raise ValueError("Optimisation bounds must be scalars" " or array scalars.")
    if x1 > x2:
        raise ValueError("The lower bound exceeds the upper bound.")

    flag = 0
    header = " Func-count     x          f(x)          Procedure"
    step = "       initial"

    sqrt_eps = sqrt(2.2e-16)
    golden_mean = 0.5 * (3.0 - sqrt(5.0))
    a, b = x1, x2
    fulc = a + golden_mean * (b - a)
    nfc, xf = fulc, fulc
    rat = e = 0.0
    x = xf
    fx = func(x, *args)
    num = 1
    fmin_data = (1, xf, fx)

    ffulc = fnfc = fx
    xm = 0.5 * (a + b)
    tol1 = sqrt_eps * numpy.abs(xf) + xtol / 3.0
    tol2 = 2.0 * tol1

    if disp > 2:
        print(" ")
        print(header)
        print("%5.0f   %12.6g %12.6g %s" % (fmin_data + (step,)))

    while numpy.abs(xf - xm) > (tol2 - 0.5 * (b - a)):
        golden = 1
        # Check for parabolic fit
        if numpy.abs(e) > tol1:
            golden = 0
            r = (xf - nfc) * (fx - ffulc)
            q = (xf - fulc) * (fx - fnfc)
            p = (xf - fulc) * q - (xf - nfc) * r
            q = 2.0 * (q - r)
            if q > 0.0:
                p = -p
            q = numpy.abs(q)
            r = e
            e = rat

            # Check for acceptability of parabola
            if (
                (numpy.abs(p) < numpy.abs(0.5 * q * r))
                and (p > q * (a - xf))
                and (p < q * (b - xf))
            ):
                rat = (p + 0.0) / q
                x = xf + rat
                step = "       parabolic"

                if ((x - a) < tol2) or ((b - x) < tol2):
                    si = numpy.sign(xm - xf) + ((xm - xf) == 0)
                    rat = tol1 * si
            else:  # do a golden section step
                golden = 1

        if golden:  # Do a golden-section step
            if xf >= xm:
                e = a - xf
            else:
                e = b - xf
            rat = golden_mean * e
            step = "       golden"

        si = numpy.sign(rat) + (rat == 0)
        x = xf + si * numpy.max([numpy.abs(rat), tol1])
        fu = func(x, *args)
        num += 1
        fmin_data = (num, x, fu)
        if disp > 2:
            print("%5.0f   %12.6g %12.6g %s" % (fmin_data + (step,)))

        if fu <= fx:
            if x >= xf:
                a = xf
            else:
                b = xf
            fulc, ffulc = nfc, fnfc
            nfc, fnfc = xf, fx
            xf, fx = x, fu
        else:
            if x < xf:
                a = x
            else:
                b = x
            if (fu <= fnfc) or (nfc == xf):
                fulc, ffulc = nfc, fnfc
                nfc, fnfc = x, fu
            elif (fu <= ffulc) or (fulc == xf) or (fulc == nfc):
                fulc, ffulc = x, fu

        xm = 0.5 * (a + b)
        tol1 = sqrt_eps * numpy.abs(xf) + xtol / 3.0
        tol2 = 2.0 * tol1

        if num >= maxfun:
            flag = 1
            break

    fval = fx
    if disp > 0:
        _endprint(x, flag, fval, maxfun, xtol, disp)

    result = Result(
        fun=fval,
        status=flag,
        success=(flag == 0),
        message={
            0: "Solution found.",
            1: "Maximum number of function calls " "reached.",
        }.get(flag, ""),
        x=xf,
        nfev=num,
    )

    return result


def is_array_scalar(x):
    """Test whether `x` is either a scalar or an array scalar."""
    return len(atleast_1d(x) == 1)


_epsilon = sqrt(numpy.finfo(float).eps)


def _endprint(x, flag, fval, maxfun, xtol, disp):
    if flag == 0:
        if disp > 1:
            print(
                "\nOptimization terminated successfully;\n"
                "The returned value satisfies the termination criteria\n"
                "(using xtol = ",
                xtol,
                ")",
            )
    if flag == 1:
        if disp:
            print(
                "\nMaximum number of function evaluations exceeded --- "
                "increase maxfun argument.\n"
            )
    return
