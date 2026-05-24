
from typing import Callable
import math
import scipy

from roots import hybrid, open_methods

def bisection(f: Callable[[float], float], a: float, b: float,
    tol: float = 1e-8, k_max: int = 1000) -> float | None:
    """Returns the root of a function in an interval using the
    bisection method.

    Args:
        f: an algebraic function
        a: the lower limit of the interval
        b: the upper limit of the interval
        tol: the numerical threshold for convergence
        k_max: the maximum number of iterations

    Returns:
        A root (float) if it exists, otherwise None.
    """

    # return scipy.optimize.bisect(f, a, b)

    if f(a)*f(b) >= 0:
        return None

    # n = math.log2((b-a)/tol)
    # print(f"Expected iterations for absolute error {tol}: {math.ceil(n)}")

    xl, xu = a, b
    xr, f_xl = a, f(xl)

    for k in range(1, k_max+1):

        xr_old = xr
        xr = (xl+xu)/2

        f_xr = f(xr)

        # print(f"i: {i}, xr: {xr}, f(xr): {f_xr}")

        err = [abs(f_xr), abs( (xr-xr_old)/(xr+1e-12) )]

        if all(e < tol for e in err):
            print(f'k = {k}. Errors (f, rel): {[f"{e:.4e}" for e in err]}')
            return xr

        if f_xl*f_xr < 0:
            xu = xr
        else:
            xl = xr
            f_xl = f_xr

    return None

def false_position(f: Callable[[float], float], a: float, b: float,
    tol: float = 1e-8, k_max: int = 1000) -> float | None:
    """Returns the root of a function in an interval using the
    false position method.

    Args:
        f: an algebraic function
        a: the lower limit of the interval
        b: the upper limit of the interval
        tol: the numerical threshold for convergence
        k_max: the maximum number of iterations

    Returns:
        A root (float) if it exists, otherwise None.
    """

    if f(a)*f(b) >= 0:
        return None

    xl, xu = a, b
    xr, f_xl, f_xu = a, f(xl), f(xu)
    il, iu = 0, 0

    for k in range(1, k_max+1):

        xr_old = xr
        try:
            xr = xu - ( f_xu*(xl-xu) )/(f_xl-f_xu)
        except ZeroDivisionError:
            return None

        f_xr = f(xr)

        # print(f"i: {i}, xr: {xr}, f(xr): {f_xr}")

        err = [abs(f_xr), abs( (xr-xr_old)/(xr+1e-12) )]

        if all(e < tol for e in err):
            print(f'k = {k}. Errors (f, rel): {[f"{e:.4e}" for e in err]}')
            return xr

        if f_xl*f_xr < 0:
            xu = xr
            f_xu = f_xr
            il += 1
            if il >= 2:
                f_xl /= 2
                il = 0
        else:
            xl = xr
            f_xl = f_xr
            iu += 1
            if iu >= 2:
                f_xu /= 2
                iu = 0

    return None

def multi_bracketing(f: Callable[[float], float],
    a: float, b: float, n: int, method: str,
    df: Callable[[Callable[[float], float], float, float, float], float] = None,
    d2f: Callable[[Callable[[float], float], float, float, float], float] = None
    ) -> list[float|None]:
    """Returns the list of roots of a function in an interval,
    divided into smaller intervals, using a root finding method.

    Args:
        f: an algebraic function
        a: the lower limit of the total interval
        b: the upper limit of the total interval
        n: the number of intervals
        method: the root finding method
        df (optional): the derivative of the function calculated with finite differences
        d2f (optional): the second derivative of the function calculated with finite differences

    Returns:
        A list of roots (float) if they exist, otherwise None.
    """

    if n <= 0:
        raise ValueError('n <= 0 in multi-bracketing!')

    roots = []
    dx = (b-a)/n
    for i in range(n):

        a_int = a + i*dx
        b_int = a_int + dx

        x0 = a_int
        if method == 'fixed_point':
            xr = open_methods.fixed_point(f, x0)
        elif method == 'secant':
            x1 = x0+dx
            xr = open_methods.secant(f, x0, x1)
        elif method == 'iqi':
            x1, x2 = x0+dx/2, x0+dx
            xr = open_methods.iqi(f, x0, x1, x2)
        elif method == 'newton_raphson':
            xr = open_methods.newton_raphson(f, x0, df)
        elif method == 'ralston_rabinowitz':
            xr = open_methods.ralston_rabinowitz(f, x0, df, d2f)
        else:
            pass

        if f(a_int)*f(b_int) > 0:
            continue

        if method == 'bisection':
            xr = bisection(f, a_int, b_int)
        elif method == 'false_position':
            xr = false_position(f, a_int, b_int)
        elif method == 'brent':
            xr = hybrid.brent(f, a_int, b_int)
        elif method == 'ridders':
            xr = hybrid.ridders(f, a_int, b_int)
        elif method == 'chandrupatla':
            xr = hybrid.chandrupatla(f, a_int, b_int)
        else:
            pass

        if xr is not None:
            roots.append(float(xr))

    return roots
