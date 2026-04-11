
from typing import Callable
import math
import numpy as np

def bisection(f: Callable[[float], float], a: float, b: float,
    eps: float = 1e-8, k_max: int = 1000) -> float|None:
    """Returns the root, f(xr) = 0, of a function f(x)
    in an interval [a, b] using the bisection method"""

    if f(a)*f(b) >= 0:
        return None

    n = math.log2((b-a)/eps)
    print(f"Expected iterations for absolute error {eps}: {math.ceil(n)}")

    xl, xu = a, b
    xr, f_xl = a, f(xl)

    for i in range(1, k_max+1):

        xr_old = xr
        xr = (xl+xu)/2

        f_xr = f(xr)

        # print(f"i: {i}, xr: {xr}, f(xr): {f_xr}")

        err = [abs(f_xr), abs(xr - xr_old), abs( (xr-xr_old)/(xr+eps) )]

        if all(e < eps for e in err):
            print(f'i = {i}. Errors (f, abs, rel): {[f"{e:.4e}" for e in err]}')
            break

        if f_xl*f_xr < 0:
            xu = xr
        else:
            xl = xr
            f_xl = f_xr

    return xr

def false_position(f: Callable[[float], float], a: float, b: float,
    eps: float = 1e-8, k_max: int = 1000) -> float|None:
    """Returns the root, f(xr) = 0, of a function f(x)
    in an interval [a, b] using the false position method"""

    if f(a)*f(b) >= 0:
        return None

    xl, xu = a, b
    xr, f_xl, f_xu = a, f(xl), f(xu)
    il, iu = 0, 0

    for i in range(1, k_max+1):

        xr_old = xr
        xr = xu - ( f_xu*(xl-xu) )/( f_xl-f_xu+eps )

        f_xr = f(xr)

        # print(f"i: {i}, xr: {xr}, f(xr): {f_xr}")

        err = [abs(f_xr), abs(xr - xr_old), abs( (xr-xr_old)/(xr+eps) )]

        if all(e < eps for e in err):
            print(f'i = {i}. Errors (f, abs, rel): {[f"{e:.4e}" for e in err]}')
            break

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

    return xr

def multi_bracketing(f: Callable[[float], float],
    a: float, b: float, n: int, method: str) -> np.ndarray:
    """Returns a list of roots of a function f(x) in an interval [a, b]
    using the bisection or regular position method across n intervals"""

    if n <= 0:
        raise ValueError('n <= 0 in multi-bracketing!')

    roots = []
    dx = (b-a)/n
    for i in range(n):
        a_int = a + i*dx
        b_int = a_int + dx
        if f(a_int)*f(b_int) > 0:
            continue

        if method == 'bisection':
            xr = bisection(f, a_int, b_int)
        elif method == 'false position':
            xr = false_position(f, a_int, b_int)
        else:
            raise ValueError('Invalid bracketing method!')
        roots.append(xr)
    return np.array(roots)
