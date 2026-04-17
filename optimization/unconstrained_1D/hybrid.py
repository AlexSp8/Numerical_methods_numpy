
from typing import Callable, Tuple, List
import numpy as np
import scipy

def brent(f: Callable[[float], float], a0: float, b0: float, mode: str = 'min',
    output: bool = False, eps: float = 1e-8, k_max: int = 1000) -> Tuple[float, float]:
    """Returns the optimum of a function f(x) in an interval [a, b]
    using Brent's hybrid method"""

    # if mode == 'min':
    #     g = f
    # else:
    #     g = lambda x: -f(x)
    # x = scipy.optimize.brent(g, brack=(a0,b0), tol=eps) # unrestricted
    # x = scipy.optimize.fminbound(g, x1=a0, x2=b0, xtol=eps)
    # r = scipy.optimize.minimize_scalar(g, bounds=(a0,b0), method='bounded', options={'xatol': eps})
    # x = r.x
    # return float(x), float(f(x))

    a, b = min(a0,b0), max(a0,b0)
    phi = (3.0 - (5.0)**0.5) / 2.0
    if mode == 'min':
        s = +1.0
    else:
        s = -1.0

    # x: best point, w: second best, v: third best
    x = w = v = a + phi*(b - a)
    fx = fw = fv = f(x)*s

    d = e = 0.0  # d: current step, e: step before last

    for k in range(1, k_max + 1):

        mid = (a + b)/2.0
        tol = eps*abs(x) + eps/10.0  # Tolerance relative to x

        err = abs(x-mid)
        if abs(x - mid) <= (2.0*tol - (b - a)/2.0):
            if output:
                print(f'k = {k}. Rel. Error: {f"{err:.4e}"}')
            break

        is_parabolic = False
        if abs(e) > tol:

            r = (x - w)*(fx - fv)
            q = (x - v)*(fx - fw)
            p = (x - v)*q - (x - w)*r #nom
            q = 2.0*(q - r) #denom

            if q > 0.0:
                p = -p
            q = abs(q)

            if abs(p) < abs(0.5*q*e) and p > q*(a - x) and p < q*(b - x):
                is_parabolic = True
                e = d
                d = p/q

        if not is_parabolic:
            if x < mid:
                e = b - x
            else:
                e = a - x
            d = phi*e

        if abs(d) < tol:
            d = float(np.sign(d)*tol)

        u = x + d
        fu = f(u)*s
        # # If u is too close to a or b
        # if (u - a) < 2*tol or (b - u) < 2*tol:
        #     if x < mid:
        #         d = tol
        #     else:
        #         d = -tol
        #     u = x + d
        #     fu = f(u)*s

        u_is_best = fu <= fx
        if u_is_best:

            if u >= x:
                a = x
            else:
                b = x

            v, fv = w, fw
            w, fw = x, fx
            x, fx = u, fu

        else:

            if u < x:
                a = u
            else:
                b = u

            u_over_w = (fu <= fw) or (w == x)
            u_over_v = (fu <= fv) or (v == x) or (v == w)
            if u_over_w:
                v, fv = w, fw
                w, fw = u, fu
            elif u_over_v:
                v, fv = u, fu

    return x, fx*s

def multi_hybrid(f: Callable[[float], float], a: float, b: float,
    n: int, method: str, mode: str = 'min', output: bool = False) -> List[Tuple[float, float]]:
    """Returns a list of extremes of a function f(x) in an interval [a, b]
    using a hybrid method across n intervals"""

    if n <= 0:
        raise ValueError('n <= 0 in multi-bracketing!')

    extremes = []
    dx = (b-a)/n
    for i in range(n):

        a_int = a + i*dx
        b_int = a_int + dx

        if method == 'brent':
            x, fx = brent(f, a_int, b_int, mode, output)
        else:
            raise ValueError('Invalid bracketing method!')

        extremes.append((x, fx))

    return extremes
