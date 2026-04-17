
from typing import Callable, List
import scipy

def brent(f: Callable[[float], float], a0: float, b0: float,
    eps: float = 1e-8, k_max: int = 1000) -> float|None:
    """Returns the root, f(x) = 0, of a function f(x)
    in an interval [a,b] using Brent's method"""

    # return scipy.optimize.brentq(f, a0, b0)

    a, b = a0, b0

    fa, fb = f(a), f(b)
    if fa*fb >= 0:
        return None

    if abs(fa) < abs(fb):
        a, b = b, a
        fa, fb = fb, fa

    c, fc = a, fa
    bflag = True
    d = 0.0

    for k in range(1, k_max+1):

        err = [abs(fb), abs( (b-a)/(b+1e-12) )]
        if all(e < eps for e in err):
            print(f'k = {k}. Errors (f, rel): {[f"{e:.4e}" for e in err]}')
            return b

        if fa != fc and fb != fc:
            # IQI
            s = (a*fb*fc)/((fa-fb)*(fa-fc)) + \
                (b*fa*fc)/((fb-fa)*(fb-fc)) + \
                (c*fa*fb)/((fc-fa)*(fc-fb))
        else:
            # Secant
            s = b - fb*(b-a)/(fb-fa)

        condition1 = not ((3*a+b)/4 <= s <= b)
        condition2 = bflag and (abs(s-b) >= abs(b-c)/2)
        condition3 = not bflag and (abs(s-b) >= abs(c-d)/2)
        condition4 = bflag and (abs(b-c) < eps)
        condition5 = not bflag and (abs(c-d) < eps)
        if any([condition1, condition2, condition3, condition4, condition5]):
            # Candidate, s, rejected
            s = (a+b)/2
            bflag = True
        else:
            bflag = False

        fs = f(s)

        d = c
        # Update previous best guess, c
        c, fc = b, fb

        if fa*fs < 0:
            # New best, b, is the new s
            b, fb = s, fs
        else:
            a, fa = s, fs

        if abs(fa) < abs(fb):
            a, b = b, a
            fa, fb = fb, fa

    return None

def ridders(f: Callable[[float], float], a: float, b: float,
    eps: float = 1e-8, k_max: int = 1000) -> float|None:
    """Returns the root, f(x) = 0, of a function f(x)
    in an interval [a,b] using Ridders' method"""

    # return scipy.optimize.ridder(f, a, b)

    x1, x2 = a, b
    f1, f2 = f(x1), f(x2)
    if f1*f2 >= 0:
        return None

    for k in range(1, k_max+1):

        m = (x1 + x2)/2
        fm = f(m)

        s = (fm**2 - f1*f2)**0.5

        dx = (m - x1)*fm/(s+1e-12)
        if f1 < f2:
            xr = m - dx
        else:
            xr = m + dx

        fr = f(xr)
        err = [abs(fr), abs( (x2-x1)/(x2+1e-12) )]
        if all(e < eps for e in err):
            print(f'k = {k}. Errors (f, rel): {[f"{e:.4e}" for e in err]}')
            return xr

        if fm*fr < 0:
            x1, x2, f1, f2 = m, xr, fm, fr
        elif f1*fr < 0:
            x2, f2 = xr, fr
        else:
            x1, f1 = xr, fr

    return None

def chandrupatla(f: Callable[[float], float], a0: float, b0: float,
    eps: float = 1e-8, k_max: int = 1000) -> float|None:
    """Returns the root, f(x) = 0, of a function f(x)
    in an interval [a,b] using Chandrupatla's method"""

    # from scipy.optimize import elementwise
    # r = elementwise.find_root(f, (a0, b0))
    # return r.x

    a, b = a0, b0
    fa, fb = f(a), f(b)
    if fa*fb >= 0:
        return None

    c, fc = a, fa
    for k in range(k_max):

        err = [abs(fb), abs( (b-a)/(b+1e-12) )]
        if all(e < eps for e in err):
            print(f'k = {k}. Errors (f, rel): {[f"{e:.4e}" for e in err]}')
            return b

        xi = (a - b)/(c - b + 1e-12)
        phi = (fa - fb)/(fc - fb + 1e-12)

        if 0.0 < phi < 1.0 and xi < 1.0 - 1.0/(1.0 - phi):
            # IQI step
            term0 = (a*fb*fc)/((fa - fb)*(fa - fc))
            term1 = (b*fa*fc)/((fb - fa)*(fb - fc))
            term2 = (c*fa*fb)/((fc - fa)*(fc - fb))
            s = term0 + term1 + term2
        else:
            s = (a + b)/2

        fs = f(s)

        c, fc = a, fa
        if fa*fs < 0:
            b, fb = s, fs
        else:
            a, fa = s, fs

        if abs(fa) < abs(fb):
            a, b = b, a
            fa, fb = fb, fa

    return None

def multi_hybrid(f: Callable[[float], float],
    a: float, b: float, n: int, method: str) -> List[float|None]:
    """Returns a list of roots of a function f(x) in an interval [a, b]
    using a hybrid method across n intervals"""

    if n <= 0:
        raise ValueError('n <= 0 in multi-bracketing!')

    roots = []
    dx = (b-a)/n
    for i in range(n):

        a_int = a + i*dx
        b_int = a_int + dx
        if f(a_int)*f(b_int) > 0:
            continue

        if method == 'brent':
            xr = brent(f, a_int, b_int)
        elif method == 'ridders':
            xr = ridders(f, a_int, b_int)
        elif method == 'chandrupatla':
            xr = chandrupatla(f, a_int, b_int)
        else:
            raise ValueError('Invalid bracketing method!')
        if xr is not None:
            roots.append(float(xr))
    return roots
