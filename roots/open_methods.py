
from typing import Callable, List
import scipy

from differentiation import forward_fd as ffd

def fixed_point(f: Callable[[float], float], x0: float,
    eps: float = 1e-8, k_max: int = 1000) -> float|None:
    """Returns the root, f(x) = 0, of a function f(x)
    around an initial guess, x0, using the fixed-point iteration method"""

    x = x0
    for k in range(1, k_max+1):

        x_old = x

        f_x = f(x)

        x = f_x + x

        # print(f"k: {k}, x: {x}, f(x): {f_x}")

        err = [abs(f_x), abs( (x-x_old)/(x+1e-12) )]

        if all(e < eps for e in err):
            print(f'k = {k}. Errors (f, rel): {[f"{e:.4e}" for e in err]}')
            return x

    return None

def secant(f: Callable[[float], float], x0: float, x1: float,
    eps: float = 1e-8, k_max: int = 1000) -> float|None:
    """Returns the root, f(x) = 0, of a function f(x)
    around an initial guess, x0, using the Secant method"""

    # return scipy.optimize.newton(f, x0=x0, x1=x1)

    f0, f1 = f(x0), f(x1)

    for k in range(1, k_max+1):

        try:
            x2 = x1 - f1*(x1-x0)/(f1-f0)
        except ZeroDivisionError:
            return None

        f2 = f(x2)
        err = [abs(f2), abs( (x2-x1)/(x2+1e-12) )]
        if all(e < eps for e in err):
            print(f'k = {k}. Errors (f, rel): {[f"{e:.4e}" for e in err]}')
            return x2

        x0, f0, x1, f1 = x1, f1, x2, f2

    return None

def iqi(f: Callable[[float], float], x0: float, x1: float, x2: float,
    eps: float = 1e-8, k_max: int = 1000) -> float|None:
    """Returns the root, f(x) = 0, of a function f(x)
    around an initial guess, x0, using the IQI method"""

    f0, f1, f2 = f(x0), f(x1), f(x2)

    for k in range(1, k_max+1):

        try:
            term0 = (x0*f1*f2)/((f0 - f1)*(f0 - f2))
            term1 = (x1*f0*f2)/((f1 - f0)*(f1 - f2))
            term2 = (x2*f0*f1)/((f2 - f0)*(f2 - f1))

            x3 = term0 + term1 + term2
        except ZeroDivisionError:
            return None

        f3 = f(x3)

        err = [abs(f3), abs( (x3-x2)/(x3+1e-12) )]
        if all(e < eps for e in err):
            print(f'k = {k}. Errors (f, rel): {[f"{e:.4e}" for e in err]}')
            return x3

        x0, f0 = x1, f1
        x1, f1 = x2, f2
        x2, f2 = x3, f3

    return None

def newton_raphson(f: Callable[[float], float], x0: float,
    eps: float = 1e-8, k_max: int = 1000, r: float = 1.0) -> float|None:
    """Returns the root, f(x) = 0, of a function f(x)
    around an initial guess, x0, using the Newton-Raphson method"""

    # return scipy.optimize.newton(f, x0=x0,
    #     fprime=lambda x: ffd.df_h(f,x))
    # return scipy.optimize.newton(f, x0=x0, fprime=lambda x: ffd.df_h(f,x),
    #     fprime2=lambda x: ffd.d2f_h(f,x)) #Halley

    x = x0
    for k in range(1, k_max+1):

        x_old = x

        f_x = f(x)
        df = ffd.df_h(f, x)

        x = x_old - r*f_x/(df+1e-12)

        # print(f"k: {k}, x: {x}, f(x): {f_x}")

        err = [abs(f_x), abs( (x-x_old)/(x+1e-12) )]

        if all(e < eps for e in err):
            print(f'k = {k}. Errors (f, rel): {[f"{e:.4e}" for e in err]}')
            return x

    return None

def ralston_rabinowitz(f: Callable[[float], float], x0: float,
    eps: float = 1e-8, k_max: int = 1000, r: float = 1.0) -> float|None:
    """Returns the root, f(x) = 0, of a function f(x)
    around an initial guess, x0, using the Ralston-Rabinowitz method"""

    x = x0
    for k in range(1, k_max+1):

        x_old = x

        f_x = f(x)

        df_x = ffd.df_h(f, x)

        d2f_x = ffd.d2f_h(f, x)

        x = x_old - r*f_x*df_x/( (df_x**2)-f_x*d2f_x+eps )

        # print(f"k: {k}, x: {x}, f(x): {f_x}")

        err = [abs(f_x), abs( (x-x_old)/(x+1e-12) )]

        if all(e < eps for e in err):
            print(f'k = {k}. Errors (f, rel): {[f"{e:.4e}" for e in err]}')
            return x

    return None

def multi_open(f: Callable[[float], float],
    a: float, b: float, n: int, method: str) -> List[float|None]:
    """Returns a list of roots of a function f(x) in an interval [a, b]
    using an open method across n intervals"""

    if n <= 0:
        raise ValueError('n <= 0 in multi-bracketing!')

    roots = []
    dx = (b-a)/n
    for i in range(n):

        x0 = a + i*dx

        if method == 'fixed_point':
            xr = fixed_point(f, x0)
        elif method == 'secant':
            x1 = x0+dx
            xr = secant(f, x0, x1)
        elif method == 'iqi':
            x1, x2 = x0+dx/2, x0+dx
            xr = iqi(f, x0, x1, x2)
        elif method == 'newton_raphson':
            xr = newton_raphson(f, x0)
        elif method == 'ralston_rabinowitz':
            xr = ralston_rabinowitz(f, x0)
        else:
            raise ValueError('Invalid open method!')
        if xr is not None:
            roots.append(float(xr))
    return roots
