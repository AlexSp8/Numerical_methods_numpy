
from typing import Callable, Tuple, List
import scipy

from differentiation import forward_fd as ffd

def parabolic_interpolation(f: Callable[[float], float],
    a: float, b: float, mode: str = 'min',
    eps: float = 1e-8, k_max: int = 1000) -> Tuple[float]:
    """Returns the extreme of a function f(x) in an interval [a, b]
    using the parabolic interpolation method"""

    if mode == 'min':
        s = +1.0
    else:
        s = -1.0

    x0, x1, x2 = a, (a+b)/2.0, b

    f0, f1, f2 = f(x0), f(x1), f(x2)

    for k in range(1, k_max+1):

        nom = ((x1 - x0)**2)*(f1 - f2) - ((x1 - x2)**2)*(f1 - f0)
        denom = 2.0*((x1 - x0)*(f1 - f2) - (x1 - x2)*(f1 - f0))
        if abs(denom) < 1e-12:
            x_opt = x1 + eps
        else:
            x_opt = x1 - (nom/denom)

        f_opt = f(x_opt)

        if x_opt > x1:
            if s*f_opt < s*f1:
                x0, f0, x1, f1 = x1, f1, x_opt, f_opt
            else:
                x2, f2 = x_opt, f_opt
        else:
            if s*f_opt < s*f1:
                x2, f2, x1, f1 = x1, f1, x_opt, f_opt
            else:
                x0, f0 = x_opt, f_opt

        x_int = x2 - x0
        err = abs( x_int/(x_opt+1e-12) )
        if err < eps:
            print(f'k = {k}. Rel. Error: {f"{err:.4e}"}')
            break

    return x_opt, f_opt

def secant(f: Callable[[float], float], x0: float,
    eps: float = 1e-8, k_max: int = 1000, r: float = 1.0) -> Tuple[float]:
    """Returns the extreme of a function f(x)
    around an initial guess, x0, using the Secant method"""

    # df = lambda x: ffd.df_h(f, x)
    # x = scipy.optimize.newton(func=df, x0=x0, tol=eps)
    # return float(x), float(f(x))

    xo = x0 + 0.01
    xc = x0
    dfo = ffd.df_h(f, xo)
    for k in range(1, k_max+1):

        df = ffd.df_h(f, xc)

        x = xc - r*df*(xc-xo)/(df-dfo+1e-12)

        # print(f"k: {k}, x: {x}, f'(x): {df_x}")

        err = [abs(df), abs( (x-xc)/(x+1e-12) )]

        if any(e < eps for e in err):
            print(f'k = {k}. Errors (df, rel): {[f"{e:.4e}" for e in err]}')
            break

        xo, dfo = xc, df
        xc = x

    return x, f(x)

def newton(f: Callable[[float], float], x0: float,
    eps: float = 1e-8, k_max: int = 1000, r: float = 1.0) -> Tuple[float]:
    """Returns the extreme of a function f(x)
    around an initial guess, x0, using the Newton-Raphson method"""

    # df = lambda x: ffd.df_h(f, x)
    # d2f = lambda x: ffd.d2f_h(f, x)

    # x = scipy.optimize.newton(func=df, x0=x0, fprime=d2f, tol=eps)
    # return float(x), float(f(x))

    x = x0
    for k in range(1, k_max+1):

        x_old = x

        df = ffd.df_h(f, x)

        d2f = ffd.d2f_h(f, x)

        x = x_old - r*df/(d2f+1e-12)

        # print(f"k: {k}, x: {x}, f'(x): {df}")

        err = [abs(df), abs( (x-x_old)/(x+1e-12) )]

        if any(e < eps for e in err):
            print(f'k = {k}. Errors (df, rel): {[f"{e:.4e}" for e in err]}')
            break

    return x, f(x)

def multi_open(f: Callable[[float], float],
    a: float, b: float, n: int, method: str,
    mode: str = 'min') -> List[Tuple[float, float]]:
    """Returns a list of extremes of a function f(x) in an interval [a, b]
    using an open method across n intervals"""

    if n <= 0:
        raise ValueError('n <= 0 in multi-bracketing!')

    extremes = []
    dx = (b-a)/n
    for i in range(n):

        a_int = a + i*dx
        b_int = a_int + dx

        if method == 'parabolic_interpolation':
            x, fx = parabolic_interpolation(f, a_int, b_int, mode)
        elif method == 'secant':
            x, fx = secant(f, a_int)
        elif method == 'newton':
            x, fx = newton(f, a_int)
        else:
            raise ValueError('Invalid bracketing method!')

        extremes.append((x, fx))

    return extremes
