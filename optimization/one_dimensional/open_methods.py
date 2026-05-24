
from typing import Callable
import scipy

def parabolic_interpolation(f: Callable[[float], float],
    a: float, b: float, tol: float = 1e-8, k_max: int = 1000) -> float:
    """Returns the optimum of a function in an interval
    using the parabolic interpolation method.

    Args:
        f: the function to optimize
        a: the lower limit of the interval
        b: the upper limit of the interval
        tol (optional): threshold for convergence
        k_max (optional): maximum iterations
    Returns:
        The point of the function optimum."""

    x0, x1, x2 = a, (a+b)/2.0, b

    f0, f1, f2 = f(x0), f(x1), f(x2)

    for k in range(1, k_max+1):

        nom = ((x1 - x0)**2)*(f1 - f2) - ((x1 - x2)**2)*(f1 - f0)
        denom = 2.0*((x1 - x0)*(f1 - f2) - (x1 - x2)*(f1 - f0))
        if abs(denom) < 1e-12:
            x_opt = x1 + tol
        else:
            x_opt = x1 - (nom/denom)

        f_opt = f(x_opt)

        if x_opt > x1:
            if f_opt < f1:
                # x_opt is better, remove x0
                x0, f0, x1, f1 = x1, f1, x_opt, f_opt
            else:
                x2, f2 = x_opt, f_opt
        else:
            if f_opt < f1:
                # x_opt is better, remove x2
                x2, f2, x1, f1 = x1, f1, x_opt, f_opt
            else:
                x0, f0 = x_opt, f_opt

        x_int = x2 - x0
        err = abs( x_int/(x_opt+1e-12) )
        if err < tol:
            print(f'k = {k}. Rel. Error: {f"{err:.4e}"}')
            break

    return x_opt

def secant(f: Callable[[float], float], x0: float,
    df: Callable[[Callable[[float], float], float, float, float], float],
    tol: float = 1e-8, k_max: int = 1000, r: float = 1.0) -> float:
    """Returns the optimum of a function around an initial guess
    using the Secant method.

    Args:
        f: the function to optimize
        x0: initial guess
        df: the function to evaluate the derivative of f
        tol (optional): threshold for convergence
        k_max (optional): maximum iterations
        r (optional): relaxation factor
    Returns:
        The point of the function optimum."""

    # df_x = lambda x: df(f, x)
    # x = scipy.optimize.newton(func=df_x, x0=x0, tol=tol)
    # return float(x)

    xo = x0 + 0.01
    xc = x0
    dfo_x = df(f, xo)
    for k in range(1, k_max+1):

        df_x = df(f, xc)

        x = xc - r*df_x*(xc-xo)/(df_x-dfo_x+1e-12)

        # print(f"k: {k}, x: {x}, f'(x): {df_x}")

        err = [abs(df_x), abs( (x-xc)/(x+1e-12) )]

        if any(e < tol for e in err):
            print(f'k = {k}. Errors (df, rel): {[f"{e:.4e}" for e in err]}')
            break

        xo, dfo_x = xc, df_x
        xc = x

    return x

def newton(f: Callable[[float], float], x0: float,
    df: Callable[[Callable[[float], float], float, float, float], float],
    d2f: Callable[[Callable[[float], float], float, float, float], float],
    tol: float = 1e-8, k_max: int = 1000, r: float = 1.0) -> float:
    """Returns the optimum of a function around an initial guess
    using the Newton-Raphson method.

    Args:
        f: the function to be optimized
        x0: the initial guess
        df (optional): the function to evaluate the derivative of f
        d2f (optional): the function to evaluate the 2nd derivative of f
        tol (optional): threshold for convergence
        k_max (optional): maximum iterations
        r (optional): relaxation factor
    Returns:
        The point of the function optimum."""

    # df_x = lambda x: df(f, x)
    # d2f_x = lambda x: d2f(f, x)

    # x = scipy.optimize.newton(func=df_x, x0=x0, fprime=d2f_x, tol=tol)
    # return float(x)

    x = x0
    for k in range(1, k_max+1):

        x_old = x

        df_x = df(f, x)

        d2f_x = d2f(f, x)

        x = x_old - r*df_x/(d2f_x+1e-12)

        # print(f"k: {k}, x: {x}, f'(x): {df_x}")

        err = [abs(df_x), abs( (x-x_old)/(x+1e-12) )]

        if any(e < tol for e in err):
            print(f'k = {k}. Errors (df, rel): {[f"{e:.4e}" for e in err]}')
            break

    return x
