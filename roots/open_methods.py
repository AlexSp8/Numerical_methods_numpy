
from typing import Callable
import scipy

def fixed_point(f: Callable[[float], float], x0: float,
    tol: float = 1e-8, k_max: int = 1000) -> float | None:
    """Returns the root of a function around an initial guess
    using the fixed-point iteration method.

    Args:
        f: an algebraic function
        x0: the initial guess
        tol: the numerical threshold for convergence
        k_max: the maximum number of iterations

    Returns:
        A root (float) if it exists, otherwise None
    """

    x = x0
    for k in range(1, k_max+1):

        x_old = x

        f_x = f(x)

        x = f_x + x

        # print(f"k: {k}, x: {x}, f(x): {f_x}")

        err = [abs(f_x), abs( (x-x_old)/(x+1e-12) )]

        if all(e < tol for e in err):
            print(f'k = {k}. Errors (f, rel): {[f"{e:.4e}" for e in err]}')
            return x

    return None

def secant(f: Callable[[float], float], x0: float, x1: float,
    tol: float = 1e-8, k_max: int = 1000) -> float | None:
    """Returns the root of a function around two initial guesses
    using the Secant method.

    Args:
        f: an algebraic function
        x0: the first initial guess
        x1: the second initial guess
        tol: the numerical threshold for convergence
        k_max: the maximum number of iterations

    Returns:
        A root (float) if it exists, otherwise None
    """

    # return scipy.optimize.newton(f, x0=x0, x1=x1)

    f0, f1 = f(x0), f(x1)

    for k in range(1, k_max+1):

        try:
            x2 = x1 - f1*(x1-x0)/(f1-f0)
        except ZeroDivisionError:
            return None

        f2 = f(x2)
        err = [abs(f2), abs( (x2-x1)/(x2+1e-12) )]
        if all(e < tol for e in err):
            print(f'k = {k}. Errors (f, rel): {[f"{e:.4e}" for e in err]}')
            return x2

        x0, f0, x1, f1 = x1, f1, x2, f2

    return None

def iqi(f: Callable[[float], float], x0: float, x1: float, x2: float,
    tol: float = 1e-8, k_max: int = 1000) -> float | None:
    """Returns the root of a function around three initial points
    using the Inverse Quadratic Interpolation (IQI) method.

    Args:
        f: an algebraic function
        x0: the first initial point
        x1: the second initial point
        x2: the third initial point
        tol: the numerical threshold for convergence
        k_max: the maximum number of iterations

    Returns:
        A root (float) if it exists, otherwise None
    """

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
        if all(e < tol for e in err):
            print(f'k = {k}. Errors (f, rel): {[f"{e:.4e}" for e in err]}')
            return x3

        x0, f0 = x1, f1
        x1, f1 = x2, f2
        x2, f2 = x3, f3

    return None

def newton_raphson(f: Callable[[float], float], x0: float,
    df: Callable[[Callable[[float], float], float, float, float], float],
    tol: float = 1e-8, k_max: int = 1000, r: float = 1.0) -> float | None:
    """Returns the root of a function around an initial guess
    using the Newton Raphson method.

    Args:
        f: an algebraic function
        df: the derivative of the function calculated with finite differences
        x0: the initial guess
        tol: the numerical threshold for convergence
        k_max: the maximum number of iterations

    Returns:
        A root (float) if it exists, otherwise None
    """

    # return scipy.optimize.newton(f, x0=x0,
    #     fprime=lambda x: df(f,x))

    x = x0
    for k in range(1, k_max+1):

        x_old = x

        f_x = f(x)
        df_x = df(f, x)

        try:
            x = x_old - r*f_x/df_x
        except ZeroDivisionError:
            return None

        # print(f"k: {k}, x: {x}, f(x): {f_x}")

        err = [abs(f_x), abs( (x-x_old)/(x+1e-12) )]

        if all(e < tol for e in err):
            print(f'k = {k}. Errors (f, rel): {[f"{e:.4e}" for e in err]}')
            return x

    return None

def ralston_rabinowitz(f: Callable[[float], float], x0: float,
    df: Callable[[Callable[[float], float], float, float, float], float],
    d2f: Callable[[Callable[[float], float], float, float, float], float],
    tol: float = 1e-8, k_max: int = 1000, r: float = 1.0) -> float|None:
    """Returns the root of a function around an initial guess
    using the Ralston Rabinowitz method.

    Args:
        f: an algebraic function
        df: the derivative of the function calculated with finite differences
        d2f: the second derivative of the function calculated with finite differences
        x0: the initial guess
        tol: the numerical threshold for convergence
        k_max: the maximum number of iterations

    Returns:
        A root (float) if it exists, otherwise None
    """

    # return scipy.optimize.newton(f, x0=x0, fprime=lambda x: df(f,x),
    #     fprime2=lambda x: d2f(f,x)) # Halley

    x = x0
    for k in range(1, k_max+1):

        x_old = x

        f_x = f(x)

        df_x = df(f, x)

        d2f_x = d2f(f, x)

        try:
            x = x_old - r*f_x*df_x/( (df_x**2)-f_x*d2f_x )
        except ZeroDivisionError:
            return None

        # print(f"k: {k}, x: {x}, f(x): {f_x}")

        err = [abs(f_x), abs( (x-x_old)/(x+1e-12) )]

        if all(e < tol for e in err):
            print(f'k = {k}. Errors (f, rel): {[f"{e:.4e}" for e in err]}')
            return x

    return None
