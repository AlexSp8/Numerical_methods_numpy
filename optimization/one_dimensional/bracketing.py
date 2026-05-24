
from typing import Callable
import scipy

from optimization.one_dimensional import hybrid, open_methods

def golden_section_search(f: Callable[[float], float],
    a: float, b: float, tol: float = 1e-8, k_max: int = 1000) -> float:
    """Returns the minimum of a function in an interval
    using the golden-section search method.

    Args:
        f: the function to optimize
        a: the lower limit of the interval
        b: the upper limit of the interval
        tol (optional): threshold for convergence
        k_max (optional): maximum iterations
    Returns:
        The point of the function minimum."""

    # x = scipy.optimize.golden(f, brack=(a, b), tol=tol)
    # return float(x)

    # r = scipy.optimize.minimize_scalar(f, bounds=(a, b), method='bounded')
    # return float(r.x)

    golden_ratio = ((5.0**0.5)-1.0)/2.0

    xl, xu = a, b

    d = golden_ratio*(xu-xl)

    x1, x2 = xu-d, xl+d
    f1, f2 = f(x1), f(x2)

    for k in range(1, k_max+1):

        if f2 < f1:
            # x2 is better, remove xl
            xl = x1
            x1, f1 = x2, f2
            x2 = xl + golden_ratio*(xu-xl)
            f2 = f(x2)
            x_opt = x2
        else:
            # x1 is better, remove xu
            xu = x2
            x2, f2 = x1, f1
            x1 = xu - golden_ratio*(xu-xl)
            f1 = f(x1)
            x_opt = x1

        x_int = xu - xl
        err = abs( x_int/(x_opt+1e-12) )
        if err < tol:
            print(f'k = {k}. Rel. Error: {f"{err:.4e}"}')
            break

    return x_opt

def multi_bracketing(f: Callable[[float], float],
    a: float, b: float, n: int, method: str,
    df: Callable[[Callable[[float], float], float, float, float], float] = None,
    d2f: Callable[[Callable[[float], float], float, float, float], float] = None) -> list[float]:
    """Returns a list of minimum values of a function in an interval
    using a bracketing optimization method.

    Args:
        f: the function
        a: the lower limit of the total interval
        b: the upper limit of the total interval
        n: the number of sub-intervals
        method: the optimization method to be used
        df (optional): the function to evaluate the derivative of f
        d2f (optional): the function to evaluate the 2nd derivative of f
    Returns:
        The list of points of the function minima."""

    if n <= 0:
        raise ValueError('n <= 0 in multi-bracketing!')

    extremes = []
    dx = (b-a)/n
    for i in range(n):

        a_int = a + i*dx
        b_int = a_int + dx

        if method == 'golden-section':
            x = golden_section_search(f, a_int, b_int)
        elif method == 'brent':
            x = hybrid.brent(f, a_int, b_int)
        elif method == 'parabolic_interpolation':
            x = open_methods.parabolic_interpolation(f, a_int, b_int)
        elif method == 'secant':
            x0 = (a_int+b_int)/2.0
            x = open_methods.secant(f, x0, df)
        elif method == 'newton':
            x0 = (a_int+b_int)/2.0
            x = open_methods.newton(f, x0, df, d2f)
        else:
            raise ValueError('Invalid bracketing method!')

        extremes.append(x)

    return extremes
