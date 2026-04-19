
from typing import Callable, List
import scipy

def golden_section_search(f: Callable[[float], float],
    a: float, b: float, mode: str = 'min',
    eps: float = 1e-8, k_max: int = 1000) -> float:
    """Returns the extreme of a function f(x) in an interval [a, b]
    using the golden-section search method"""

    if mode == 'min':
        s = +1.0
    else:
        s = -1.0

    # g = lambda x: s*f(x)
    # x = scipy.optimize.golden(g, brack=(a, b), tol=eps)
    # return float(x)

    # r = scipy.optimize.minimize_scalar(g, bounds=(a, b), method='bounded')
    # return float(r.x)

    golden_ratio = ((5.0**0.5)-1.0)/2.0

    xl, xu = a, b

    d = golden_ratio*(xu-xl)

    x1, x2 = xl+d, xu-d
    f1, f2 = f(x1), f(x2)

    for k in range(1, k_max+1):

        if s*f1 < s*f2:
            xl = x2
            x2, f2 = x1, f1
            x1 = xl + golden_ratio*(xu-xl)
            f1 = f(x1)
            x_opt = x1
            f_opt = f1
        else:
            xu = x1
            x1, f1 = x2, f2
            x2 = xu - golden_ratio*(xu-xl)
            f2 = f(x2)
            x_opt = x2
            f_opt = f2

        x_int = xu - xl
        err = abs( x_int/(x_opt+1e-12) )
        if err < eps:
            print(f'k = {k}. Rel. Error: {f"{err:.4e}"}')
            break

    return x_opt

def multi_bracketing(f: Callable[[float], float],
    a: float, b: float, n: int, method: str,
    mode: str = 'min') -> List[float]:
    """Returns a list of extremes of a function f(x) in an interval [a, b]
    using a bracketing method across n intervals"""

    if n <= 0:
        raise ValueError('n <= 0 in multi-bracketing!')

    extremes = []
    dx = (b-a)/n
    for i in range(n):

        a_int = a + i*dx
        b_int = a_int + dx

        if method == 'golden-section':
            x = golden_section_search(f, a_int, b_int, mode)
        else:
            raise ValueError('Invalid bracketing method!')

        extremes.append(x)

    return extremes
