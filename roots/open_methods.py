
from typing import Callable

from differentiation import forward_fd as ffd

def fixed_point(f: Callable[[float], float], x0: float,
    eps: float = 1e-8, k_max: int = 1000) -> float:
    """Returns the root, f(x) = 0, of a function f(x)
    around an initial guess, x0, using the fixed-point iteration method"""

    x = x0
    for k in range(1, k_max+1):

        x_old = x

        f_x = f(x)

        x = f_x + x

        # print(f"k: {k}, x: {x}, f(x): {f_x}")

        err = [abs(f_x), abs(x - x_old), abs( (x-x_old)/(x+eps) )]

        if all(e < eps for e in err):
            print(f'k = {k}. Errors (f, abs, rel): {[f"{e:.4e}" for e in err]}')
            break
    return x

def newton_raphson(f: Callable[[float], float], x0: float,
    eps: float = 1e-8, k_max: int = 1000, r: float = 1.0) -> float:
    """Returns the root, f(x) = 0, of a function f(x)
    around an initial guess, x0, using the Newton-Raphson method"""

    x = x0
    for k in range(1, k_max+1):

        x_old = x

        f_x = f(x)
        df_x = ffd.df_h(f, x)

        x = x_old - r*f_x/(df_x+eps)

        # print(f"k: {k}, x: {x}, f(x): {f_x}")

        err = [abs(f_x), abs(x - x_old), abs( (x-x_old)/(x+eps) )]

        if all(e < eps for e in err):
            print(f'k = {k}. Errors (f, abs, rel): {[f"{e:.4e}" for e in err]}')
            break

    return x

def ralston_rabinowitz(f: Callable[[float], float], x0: float,
    eps: float = 1e-8, k_max: int = 1000, r: float = 1.0) -> float:
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

        err = [abs(f_x), abs(x - x_old), abs( (x-x_old)/(x+eps) )]

        if all(e < eps for e in err):
            print(f'k = {k}. Errors (f, abs, rel): {[f"{e:.4e}" for e in err]}')
            break

    return x
