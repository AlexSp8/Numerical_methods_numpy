
from typing import Callable

import numpy as np

def richardson_extrapolation(df: Callable[[Callable[[float], float], float, float], float],
    f: Callable[[float], float], x: float, h: float = 1e-8,
    k_max: int = 10, eps: float = 1e-8) -> float:
    """Returns the value of the integral of function, f, in [a,b]
    using the Romberg's formula."""

    d = np.zeros((k_max+1,k_max+1))

    d[0,0] = df(f, x, h)

    for k in range(1, k_max+1):

        d[k,0] = df(f, x, h/(2**k))

        for j in range(1,k+1):
            w = 4**j
            d[k,j] = ( w*d[k,j-1] - d[k-1,j-1] )/(w-1)

        err = abs( (d[k,k]-d[k-1,k-1])/d[k,k] )
        if err < eps:
            # print(f"Richardson converged at k = {k}")
            break

    return d[k][k]
