
from typing import Callable
import numpy as np
import numpy.typing as npt

from optimization.one_dimensional import hybrid
from optimization.unconstrained.line_search import LineSearch

def powell(f: Callable[[npt.NDarray[np.float64]], float],
    x0: npt.NDarray[np.float64] = None,
    d: npt.NDarray[np.float64] = None, mode: str = 'min',
    eps: float = 1e-8, k_max = 1000) -> npt.NDarray[np.float64]:
    """Returns the point x(x1, x2, ..., xn) where the extreme of a
    multi-variable function, f, is found using Powell's method.
    A starting guess point, x0, and the initial direction vectors, d, should be given."""

    n = len(x0)

    x = x0.copy()

    if d is not None:
        d_vectors = d.copy()
    else:
        d_vectors = np.eye(n)

    l_min, l_max = -1.0, 1.0

    line_search = LineSearch(f)

    for k in range(1, k_max+1):

        x_start = x.copy()
        max_diff, max_diff_dir = 0.0, 0

        # Line search along i direction
        for i in range(n):

            line_search.update(x, d_vectors[i])

            # Find h_opt (step size along i direction)
            h_opt = hybrid.brent(line_search.f_line, l_min, l_max, mode)

            x_old = x.copy()

            x = x + h_opt*d_vectors[i]

            # Track the direction in which f changes the most
            diff = abs(f(x) - f(x_old))
            if diff > max_diff:
                max_diff, max_diff_dir = diff, i

        f_x = f(x)
        if abs( (f_x - f(x_start))/(f_x+eps) ) < eps:
            # print('k =', k)
            break

        # New (normalized) direction dn = (Pn - P0)
        d_new = x - x_start
        norm = np.linalg.norm(d_new)
        if norm > 1e-12:
            d_vectors = np.delete(d_vectors, max_diff_dir, axis=0)
            d_vectors = np.vstack([d_vectors, d_new / norm])

    return x
