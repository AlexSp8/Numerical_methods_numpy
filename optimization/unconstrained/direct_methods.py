
from typing import Callable
import numpy as np

from optimization.one_dimensional import hybrid
from optimization.unconstrained.line_search import LineSearch

def powell(f: Callable[[np.ndarray[tuple[int], np.dtype[np.float64]]], float],
    x0: np.ndarray[tuple[int], np.dtype[np.float64]],
    d: np.ndarray[tuple[int], np.dtype[np.float64]] = None,
    tol: float = 1e-8, k_max = 100) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """Returns the point of minimum of a multi-variable function using Powell's method.

    Args:
        f: the function to be optimized
        x0: starting guess
        d (optional): starting direction vectors
        tol (optional): threshold for convergence
        k_max (optional): maximum iterations
    Returns:
        The point of minimum of f."""

    n = len(x0)

    x = x0.copy()

    if d is not None:
        d_vectors = d.copy()
    else:
        d_vectors = np.eye(n)

    line_search = LineSearch(f)

    for k in range(1, k_max+1):

        x_start = x.copy()
        max_diff, max_diff_dir = 0.0, 0

        # Line search along i direction
        for i in range(n):

            line_search.update(x, d_vectors[i])

            h_min, h_max = line_search.h_interval(k_max=10)

            # Find h_opt (step size along i direction)
            h_opt = hybrid.brent(line_search.f_line, h_min, h_max)

            x_old = x.copy()

            x = x + h_opt*d_vectors[i]

            # Track the direction in which f changes the most
            diff = np.abs(f(x) - f(x_old))
            if diff > max_diff:
                max_diff, max_diff_dir = diff, i

        f_x = f(x)
        if np.abs( (f_x - f(x_start))/(f_x+tol) ) < tol:
            # print('k =', k)
            break

        # New (normalized) direction dn = (Pn - P0)
        d_new = x - x_start
        norm = np.linalg.norm(d_new)
        if norm > 1e-12:
            d_vectors = np.delete(d_vectors, max_diff_dir, axis=0)
            d_vectors = np.vstack([d_vectors, d_new/norm])

    return x
