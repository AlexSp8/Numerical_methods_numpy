
from typing import Callable
import numpy as np

class LineSearch:
    """A utility class to transform a multi-variable function into a 1D function.
    This class 'freezes' a starting point and a direction vector, allowing
    1D optimization algorithms (like Golden Section or Brent's method) to
    search for an optimal step size (lambda) along that specific vector.

    Attributes:
        f: The multi-variable function to be optimized.
        x: The current position in the n-dimensional space.
        d: The direction vector along which to search."""

    def __init__(self, f: Callable[[np.ndarray[tuple[int], np.dtype[np.float64]]], float]):
        """Initializes the LineSearch wrapper with an function.

        Args:
            f: A function that takes a NumPy array and returns a float."""

        self.f = f
        self.x = np.array([])
        self.d = np.array([])

    def update(self, x: np.ndarray[tuple[int], np.dtype[np.float64]],
        d: np.ndarray[tuple[int], np.dtype[np.float64]]) -> None:
        """Updates the current position and direction for the next line search.

        Args:
            x: The current starting point (origin) for the search.
            d: The direction vector to explore."""

        self.x = x
        self.d = d

    def f_line(self, h: float) -> float:
        """Evaluates the objective function at a distance, h, along direction, d.

        Args:
            h: The step size (scalar) to move from the starting point.

        Returns:
            float: The value of the objective function at x + h*d."""

        return self.f(self.x + h*self.d)

    def h_interval(self, k_max = 100) -> tuple[float, float]:
        """Returns the interval for 1D Line Search.

        Args:
            k_max (optional): maximum iterations
        Returns:
            The lower and upper limit of the interval to search."""

        h_min, h_max = 0.0, 0.1

        f_min = self.f_line(h_min)

        for k in range(k_max):

            f_max = self.f_line(h_max)

            if f_min > f_max:
                h_min = h_max
                f_min = f_max
                h_max *= 10.0
            else:
                # print('h_interval iter k =', k)
                break

        return h_min, h_max
