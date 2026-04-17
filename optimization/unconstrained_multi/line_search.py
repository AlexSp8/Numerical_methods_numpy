
from typing import Callable
import numpy as np
import numpy.typing as npt

class LineSearch:
    """
    A utility class to transform a multivariable function into a 1D function.

    This class 'freezes' a starting point and a direction vector, allowing
    1D optimization algorithms (like Golden Section or Brent's method) to
    search for an optimal step size (lambda) along that specific vector.

    Attributes:
        f (Callable): The multivariable objective function to be minimized.
        x (npt.NDArray[np.float64]): The current position in the n-dimensional space.
        d (npt.NDArray[np.float64]): The direction vector along which to search.
    """

    def __init__(self, f: Callable[[npt.NDArray[np.float64]], float]):
        """
        Initializes the LineSearch wrapper with an objective function.

        Args:
            f: A function that takes a NumPy array and returns a float.
        """
        self.f = f
        self.x: npt.NDArray[np.float64] = np.array([])
        self.d: npt.NDArray[np.float64] = np.array([])

    def update(self, x: npt.NDArray[np.float64], d: npt.NDArray[np.float64]) -> None:
        """
        Updates the current position and direction for the next line search.

        Args:
            x: The current starting point (origin) for the search.
            d: The direction vector to explore.
        """
        self.x = x
        self.d = d

    def f_line(self, h: float) -> float:
        """
        Evaluates the objective function at a distance 'h' along direction 'd'.

        This method acts as the g(h) function required by 1D optimizers.

        Args:
            h: The step size (scalar) to move from the starting point.

        Returns:
            float: The value of the objective function at x + h*d.
        """
        return self.f(self.x + h*self.d)
