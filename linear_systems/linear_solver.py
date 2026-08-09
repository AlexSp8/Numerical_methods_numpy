
from abc import ABC, abstractmethod
import numpy as np
from numpy.typing import NDArray

class LinearSolver(ABC):

    def __init__(self):
        pass

    @abstractmethod
    def solve(self):
        pass

    @staticmethod
    def get_residual(A: NDArray[np.float64],
        b: NDArray[np.float64],
        x: NDArray[np.float64]
        ) -> NDArray[np.float64]:
        """Returns the residual vector of a linear system, r = b - Ax.

        Args:
            A: the matrix of the linear system
            b: the right-hand side vector of the linear system
            x: the solution vector approximation of the linear system

        Returns:
            The residual vector of the linear system."""

        return b - np.dot(A, x)
