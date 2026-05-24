
from abc import ABC, abstractmethod
import numpy as np

class LinearSolver(ABC):

    def __init__(self):
        pass

    @abstractmethod
    def solve(self):
        pass

    @staticmethod
    def get_residual(A: np.ndarray[tuple[int, int], np.dtype[np.float64]],
        b: np.ndarray[tuple[int], np.dtype[np.float64]],
        x: np.ndarray[tuple[int], np.dtype[np.float64]]
        ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        """Returns the residual vector of a linear system, r = b - Ax.

        Args:
            A: the matrix of the linear system
            b: the right-hand side vector of the linear system
            x: the solution vector approximation of the linear system

        Returns:
            The residual vector of the linear system."""

        return b - np.dot(A, x)
