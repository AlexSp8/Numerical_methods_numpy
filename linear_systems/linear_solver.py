
from abc import ABC, abstractmethod
import numpy as np
import numpy.typing as npt

class LinearSolver(ABC):

    def __init__(self):
        pass

    @abstractmethod
    def solve(self):
        pass

    @staticmethod
    def get_residual(A: npt.NDArray[np.float64], b: npt.NDArray[np.float64],
        x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return b - np.dot(A, x)
