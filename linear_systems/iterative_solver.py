
from typing import Tuple

import numpy as np
import numpy.typing as npt

from linear_systems.linear_solver import LinearSolver

class IterativeSolver(LinearSolver):

    def __init__(self, *args, k_max: int = 1000, tol: float = 1e-8,
        x0: npt.NDArray[np.float64] = None, **kwargs):

        super().__init__(*args, **kwargs)
        self.k_max = k_max
        self.tol = tol
        self.x0 = x0

    def set_initial_guess(self, A: npt.NDArray[np.float64],
        b: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:

        if self.x0 is not None:
            return self.x0.copy()

        n = A.shape[0]
        diagA = np.diag(A)

        if np.any(diagA == 0):
            return np.ones(n)
        else:
            return b/diagA

    def get_norms(self, dx: npt.NDArray[np.float64], res: npt.NDArray[np.float64]
        ) -> Tuple[float, float]:
        """Returns the correction and residual norms."""
        cor_norm = np.linalg.norm(dx)
        res_norm = np.linalg.norm(res)
        return cor_norm, res_norm

    def is_converged(self, cor_norm: float, res_norm: float, mode: str = 'and') -> bool:
        """Checks the correction and residual norm convergence criteria."""
        if mode == 'and':
            return (res_norm < self.tol) and (cor_norm < self.tol)
        else:
            return (res_norm < self.tol) or (cor_norm < self.tol)


class JacobiSolver(IterativeSolver):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

    def solve(self, A: npt.NDArray[np.float64], b: npt.NDArray[np.float64],
        output: bool = False) -> npt.NDArray[np.float64]:
        """Returns the solution, x, to a system of linear algebraic equations,
        Ax = b, using the Jacobi iterative method."""

        diagA = np.diag(A)
        if any(diagA == 0):
            raise ValueError("Diagonal contains zero")

        x = self.set_initial_guess(A, b)

        for k in range(1, self.k_max+1):

            x_old = x.copy()

            res = self.get_residual(A, b, x_old)

            x = x_old + res/diagA

            cor_norm, res_norm = self.get_norms(x-x_old, res)

            if self.is_converged(cor_norm, res_norm, 'and'):
                break

        if output:
            print(f'k = {k}, Res Norm = {res_norm:.4e}, Cor Norm = {cor_norm:.4e}')

        return x


class SORSolver(IterativeSolver):

    def __init__(self, *args, w=1.5, **kwargs):

        super().__init__(*args, **kwargs)
        self.w = w

    def solve(self, A: npt.NDArray[np.float64], b: npt.NDArray[np.float64],
        output: bool = False) -> npt.NDArray[np.float64]:
        """Returns the solution, x, to a system of linear algebraic equations,
        Ax = b, using the SOR iterative method. The relaxation can be optionally
        given. For ω = 1, the Gauss-Seidel method is used."""

        n = A.shape[0]

        diagA = np.diag(A)
        if any(diagA == 0):
            raise ValueError("Diagonal contains zero")

        x = self.set_initial_guess(A, b)

        for k in range(1, self.k_max+1):

            x_old = x.copy()

            res = np.zeros(n)
            for i in range(n):
                res[i] = b[i] - np.dot(A[i], x)
                x[i] += self.w*res[i]/diagA[i]

            res = self.get_residual(A, b, x)
            cor_norm, res_norm = self.get_norms(x-x_old, res)

            if self.is_converged(cor_norm, res_norm, 'and'):
                break

        if output:
            print(f'k = {k}, Res Norm = {res_norm:.4e}, Cor Norm = {cor_norm:.4e}')

        return x


class SDSolver(IterativeSolver):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def solve(self, A: npt.NDArray[np.float64], b: npt.NDArray[np.float64],
        output: bool = False) -> npt.NDArray[np.float64]:
        """Returns the solution, x, to a system of linear algebraic equations,
        Ax = b, using the Steepest Descent iterative method."""

        x = self.set_initial_guess(A, b)

        r = self.get_residual(A, b, x)

        for k in range(1, self.k_max+1):

            rr = np.dot(r, r)

            Ar = np.dot(A, r)

            rAr = np.dot(r, Ar)

            alpha = rr/(rAr+1e-12)

            x_old = x.copy()

            x = x + alpha*r
            r = r - alpha*Ar

            # res_norm = np.sqrt(rr)
            cor_norm, res_norm = self.get_norms(x-x_old, r)

            if self.is_converged(cor_norm, res_norm, 'or'):
                break

        if output:
            print(f'k = {k}, Res Norm = {res_norm:.4e}, Cor Norm = {cor_norm:.4e}')

        return x


class CGSolver(IterativeSolver):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def solve(self, A: npt.NDArray[np.float64], b: npt.NDArray[np.float64],
        output: bool = False) -> npt.NDArray[np.float64]:
        """Returns the solution, x, to a system of linear algebraic equations,
        Ax = b, using the Conjugate Gradient iterative method."""

        # x, info = scipy.sparse.linalg.cg(A, b, tol=self.tol)
        # return x

        x = self.set_initial_guess(A, b)

        r = self.get_residual(A, b, x)
        rr_old = np.dot(r, r)

        p = r.copy()
        for k in range(1, self.k_max+1):

            Ap = np.dot(A, p)
            pAp = np.dot(p, Ap)
            alpha = rr_old/(pAp+1e-12)

            x_old = x.copy()

            x = x + alpha*p
            r = r - alpha*Ap

            rr = np.dot(r,r)

            # res_norm = np.sqrt(rr)
            cor_norm, res_norm = self.get_norms(x-x_old, r)

            if self.is_converged(cor_norm, res_norm, 'or'):
                break

            beta = rr/rr_old
            p = r + beta*p

            rr_old = rr

        if output:
            print(f'k = {k}, Res Norm = {res_norm:.4e}, Cor Norm = {cor_norm:.4e}')

        return x
