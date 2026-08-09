
import numpy as np
from numpy.typing import NDArray

from linear_systems.linear_solver import LinearSolver

class IterativeSolver(LinearSolver):

    def __init__(self, *args, k_max: int = 1000, tol: float = 1e-8,
        x0: NDArray[np.float64] = None, **kwargs):

        super().__init__(*args, **kwargs)
        self.k_max = k_max
        self.tol = tol
        self.x0 = x0

    def set_initial_guess(self, A: NDArray[np.float64], b: NDArray[np.float64]
        ) -> NDArray[np.float64]:
        """Sets the initial guess for the solution vector.

        Args:
            A: the matrix of the linear system
            b: the right-hand side vector of the linear system
        Returns:
            The initial guess for the solution vector, x0."""

        if self.x0 is not None:
            return self.x0.copy()

        diagA = np.diag(A)

        if np.any(diagA == 0):
            n = A.shape[0]
            return np.ones(n)
        else:
            return b/diagA

    def get_norms(self, dx: NDArray[np.float64], res: NDArray[np.float64]
        ) -> tuple[float, float]:
        """Returns the norms of the correction and residual vectors.

        Args:
            dx: the correction vector
            res: the residual vector
        Returns:
            The correction and residual norms."""

        cor_norm = np.linalg.norm(dx)
        res_norm = np.linalg.norm(res)

        return cor_norm, res_norm


class JacobiSolver(IterativeSolver):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

    def solve(self, A: NDArray[np.float64], b: NDArray[np.float64],
        output: bool = False) -> NDArray[np.float64]:
        """Returns the solution vector of a linear system of algebraic equations
        (Ax = b) with the Jacobi method.

        Args:
            A: the matrix of the linear system
            b: the right-hand side vector of the linear system
            output: flag for iteration info output
        Returns:
            The solution vector of the linear system, x."""

        diagA = np.diag(A)
        if any(diagA == 0):
            raise ValueError("Diagonal contains zero")

        x = self.set_initial_guess(A, b)

        for k in range(1, self.k_max+1):

            x_old = x.copy()

            res = self.get_residual(A, b, x)

            x = x_old + res/diagA

            cor_norm, res_norm = self.get_norms(x-x_old, res)

            if output:
                print(f'k = {k}, Res Norm = {res_norm:.4e}, Cor Norm = {cor_norm:.4e}')

            if (res_norm < self.tol) and (cor_norm < self.tol):
                break

        return x


class SORSolver(IterativeSolver):

    def __init__(self, *args, w: float = 1.5, **kwargs):

        super().__init__(*args, **kwargs)
        self.w = w

    def solve(self, A: NDArray[np.float64], b: NDArray[np.float64],
        output: bool = False) -> NDArray[np.float64]:
        """Returns the solution vector of a linear system of algebraic equations
        (Ax = b) with the SOR method.

        Args:
            A: the matrix of the linear system
            b: the right-hand side vector of the linear system
            output: flag for iteration info output
        Returns:
            The solution vector of the linear system, x."""

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

            if output:
                print(f'k = {k}, Res Norm = {res_norm:.4e}, Cor Norm = {cor_norm:.4e}')

            if (res_norm < self.tol) and (cor_norm < self.tol):
                break

        return x


class SDSolver(IterativeSolver):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def solve(self, A: NDArray[np.float64], b: NDArray[np.float64],
        output: bool = False) -> NDArray[np.float64]:
        """Returns the solution vector of a linear system of algebraic equations
        (Ax = b) with the Steepest Descent method.

        Args:
            A: the matrix of the linear system
            b: the right-hand side vector of the linear system
            output: flag for iteration info output
        Returns:
            The solution vector of the linear system, x."""

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

            if output:
                print(f'k = {k}, Res Norm = {res_norm:.4e}, Cor Norm = {cor_norm:.4e}')

            if (res_norm < self.tol) or (cor_norm < self.tol):
                break

        return x


class CGSolver(IterativeSolver):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def solve(self, A: NDArray[np.float64], b: NDArray[np.float64],
        output: bool = False) -> NDArray[np.float64]:
        """Returns the solution vector of a linear system of algebraic equations
        (Ax = b) with the Conjugate Gradient method.

        Args:
            A: the matrix of the linear system
            b: the right-hand side vector of the linear system
            output: flag for iteration info output
        Returns:
            The solution vector of the linear system, x."""

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

            if output:
                print(f'k = {k}, Res Norm = {res_norm:.4e}, Cor Norm = {cor_norm:.4e}')

            if (res_norm < self.tol) or (cor_norm < self.tol):
                break

            beta = rr/rr_old
            p = r + beta*p

            rr_old = rr

        return x
