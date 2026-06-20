
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

from utilities import matrix_operations
from linear_systems.linear_solver import LinearSolver

class DirectSolver(LinearSolver):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

class SparseSolver(DirectSolver):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def solve(self, A: np.ndarray[tuple[int, int], np.dtype[np.float64]],
        b: np.ndarray[tuple[int], np.dtype[np.float64]]
        ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        A_sparse = csr_matrix(A)
        return spsolve(A_sparse, b)


class LUSolver(DirectSolver):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def solve(self, A: np.ndarray[tuple[int, int], np.dtype[np.float64]],
        b: np.ndarray[tuple[int], np.dtype[np.float64]]
        ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        """Returns the solution vector of a linear system of algebraic equations
        (Ax = b) with LU decomposition.
        First, A = LU (decomposition, including partial pivoting) is performed.
        Then, Ld = b (forward substitution). Finally, Ux = d (back substitution).

        Args:
            A: the matrix of the linear system
            b: the right-hand side vector of the linear system
        Returns:
            The solution vector of the linear system, x."""

        # A_lu, rows_order = scipy.linalg.lu_factor(A)
        # n = A.shape[0]
        # U, L = np.triu(A_lu), np.tril(A_lu,-1) + np.eye(n)
        # return scipy.linalg.lu_solve((A_lu, rows_order), b)

        A_lu, rows_order = matrix_operations.lu_decomposition(A)

        # print(matrix_operations.determinant_upper_diagonal(A_lu))
        # print(matrix_operations.log_determinant_upper_diagonal(A_lu))

        b_lu = b[rows_order]

        d = matrix_operations.forward_substitution(A_lu, b_lu)

        return matrix_operations.back_substitution(A_lu, d)


class GaussSolver(DirectSolver):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def solve(self, A: np.ndarray[tuple[int, int], np.dtype[np.float64]],
        b: np.ndarray[tuple[int], np.dtype[np.float64]]
        ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        """Returns the solution vector of a linear system of algebraic equations
        (Ax = b) with Gauss elimination.

        Args:
            A: the matrix of the linear system
            b: the right-hand side vector of the linear system
        Returns:
            The solution vector of the linear system, x.
        """

        A_aug = self.forward_elimination(A, b)

        n = A_aug.shape[0]

        b_up = A_aug[:,n]

        return matrix_operations.back_substitution(A_aug[:,:-1], b_up)

    def forward_elimination(self, A: np.ndarray[tuple[int, int], np.dtype[np.float64]],
        b: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        """Returns the augmented matrix that result from
        forward elimination of a square matrix A.

        Args:
            A: the matrix of the linear system
            b: the right-hand side vector of the linear system
        Returns:
            The augmented matrix that result after forward elimination."""

        n = A.shape[0]

        n_swaps = 0

        A_aug = np.column_stack((A,b))
        # A_aug = np.zeros((n, n+1))
        # A_aug[:n,:n] = A.copy()
        # A_aug[:,n] = b.copy()

        for k in range(n):

            i_max = matrix_operations.pivot_row(A_aug, k)
            if i_max != k:
                A_aug[[k,i_max]] = A_aug[[i_max,k]]
                n_swaps += 1

            for i in range(k+1,n):
                f = A_aug[i,k]/A_aug[k,k]
                A_aug[i,k:] -= f*A_aug[k,k:]

        # print(matrix_operations.determinant_upper_diagonal(A_aug, n_swaps))
        # print(matrix_operations.log_determinant_upper_diagonal(A_aug, n_swaps))

        return A_aug


class CramerSolver(DirectSolver):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def solve(self, A: np.ndarray[tuple[int, int], np.dtype[np.float64]],
        b: np.ndarray[tuple[int], np.dtype[np.float64]]
        ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        """Returns the solution vector of a linear system of algebraic equations
        (Ax = b) with Cramer's rule.

        Args:
            A: the matrix of the linear system
            b: the right-hand side vector of the linear system
        Returns:
            The solution vector of the linear system, x."""

        detA = np.linalg.det(A)
        if abs(detA) < 1e-12:
            raise ValueError("Matrix is singular! Cannot perform Cramer's rule.")

        n = A.shape[0]

        x = np.zeros(n)

        for i in range(n):
            Ai = matrix_operations.replace_column(A, b, i)
            detAi = np.linalg.det(Ai)
            x[i] = detAi/detA

        return x


def thomas_tri_diagonal(l:np.ndarray[tuple[int], np.dtype[np.float64]],
    d: np.ndarray[tuple[int], np.dtype[np.float64]],
    u: np.ndarray[tuple[int], np.dtype[np.float64]],
    b: np.ndarray[tuple[int], np.dtype[np.float64]]
    ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """Returns the solution vector of a linear system of algebraic equations
    (Ax = b) with the Thomas algorithm for tri-diagonal matrices.

    Args:
        l: the sub-diagonal elements of the matrix
        d: the diagonal elements of the matrix
        u: the upper-diagonal elements of the matrix
        b: the right-hand side vector of the linear system
    Returns:
        The solution vector of the linear system, x."""

    # udl = np.array([
    #     np.insert(u, 0, 0),  # Upper: with 0 at the start
    #     d,                   # Main
    #     np.append(l, 0)      # Lower: with 0 at the end
    # ])
    # x = scipy.linalg.solve_banded((1,1), udl, b)
    # return x

    n = len(d)

    d_lu = d.copy()

    b_lu = b.copy()

    for k in range(1,n):
        # Decomposition
        f = l[k-1]/d_lu[k-1]
        d_lu[k] -= f*u[k-1]
        # Forward substitution
        b_lu[k] -= f*b_lu[k-1]

    # Back substitution
    x = np.zeros(n)
    x[n-1] = b_lu[n-1]/d_lu[n-1]
    for k in range(n-2, -1, -1):
        x[k] = ( b_lu[k] - u[k]*x[k+1] )/d_lu[k]

    return x
