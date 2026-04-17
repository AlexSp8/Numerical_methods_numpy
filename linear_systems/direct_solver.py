
import numpy as np
import numpy.typing as npt

from utilities import matrix_operations
from linear_systems.linear_solver import LinearSolver

class DirectSolver(LinearSolver):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class LUSolver(DirectSolver):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def solve(self, A: npt.NDArray[np.float64], b: npt.NDArray[np.float64]
        ) -> npt.NDArray[np.float64]:
        """Returns the solution vector, x, of a linear system of algebraic equations (Ax = b)
        with LU decomposition. First, A = LU (decomposition, including partial pivoting).
        Then, Ld = b (forward substitution). Finally, Ux = d (back substitution)."""

        # A_lu, rows_order = scipy.linalg.lu_factor(A)
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

    def solve(self, A: npt.NDArray[np.float64], b: npt.NDArray[np.float64]
        ) -> npt.NDArray[np.float64]:
        """Solves a linear system of algebraic equations with Gauss elimination.
        Partial pivoting is performed during forward elimination."""

        A_aug, n_swaps = self.forward_elimination(A, b)

        n = A_aug.shape[0]

        # print(matrix_operations.determinant_upper_diagonal(A_aug, n_swaps))
        # print(matrix_operations.log_determinant_upper_diagonal(A_aug, n_swaps))

        b_up = A_aug[:,n]

        return matrix_operations.back_substitution(A_aug[:,:-1], b_up)

    def forward_elimination(self, A: npt.NDArray[np.float64], b: npt.NDArray[np.float64]
        ) -> npt.NDArray[np.float64]:
        """Returns the augmented upper diagonal matrix that result from
        forward elimination of unknowns of a square matrix A"""

        n = A.shape[0]
        n_swaps = 0
        A_aug = np.column_stack((A,b))
        for k in range(n):

            i_max = matrix_operations.partial_pivot(A_aug, k)
            if i_max != k:
                A_aug[[k,i_max]] = A_aug[[i_max,k]]
                n_swaps += 1

            for i in range(k+1,n):
                f = A_aug[i,k]/A_aug[k,k]
                A_aug[i,k:] -= f*A_aug[k,k:]

        return A_aug, n_swaps


class CramerSolver(DirectSolver):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def solve(self, A: npt.NDArray[np.float64], b: npt.NDArray[np.float64]
        ) -> npt.NDArray[np.float64]:
        """Solves a linear system of algebraic equations with Cramer's rule"""

        detA = np.linalg.det(A)
        if np.isclose(detA, 0):
            raise ValueError("Matrix is singular! Cannot perform Cramer's rule.")

        n = A.shape[0]
        x = np.zeros(n)
        for i in range(n):
            Ai = matrix_operations.replace_column(A, b, i)
            detAi = np.linalg.det(Ai)
            x[i] = detAi/detA

        return x


def thomas_tri_diagonal(l: npt.NDArray[np.float64], d: npt.NDArray[np.float64],
    u: npt.NDArray[np.float64], b: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Returns the solution, x, to a system of linear algebraic equations,
    Ax = b, using the Thomas algorithm for tri-diagonal matrices.
    Matrix elements are stored in l (sub-diagonal), d (diagonal),
    and u (upper-diagonal). """

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
