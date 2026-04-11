
from typing import List, Tuple
import math
import numpy as np
import numpy.typing as npt
import scipy

from linear_systems import direct_solvers

def addition(A: List[List[float]],
    B: List[List[float]], v: float = 1.0) -> List[List[float]]:
    """Returns the resulting matrix, C(n x m) from the
    generalized addition of two matrices,
    A(n x m), B(n x m): C = A + v*B"""

    n, m = len(A), len(A[0])

    C = [ [0.0 for _ in range(m)] for _ in range(n) ]
    for i in range(n):
        for j in range(m):
            C[i][j] = A[i][j] + v*B[i][j]

    return C

def scalar_multiplication(A: List[List[float]],
    v: float) -> List[List[float]]:
    """Returns the resulting matrix, B(n x m), from the
    multiplication of a matrix A(n x m) with a scalar, v: B = v*A"""

    n, m = len(A), len(A[0])

    B = [ [0.0 for _ in range(m)] for _ in range(n) ]
    for i in range(n):
        for j in range(m):
            B[i][j] = v*A[i][j]

    return B

def element_multiplication(A: List[List[float]],
    B: List[List[float]]) -> List[List[float]]:
    """Returns the resulting matrix, C(n x m), from the element-wise
    multiplication of two matrices A(n x m), B(n x m): C = A*B"""

    n, m = len(A), len(A[0])

    C = [ [0.0 for _ in range(m)] for _ in range(n) ]
    for i in range(n):
        for j in range(m):
            C[i][j] = A[i][j]*B[i][j]

    return C

def multiplication(A: List[List[float]],
    B: List[List[float]]) -> List[List[float]]:
    """Returns the resulting matrix, C(n x l), from the
    multiplication of two matrices A(n x m), B(m x l): C = A@B"""

    n, m, l = len(A), len(A[0]), len(B[0])

    C = [ [0.0 for _ in range(l)] for _ in range(n) ]
    for i in range(n):
        for j in range(l):
            s = 0
            for k in range(m):
                s += A[i][k]*B[k][j]
            C[i][j] = s

    return C

def matrix_vector_multiplication(A: List[List[float]],
    x: List[float]) -> List[float]:
    """Multiplies an (m x n) matrix, A, by an (n) vector, x."""

    m = len(A)
    n = len(A[0])
    result = [0.0]*m
    for i in range(m):
        for j in range(n):
            result[i] += A[i][j]*x[j]

    return result

def transpose(A: List[List[float]]) -> List[List[float]]:
    """Returns the transpose of a matrix A(n x m): At(m x n)"""

    n, m = len(A), len(A[0])

    At = [ [0.0 for _ in range(n)] for _ in range(m) ]
    for i in range(n):
        for j in range(m):
            At[j][i] = A[i][j]

    return At

def trace(A: List[List[float]]) -> float:
    """Returns the trace of a matrix A(n x n): tr(A)"""

    trA = 0.0
    for i in range(len(A)):
        trA += A[i][i]

    return trA

def minor_matrix(A: List[List[float]], i_ex: int, j_ex: int) -> List[List[float]]:
    """Returns the minors of matrix A by excluding row i_ex and columns j_ex"""

    minor = A[:i_ex] + A[i_ex+1:]
    for i in range(len(minor)):
        minor[i] = minor[i][:j_ex] + minor[i][j_ex+1:]

    return minor

def determinant2D(A: List[List[float]]) -> float:
    """Returns the determinant of a matrix A(2 x 2): det(A)"""

    n, m = len(A), len(A[0])
    if n != 2 or m != 2:
        raise ValueError('Dimensions must be 2x2 for determinant 2D')

    return A[0][0]*A[1][1] - A[0][1]*A[1][0]

def determinant3D(A: List[List[float]]) -> float:
    """Returns the determinant of a matrix A(3 x 3): det(A)"""

    n, m = len(A), len(A[0])
    if n != 3 or m != 3:
        raise ValueError('Dimensions must be 3x3 for determinant 3D')

    a, b, c = A[0][0], A[0][1], A[0][2]

    minor_a = minor_matrix(A, 0, 0)
    minor_b = minor_matrix(A, 0, 1)
    minor_c = minor_matrix(A, 0, 2)

    detA = ( a*determinant2D(minor_a) - b*determinant2D(minor_b)
            + c*determinant2D(minor_c) )

    return detA

def inverse2D(A: List[List[float]]) -> List[List[float]]:
    """Returns the inverse of a matrix A(2 x 2): A^-1"""

    detA = determinant2D(A)
    if detA == 0:
        raise ValueError("2x2 matrix is singular! Cannot be inverted.")

    n, m = len(A), len(A[0])
    if n != 2 or m != 2:
        raise ValueError('Dimensions must be 2x2 for inverse 2D')

    invA = [
        [ A[1][1], -A[0][1] ],
        [ -A[1][0], A[0][0]]
    ]

    return scalar_multiplication(invA, 1.0/detA)

def inverse3D(A: List[List[float]]) -> List[List[float]]:
    """Returns the inverse of a matrix A(3 x 3): A^-1"""

    n, m = len(A), len(A[0])
    if n != 3 or m != 3:
        raise ValueError('Dimensions must be 3x3 for inverse 3D')

    detA = determinant3D(A)
    if detA == 0:
        raise ValueError("3x3 matrix is singular! Cannot be inverted.")

    cofactors = []
    for i in range(n):
        row = []
        for j in range(m):
            minor = minor_matrix(A, i, j)
            sign = (-1)**(i + j)
            row.append(sign*determinant2D(minor))
        cofactors.append(row)

    adjA = transpose(cofactors)

    return scalar_multiplication(adjA, 1.0/detA)

def inverse(A: List[List[float]]) -> List[List[float]]:
    """Returns the inverse of a matrix A using LU decomposition."""

    n = len(A)

    A_lu, rows_order = lu_decomposition(A)

    Ai = [ [0.0 for _ in range(n)] for _ in range(n) ]
    for j in range(n):

        bi = [0.0]*n
        bi[j] = 1.0

        bi_lu = [ bi[rows_order[k]] for k in range(n) ]

        d = direct_solvers.forward_substitution(A_lu, bi_lu)

        x = direct_solvers.back_substitution(A_lu, d)

        for i in range(n):
            Ai[i][j] = x[i]

    return Ai

def euclidean_norm(A: List[List[float]]) -> float:
    """Returns the Euclidean (Frobenius) norm of a matrix A."""

    n, m = len(A), len(A[0])
    normA = 0.0
    for i in range(n):
        for j in range(m):
            normA += A[i][j]**2
    return math.sqrt(normA)

def row_sum_norm(A: List[List[float]]) -> float:
    """Returns the row-sum (infinity) norm of a matrix A."""

    n, m = len(A), len(A[0])
    normA = 0.0
    for i in range(n):
        s_row = 0.0
        for j in range(m):
            s_row += abs(A[i][j])
        normA = max(normA, s_row)
    return normA

def column_sum_norm(A: List[List[float]]) -> float:
    """Returns the column-sum 1-norm of a matrix A."""

    n, m = len(A), len(A[0])
    normA = 0.0
    for j in range(m):
        s_col = 0.0
        for i in range(n):
            s_col += abs(A[i][j])
        normA = max(normA, s_col)
    return normA

def condition_number(A: List[List[float]]) -> float:
    """Returns the condition number of a matrix A."""

    normA = row_sum_norm(A)
    Ai = inverse(A)
    normAi = row_sum_norm(Ai)

    return normA*normAi

def replace_column(A: npt.NDArray[np.float64],
    b: npt.NDArray[np.float64], j: int) -> npt.NDArray[np.float64]:
    """Returns a matrix that has a column, j, replaced by a vector b"""

    A_new = A.copy()
    A_new[:,j] = b
    return A_new

def determinant_upper_diagonal(A: npt.NDArray[np.float64],
    n_swaps: int = 0) -> float:
    """Returns the determinant of an upper diagonal matrix"""

    d = np.diag(A)
    detA = ((-1)**n_swaps)*np.prod(d)
    return detA

def log_determinant_upper_diagonal(A: npt.NDArray[np.float64],
    n_swaps: int = 0) -> Tuple[float, float]:
    """Returns the log-determinant of an upper diagonal matrix"""

    d = np.diag(A)
    log_detA = np.sum(np.log(np.abs(d)))
    n_neg = np.count_nonzero(d < 0)
    sign = (-1)**(n_swaps + n_neg)
    # n = A.shape[0]
    # sign, log_detA = np.linalg.slogdet(A[:,:n])
    return float(sign), float(log_detA)

def partial_pivot(A: npt.NDArray[np.float64], k: int) -> int:
    """Returns the matrix after partial pivoting where
    row k is swapped with the row that has the largest element in column k"""

    i_max = np.argmax(np.abs(A[k:,k])) + k

    if np.abs(A[i_max,k]) < 1e-12:
        raise ValueError('Matrix is nearly singular!')

    return i_max

def lu_decomposition(A: npt.NDArray[np.float64]
    ) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Returns a new matrix that is the LU decomposition of a matrix, A.
    L is the lower and U is the upper diagonal part of the new matrix.
    Also returns the new order of rows from partial pivoting for
    reordering the right-hand side vector later."""

    n = A.shape[0]
    A_lu = A.copy()
    rows_order = np.arange(n)
    for k in range(n):

        i_max = partial_pivot(A_lu, k)
        if i_max != k:
            A_lu[[k, i_max]] = A_lu[[i_max, k]]
            rows_order[[k, i_max]] = rows_order[[i_max, k]]

        for i in range(k+1,n):
            f = A_lu[i][k]/A_lu[k][k]
            A_lu[i,k] = f
            A_lu[i,k+1:] -= f*A_lu[k,k+1:]

    return A_lu, rows_order

def cholesky_decomposition(A: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Returns the lower diagonal matrix, L, from the Cholesky
    decomposition of a symmetric, positive-definite matrix, A."""

    # L = scipy.linalg.cholesky(A)

    n = A.shape[0]
    L = np.zeros((n,n))
    for k in range(n):
        for i in range(k):
            s = np.dot(L[i,:i], L[k,:i])
            L[k,i] = ( A[k,i] - s )/L[i,i]
        s_d = np.dot(L[k,:k], L[k,:k])
        diff = A[k,k] - s_d
        if diff <= 0:
            raise ValueError("Matrix is not positive-definite!")
        L[k,k] = np.sqrt(diff)

    return L

def tri_diagonal(A: npt.NDArray[np.float64]) -> Tuple[npt.NDArray[np.float64]]:
    """Returns the tri-diagonal part of a matrix.
    l, d, u are the lower, diagonal, and upper parts, respectively."""

    l = np.diag(A, k=-1)
    d = np.diag(A, k=0)
    u = np.diag(A, k=1)

    return l, d, u

def qr_decomposition(A: List[List[float]]
    ) -> Tuple[List[List[float]], List[List[float]]]:
    """Returns the QR decomposition of a matrix A using Gram-Schmidt process.
    Q is the orthogonal rotation matrix and R is the upper-triangular."""

    n = len(A)
    Q = [ [0.0]*n for _ in range(n) ]
    R = [ [0.0]*n for _ in range(n) ]

    for j in range(n):
        v = [A[i][j] for i in range(n)]
        for i in range(j):
            R[i][j] = sum(Q[k][i]*A[k][j] for k in range(n))
            for k in range(n):
                v[k] -= R[i][j]*Q[k][i]

        R[j][j] = math.sqrt(sum(x**2 for x in v))
        for i in range(n):
            Q[i][j] = v[i]/R[j][j]
    return Q, R
