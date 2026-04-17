
from typing import Tuple
import numpy as np
import numpy.typing as npt
import scipy

def minor_matrix_numpy(A: npt.NDArray[np.float64],
    i: int, j: int) -> npt.NDArray[np.float64]:
    """Returns the minors of matrix A by excluding row i_ex and columns j_ex"""
    A_mod = np.delete(A, i, axis=0)
    return np.delete(A_mod, j, axis=1)

def inverse(A: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Returns the inverse of a matrix A using LU decomposition."""

    n = len(A)

    A_lu, rows_order = lu_decomposition(A)

    Ai = np.zeros((n,n))
    for j in range(n):

        bi = [0.0]*n
        bi[j] = 1.0

        bi_lu = [ bi[rows_order[k]] for k in range(n) ]

        d = forward_substitution(A_lu, bi_lu)

        Ai[:,j] = back_substitution(A_lu, d)

    return Ai

def euclidean_norm(A: npt.NDArray[np.float64]) -> float:
    """Returns the Euclidean (Frobenius) norm of a matrix A."""

    return np.linalg.norm(A, ord='fro')
    # return np.sqrt(np.sum(A**2))

def row_sum_norm(A: npt.NDArray[np.float64]) -> float:
    """Returns the row-sum (infinity) norm of a matrix A."""

    # normA = np.linalg.norm(A, ord=np.inf)
    # normA = np.sum(np.abs(A), axis=1).max()

    n = A.shape[0]
    normA = 0
    for i in range(n):
        s_row = np.sum(np.abs(A[i,:]))
        normA = np.maximum(normA, s_row)
    return normA

def column_sum_norm(A: npt.NDArray[np.float64]) -> float:
    """Returns the column-sum 1-norm of a matrix A."""

    # normA = np.sum(np.abs(A), axis=0).max()

    m = A.shape[1]
    normA = 0
    for j in range(m):
        s_col = np.sum(np.abs(A[:,j]))
        normA = np.maximum(normA, s_col)
    return normA

def condition_number(A: npt.NDArray[np.float64]) -> float:
    """Returns the condition number of a matrix A."""

    # cond = np.linalg.cond(A, p=np.inf) # norm_oo: row-sum norm

    # cond = np.linalg.cond(A, p=2) # norm-2 from SVD
    # s = np.linalg.svd(A, compute_uv=False)
    # # Condition number is the ratio of max to min
    # cond = s[0]/s[-1]

    normA = row_sum_norm(A)
    normAi = row_sum_norm(inverse(A))
    cond = normA*normAi

    return cond

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

def forward_substitution(A: npt.NDArray[np.float64],
    b: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Returns the intermediate vector, d, from Ld = b.
    L is the lower diagonal part of matrix A and b is the right-hand side vector."""

    n = A.shape[0]
    d = np.zeros(n)
    for i in range(n):
        s = np.dot(A[i,:i],d[:i])
        d[i] = b[i]-s
    return d

def back_substitution(A: npt.NDArray[np.float64],
    b: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Returns the solution vector, x, of a linear system Ux = b.
    U is the upper diagonal part of matrix A and b is the right-hand side vector."""

    n = A.shape[0]
    x = np.zeros(n)
    for i in range(n-1,-1,-1):
        s = np.dot(A[i, i+1:], x[i+1:])
        x[i] = (b[i]-s)/A[i,i]

    return x

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
            f = A_lu[i,k]/A_lu[k,k]
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

def qr_decomposition(A: npt.NDArray[np.float64]
    ) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Returns the QR decomposition of a matrix A using Gram-Schmidt process.
    Q is the orthogonal rotation matrix and R is the upper-triangular."""

    n = A.shape[1]
    Q = np.zeros_like(A)
    R = np.zeros((n, n))

    for j in range(n):
        v = A[:,j]
        for i in range(j):
            R[i][j] = np.dot(Q[:,i], A[:,j])
            v = v - R[i,j]*Q[:,i]

        R[j][j] = np.linalg.norm(v)
        Q[:,j] = v/R[j,j]
    return Q, R
