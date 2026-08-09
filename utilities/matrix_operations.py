
import numpy as np
from numpy.typing import NDArray
import scipy

def minor_matrix_numpy(A: NDArray[np.float64],
    i: int, j: int) -> NDArray[np.float64]:
    """Returns the minors of matrix A by excluding row i_ex and columns j_ex"""
    A_mod = np.delete(A, i, axis=0)
    return np.delete(A_mod, j, axis=1)

def inverse(A: NDArray[np.float64]
    ) -> NDArray[np.float64]:
    """Returns the inverse of a matrix A using LU decomposition.

    Args:
        A: the matrix
    Returns:
        The inverse of the matrix."""

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

def euclidean_norm(A: NDArray[np.float64]) -> float:
    """Returns the Euclidean (Frobenius) norm of a matrix A.

    Args:
        A: the matrix
    Returns:
        The Euclidean norm of the matrix."""

    # return np.linalg.norm(A, ord='fro')
    return np.sqrt(np.sum(A**2))

def row_sum_norm(A: NDArray[np.float64]) -> float:
    """Returns the row-sum (infinity) norm of a matrix A.

    Args:
        A: the matrix
    Returns:
        The row-sum (infinity) norm of the matrix."""

    # normA = np.linalg.norm(A, ord=np.inf)
    # normA = np.sum(np.abs(A), axis=1).max()

    n = A.shape[0]
    normA = 0
    for i in range(n):
        s_row = np.sum(np.abs(A[i,:]))
        normA = np.maximum(normA, s_row)
    return normA

def column_sum_norm(A: NDArray[np.float64]) -> float:
    """Returns the column-sum 1-norm of a matrix A.

    Args:
        A: the matrix
    Returns:
        The column-sum 1-norm of the matrix."""

    # normA = np.sum(np.abs(A), axis=0).max()

    m = A.shape[1]
    normA = 0
    for j in range(m):
        s_col = np.sum(np.abs(A[:,j]))
        normA = np.maximum(normA, s_col)
    return normA

def condition_number(A: NDArray[np.float64]) -> float:
    """Returns the condition number of a matrix A.

    Args:
        A: the matrix
    Returns:
        The condition number of the matrix."""

    # cond = np.linalg.cond(A, p=np.inf) # norm_oo: row-sum norm

    # cond = np.linalg.cond(A, p=2) # norm-2 from SVD
    # s = np.linalg.svd(A, compute_uv=False)
    # # Condition number is the ratio of max to min
    # cond = s[0]/s[-1]

    normA = row_sum_norm(A)
    normAi = row_sum_norm(inverse(A))
    cond = normA*normAi

    return cond

def replace_column(A: NDArray[np.float64],
    b: NDArray[np.float64], j: int
    ) -> NDArray[np.float64]:
    """Returns a matrix that has a column, j, replaced by a vector b.

    Args:
        A: the matrix
        b: the vector to be used for the replacement
        j: the column of the matrix to be replaced

    Returns:
        The new matrix containing vector b in column j."""

    A_new = A.copy()

    A_new[:,j] = b.copy()

    return A_new

def determinant_upper_diagonal(A: NDArray[np.float64],
    n_swaps: int = 0) -> float:
    """Returns the determinant of an upper diagonal matrix.

    Args:
        A: the matrix
        n_swaps: the number of row swaps during the partial pivot phase

    Returns:
        The determinant of matrix A."""

    d = np.diag(A)

    p = np.prod(d)

    # p = 1
    # n = A.shape[0]
    # for i in range(n):
    #     p *= A[i,i]

    detA = ((-1)**n_swaps)*p

    return detA

def log_determinant_upper_diagonal(A: NDArray[np.float64],
    n_swaps: int = 0) -> tuple[float, float]:
    """Returns the log-determinant of an upper diagonal matrix.

    Args:
        A: the matrix
        n_swaps: the number of row swaps during the partial pivot phase

    Returns:
        The log-determinant of matrix A."""

    # n = A.shape[0]
    # sign, log_detA = np.linalg.slogdet(A[:,:n])

    d = np.diag(A)

    log_detA = np.sum(np.log(np.abs(d)))

    n_neg = np.count_nonzero(d < 0)

    # n = d.shape[0]
    # n_neg = 0
    # for i in range(n):
    #     if d[i] < 0:
    #         n_neg += 1

    sign = (-1)**(n_swaps + n_neg)

    return float(sign), float(log_detA)

def pivot_row(A: NDArray[np.float64],
    k: int) -> int:
    """Returns the row of a matrix below the diagonal
    that has the largest element in column k.

    Args:
        A: the matrix
        k: the column

    Returns:
        The row of the matrix below the diagonal
        which has the largest element in column k."""

    i_max = np.argmax(np.abs(A[k:,k])) + k

    # i_max = k
    # max_val = abs(A[k,k])
    # n = A.shape[0]
    # for i in range(k+1, n):
    #     val = abs(A[i,k])
    #     if val > max_val:
    #         max_val = val
    #         i_max = i

    if np.abs(A[i_max,k]) < 1e-12:
        raise ValueError('Matrix is nearly singular!')

    return i_max

def forward_substitution(A: NDArray[np.float64],
    b: NDArray[np.float64]
    ) -> NDArray[np.float64]:
    """Returns the intermediate vector, d, from Ld = b after forward substitution.
    L is the lower diagonal part of matrix A and b is the right-hand side vector.

    Args:
        A: the matrix of the linear system
        b: the right-hand side vector of the linear system
    Returns:
        The intermediate vector, d, after forward substitution."""

    n = A.shape[0]

    d = np.zeros(n)

    for i in range(n):
        s = np.dot(A[i,:i], d[:i])
        d[i] = b[i]-s

    return d

def back_substitution(A: NDArray[np.float64],
    b: NDArray[np.float64]
    ) -> NDArray[np.float64]:
    """Returns the vector x of a linear system Ax = b after back substitution.
    This is the solution vector if forward substitution has been performed in
    A to convert it into upper diagonal.

    Args:
        A: the matrix of the linear system
        b: the right-hand side vector of the linear system
    Returns:
        The solution vector, x, after back substitution."""

    n = A.shape[0]

    x = np.zeros(n)

    for i in range(n-1,-1,-1):
        s = np.dot(A[i, i+1:], x[i+1:])
        x[i] = (b[i]-s)/A[i,i]

    return x

def lu_decomposition(A: NDArray[np.float64],
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Returns the LU decomposition of a matrix, A and the order of rows after partial pivot.
    L is the lower and U is the upper diagonal part of the new matrix.
    The new order of rows can be used for reordering the right-hand side vector later.

    Args:
        A: the matrix
    Returns:
        The LU decomposition of matrix A and the new order of rows of matrix A after partial pivot."""

    n = A.shape[0]

    A_lu = A.copy()

    rows_order = np.arange(n)
    # rows_order = np.zeros(n, dtype=np.intp)
    # for i in range(n):
    #     rows_order[i] = i

    for k in range(n):

        i_max = pivot_row(A_lu, k)
        if i_max != k:
            A_lu[[k, i_max]] = A_lu[[i_max, k]]
            rows_order[[k, i_max]] = rows_order[[i_max, k]]

        for i in range(k+1, n):
            f = A_lu[i,k]/A_lu[k,k]
            A_lu[i,k] = f
            A_lu[i,k+1:] -= f*A_lu[k,k+1:]

    return A_lu, rows_order

def tri_diagonal(A: NDArray[np.float64]) -> tuple[NDArray[np.float64]]:
    """Returns the tri-diagonal part of a matrix.

    Args:
        A: the matrix
    Returns:
        l: the sub-diagonal elements of the matrix
        d: the diagonal elements of the matrix
        u: the upper-diagonal elements of the matrix."""

    # l = np.diag(A, k=-1)
    # d = np.diag(A, k=0)
    # u = np.diag(A, k=1)

    n = A.shape[0]
    l = np.zeros(n-1)
    d = np.zeros(n)
    u = np.zeros(n-1)
    for i in range(n):
        d[i] = A[i,i]
        if i > 0:
            l[i-1] = A[i,i-1]
        if i < n-1:
            u[i] = A[i,i+1]

    return l, d, u

def cholesky_decomposition(A: NDArray[np.float64]
    ) -> NDArray[np.float64]:
    """Returns the lower diagonal matrix, L, from the Cholesky
    decomposition of a symmetric, positive-definite matrix, A.

    Args:
        A: the matrix
    Returns:
        The lower diagonal matrix from the Cholesky decomposition."""

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

def qr_decomposition(A: NDArray[np.float64]
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Returns the QR decomposition of a matrix A using Gram-Schmidt process.
    Q is the orthogonal rotation matrix and R is the upper-triangular.

    Args:
        A: the matrix on which QR decomposition is performed
    Returns:
        A tuple of the Q and R matrices."""

    n = A.shape[1]
    Q = np.zeros_like(A)
    R = np.zeros((n, n))

    for j in range(n):
        v = A[:,j]
        for i in range(j):
            R[i,j] = np.dot(Q[:,i], A[:,j])
            v = v - R[i,j]*Q[:,i]

        R[j,j] = np.linalg.norm(v)
        Q[:,j] = v/R[j,j]

    return Q, R
