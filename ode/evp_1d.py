
from typing import List, Callable, Tuple

from utilities import matrix_operations

def solve_sturm_liouville(x: List[float], p: Callable[[List[float]], float],
    q: Callable[[List[float]], float], w: Callable[[List[float]], float]
    ) -> Tuple[List[float], List[List[float]]]:
    """Returns the solution (eigenvalues and eigenvectors) of a general
    linear 2nd order ODE using Finite Differences. The domain, x, and the
    functions p (stiffness/diffusion), q (potential/restoration),
    w (weight/mass density) must be given."""

    n = len(x)
    n_int = n-2  # Internal nodes only (Dirichlet BCs: y=0)

    A = [ [0.0 for _ in range(n_int)] for _ in range(n_int) ]
    B = [0.0 for _ in range(n_int)]

    for i in range(n_int):

        idx = i + 1
        hb = x[idx] - x[idx-1]
        hf = x[idx+1] - x[idx]
        dx_avg = (hb + hf)/2.0

        pb = p((x[idx] + x[idx-1])/2.0)
        pf = p((x[idx+1] + x[idx])/2.0)

        if i > 0:
            A[i][i-1] = -pb/(dx_avg*hb)

        if i < n_int - 1:
            A[i][i+1] = -pf/(dx_avg*hf)

        A[i][i] = (1.0/dx_avg)*(pf/hf + pb/hb) - q(x[idx])

        B[i] = w(x[idx])

    eigenvals, eigenvecs = solve_generalized_evp(A, B)

    return eigenvals, eigenvecs

def solve_generalized_evp(A: List[List[float]], B: List[float],
    iterations = 100) -> Tuple[List[float], List[List[float]]]:

    n = len(A)

    B_inv_sqrt = [1.0/(w**0.5) for w in B]

    A_new = [ [0.0]*n for _ in range(n) ]
    for i in range(n):
        for j in range(n):
            A_new[i][j] = A[i][j]*B_inv_sqrt[i]*B_inv_sqrt[j]

    lambdas, Z_vectors = solve_evp(A_new, iterations)

    # Convert Z back to physical eigenvectors Y
    # Y = B_inv_sqrt * Z
    Y_vectors = [[0.0]*n for _ in range(n)]
    for j in range(n):
        for i in range(n):
            Y_vectors[i][j] = Z_vectors[i][j]*B_inv_sqrt[i]

    # Transpose Y_vectors so each element in the list is a full eigenvector
    modes = [[Y_vectors[i][j] for i in range(n)] for j in range(n)]

    # Pair each eigenvalue with its mode, sort by the eigenvalue, then unpack
    combined = sorted(zip(lambdas, modes), key=lambda pair: pair[0])

    sorted_eigenvals, sorted_modes = zip(*combined)

    return list(sorted_eigenvals), list(sorted_modes)

def solve_evp(A: List[List[float]],
    iterations = 100) -> Tuple[List[float], List[List[float]]]:

    n = len(A)
    A_new = [row[:] for row in A]
    V = [ [1.0 if i == j else 0.0 for j in range(n)] for i in range(n) ]

    for _ in range(iterations):
        Q, R = matrix_operations.qr_decomposition(A_new)
        A_new = matrix_operations.multiplication(R, Q)
        V = matrix_operations.multiplication(V, Q)

    # The eigenvalues are the diagonal elements
    eigenvals = [A_new[i][i] for i in range(n)]
    return eigenvals, V
