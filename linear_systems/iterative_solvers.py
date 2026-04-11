
from typing import List

from utilities import vector_operations

def jacobi(A: List[List[float]], b: List[float],
    x0: List[float] = None, output: bool = False,
    eps: float = 1e-8, k_max: int = 1000) -> List[float]:
    """Returns the solution, x, to a system of linear algebraic equations,
    Ax = b, using the Jacobi iterative method."""

    n = len(A)
    if x0 is not None:
        x = [xi for xi in x0]
    else:
        x = [0.0]*n
        for i in range(n):
            x[i] = b[i]/(A[i][i]+eps)

    res = [0.0]*n

    for k in range(1, k_max+1):

        x_old = x[:]
        for i in range(n):

            res[i] = b[i] - vector_operations.dot_product(A[i], x_old)

            x[i] = x_old[i] + res[i]/(A[i][i]+eps)

        # for i in range(n):
        #     res[i] = b[i] - vector_operations.dot_product(A[i], x)
        res_norm = vector_operations.norm2(res)

        dx = [ x[i]-x_old[i] for i in range(n) ]
        cor_norm = vector_operations.norm2(dx)

        converged = (res_norm < eps) and (cor_norm < eps)
        if converged:
            if output:
                print(f'Jacobi: k = {k}, Res Norm = {res_norm:.4e}, Cor Norm = {cor_norm:.4e}')
            break

    return x

def sor_solver(A: List[List[float]], b: List[float], x0: List[float] = None,
    w: float = 1.0, output: bool = False,
    eps: float = 1e-8, k_max: int = 1000) -> List[float]:
    """Returns the solution, x, to a system of linear algebraic equations,
    Ax = b, using the SOR iterative method. The relaxation can be optionally
    given. For ω = 1, the Gauss-Seidel method is used."""

    n = len(A)
    if x0 is not None:
        x = [xi for xi in x0]
    else:
        x = [1.0]*n
        for i in range(n):
            x[i] = b[i]/(A[i][i]+eps)

    res = [0.0]*n

    for k in range(1, k_max+1):

        x_old = x[:]
        for i in range(n):

            # res[i] = b[i] - sum(A[i][j]*x[j] for j in range(n))
            res[i] = b[i] - vector_operations.dot_product(A[i], x)

            x[i] += w*res[i]/(A[i][i]+eps)

        res_norm = vector_operations.norm2(res)

        dx = [ x[i]-x_old[i] for i in range(n) ]
        cor_norm = vector_operations.norm2(dx)

        converged = (res_norm < eps) and (cor_norm < eps)
        if converged:
            if output:
                print(f'SOR: k = {k}, Res Norm = {res_norm:.4e}, Cor Norm = {cor_norm:.4e}')
            break

    return x

def steepest_descent(A: List[List[float]], b: List[float], x0: List[float] = None,
    output: bool = False, eps: float = 1e-8, k_max: int = 1000) -> List[float]:
    """Returns the solution, x, to a system of linear algebraic equations,
    Ax = b, using the Steepest Descent iterative method."""

    n = len(A)

    if x0 is not None:
        x = [xi for xi in x0]
    else:
        x = [0.0]*n
        for i in range(n):
            x[i] = b[i]/(A[i][i]+eps)

    r = [bi for bi in b]
    for i in range(n):
        r[i] -= vector_operations.dot_product(A[i],x)

    Ar = [0.0]*n
    dx = [0.0]*n
    for k in range(1, k_max+1):

        rr = vector_operations.dot_product(r, r)

        for i in range(n):
            Ar[i] = vector_operations.dot_product(A[i], r)

        rAr = vector_operations.dot_product(r, Ar)

        alpha = rr/(rAr+eps)

        x_old = x[:]

        for i in range(n):
            x[i] += alpha*r[i]
            r[i] -= alpha*Ar[i]
            dx[i] = x[i] - x_old[i]

        res_norm = rr**0.5
        cor_norm = vector_operations.norm2(dx)

        converged = (res_norm < eps) or (cor_norm < eps)
        if converged:
            if output:
                print(f'Steepest Descent: k = {k}, Res Norm = {res_norm:.4e}, Cor Norm = {cor_norm:.4e}')
            break

    return x

def conjugate_gradient(A: List[List[float]], b: List[float], x0: List[float] = None,
    output: bool = False, eps: float = 1e-8, k_max: int = 1000) -> List[float]:
    """Returns the solution, x, to a system of linear algebraic equations,
    Ax = b, using the Conjugate Gradient iterative method."""

    n = len(A)

    if x0 is not None:
        x = [xi for xi in x0]
    else:
        x = [0.0]*n
        for i in range(n):
            x[i] = b[i]/(A[i][i]+eps)

    r = b[:]
    for i in range(n):
        r[i] -= vector_operations.dot_product(A[i],x)

    rr_old = vector_operations.dot_product(r, r)

    p = r[:]

    Ap = [0.0]*n
    dx = [0.0]*n
    for k in range(1, k_max+1):

        for i in range(n):
            Ap[i] = vector_operations.dot_product(A[i], p)

        pAp = vector_operations.dot_product(p, Ap)

        alpha = rr_old/(pAp+eps)

        x_old = x[:]

        for i in range(n):
            x[i] += alpha*p[i]
            r[i] -= alpha*Ap[i]
            dx[i] = x[i] - x_old[i]

        rr = vector_operations.dot_product(r,r)

        res_norm = rr**0.5
        cor_norm = vector_operations.norm2(dx)

        converged = (res_norm < eps) or (cor_norm < eps)
        if converged:
            if output:
                print(f'Conjugate Gradient: k = {k}, Res Norm = {res_norm:.4e}, Cor Norm = {cor_norm:.4e}')
            break

        beta = rr/rr_old
        for i in range(n):
            p[i] = r[i] + beta*p[i]

        rr_old = rr

    return x
