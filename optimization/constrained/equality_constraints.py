
from typing import Callable, List, Tuple

from differentiation import forward_fd as ffd
from utilities import vector_operations
from linear_systems import direct_solvers
from non_linear_systems import iterative_solvers

def lagrange_multipliers(f: Callable[[List[float]], float],
    g: Callable[[List[float]], List[float]], x0: List[float] = None, l0: List[float] = None,
    mode: str = 'min', eps: float = 1e-8, k_max = 1000) -> Tuple[List[float], float]:
    """Returns the optimum point, x, of a multi-variable function, f,
    under equality constraints, g, using the Lagrange multipliers method.
    A starting guess point, x0, and Lagrange multipliers, l0, should be given."""

    x = [xi for xi in x0] if x0 else [0.0]*len(x0)
    l = [li for li in l0] if l0 else [0.0]*len(l0)

    nx = len(x)
    nl = len(l)

    for k in range(1, k_max+1):

        res = residual(f,g,x,l)
        res_norm = vector_operations.norm2(res)

        jac = jacobian(residual, f, g, x, l, res)

        b = [-r for r in res]

        dx = direct_solvers.lu_decomposition_solve(jac,b)

        cor_norm = vector_operations.norm2(dx)

        print(f'k = {k}, Res Norm: {res_norm:.4e}, Cor Norm: {cor_norm:.4e}')
        if cor_norm < eps and res_norm < eps:
            break

        for i in range(nx):
            x[i] += dx[i]
        for i in range(nl):
            l[i] += dx[nx+i]

    return x, l, f(x)

def residual(f: Callable[[List[float]], float],
    g: Callable[[List[float]], List[float]], x: List[float] = None,
    l: List[float] = None) -> List[float]:
    """Returns the residual of the unconstrained optimization problem."""

    nx = len(x)
    nl = len(l)

    dfdx = ffd.grad_f_fd(f,x)
    dLdx = [dfdx[i] for i in range(nx)]

    dgdx = iterative_solvers.jacobian(g, x, g(x), h=1e-4)

    for j in range(nl):
        for i in range(nx):
            dLdx[i] += l[j]*dgdx[j][i]

    return dLdx + g(x)

def jacobian(F: Callable[[List[float]], List[float]],
    f: Callable[[List[float]], List[float]], g: Callable[[List[float]], List[float]],
    x: List[float], l: List[float], res: List[float], h: float = 1e-4) -> List[List[float]]:
    """Returns the Jacobian of a system of non-linear algebraic equations
    around the values, x, of the unknowns."""

    nx, nl, m = len(x), len(l), len(res)

    jac = [ [0.0 for _ in range(nx+nl)] for _ in range(m) ]
    for j in range(nx):

        x_f = x[:]
        x_f[j] += h

        res_f = F(f, g, x_f, l)

        for i in range(m):
            jac[i][j] = (res_f[i]-res[i])/h

    for j in range(nl):

        l_f = l[:]
        l_f[j] += h

        res_f = F(f, g, x, l_f)

        for i in range(m):
            jac[i][nx+j] = (res_f[i]-res[i])/h

    return jac
