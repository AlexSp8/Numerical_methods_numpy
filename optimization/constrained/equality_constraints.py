
from typing import Callable, Tuple
import numpy as np
import numpy.typing as npt

from differentiation import forward_fd as ffd, partial_derivatives
from differentiation import backward_fd as bfd
from differentiation import central_fd as cfd
from linear_systems import direct_solver
from non_linear_systems import newton_solver, non_linear_problem
from optimization.unconstrained import gradient_methods

def lagrange_multipliers(f: Callable[[npt.NDArray[np.float64]], float],
    g: Callable[[npt.NDArray[np.float64]], npt.NDArray[np.float64]],
    x0: npt.NDArray[np.float64], l0: npt.NDArray[np.float64] = None,
    fd_type: str = 'ffd', output: bool = False
    ) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Returns the optimum point, x, of a multi-variable function, f,
    under equality constraints, g, using the Lagrange multipliers method.
    Starting point, x0, and (optionally) Lagrange multipliers, l0, should be given."""

    ls_solver = direct_solver.LUSolver()

    m = len(g(x0))

    if l0 is None:
        l = np.zeros(m)
    else:
        l = l0

    u0 = np.concatenate([x0, l])

    nr_solver = newton_solver.NewtonSolver(ls_solver=ls_solver, u0=u0, k_max=100, tol=1e-8, r=1.0)

    if fd_type == 'ffd':
        df = ffd.df_h
    elif fd_type == 'bfd':
        df = bfd.df_h
    else:
        df = cfd.df_h

    nu = x0.shape[0]
    problem = non_linear_problem.LagrangeMultiplier(
        f=f, nu=nu, g=g, df=df, grad_f=partial_derivatives.grad_f_fd)

    u = nr_solver.solve(problem, output=output)

    x = u[:nu]
    l = u[nu:]

    return x, l

def augmented_lagrangian(f: Callable[[npt.NDArray[np.float64]], float],
    g_con: Callable[[npt.NDArray[np.float64]], npt.NDArray[np.float64]],
    x0: npt.NDArray[np.float64], mode: str = 'min', output: bool = False,
    l0: npt.NDArray[np.float64] = None, tol: float = 1e-8, k_max: int = 100):
    """Returns the optimum point, x, of a multi-variable function, f,
    under equality constraints, g, using the Augmented Lagrangian method.
    Starting point, x0, and (optionally) Lagrange multipliers, l0, should be given."""

    x = x0.copy()

    g = g_con(x0)
    m = len(g)

    if l0 is None:
        l = np.zeros(m)
    else:
        l = l0

    p = 20.0

    if mode == 'min':
        s = +1.0
    else:
        s = -1.0

    for k in range(1, k_max+1):

        def current_lagrangian(x_c: npt.NDArray[np.float64]) -> float:
            v = g_con(x_c)
            return f(x_c) + s*np.dot(l,v) + s*p*np.sum(v**2)

        x = gradient_methods.bfgs(f=current_lagrangian, x0=x, fd_type='ffd', mode=mode)

        g = g_con(x)
        norm_g = np.linalg.norm(g)

        l_old = l
        l = l + p*g

        cor_norm = np.linalg.norm(l-l_old)

        if output:
            print(f'k = {k}, Norm g: {norm_g:.4e}, Cor Norm: {cor_norm:.4e}')

        converged = norm_g < tol or cor_norm < tol
        if converged:
            print(f'k = {k}, Norm g: {norm_g:.4e}, Cor Norm: {cor_norm:.4e}')
            break

        if p < 1e3:
            mult = 1.2
        else:
            mult = 1.02
        p *= mult
        p = min(p, 1e8)

    return x, l
