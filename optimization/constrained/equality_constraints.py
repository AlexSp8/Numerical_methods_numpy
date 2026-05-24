
from typing import Callable
import numpy as np
import numpy.typing as npt

from differentiation import partial_derivatives
from linear_systems import direct_solver
from non_linear_systems import newton_solver, non_linear_problem
from optimization.unconstrained import gradient_methods

def lagrange_multipliers(f: Callable[[npt.NDArray[np.float64]], float],
    g: Callable[[npt.NDArray[np.float64]], npt.NDArray[np.float64]],
    x0: npt.NDArray[np.float64],
    df: Callable[[Callable[[float], float], float, float, float], float],
    l0_m: npt.NDArray[np.float64] = None, output: bool = False
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Returns the optimum point and corresponding Lagrange multipliers
    of a multi-variable function under equality constraints using the
    Lagrange multipliers method.

    Args:
        f: the function to be optimized
        g: the constraints function
        x0: initial guess
        fd: the function to calculate the derivative of f
        l0_m (optional): initial guess for the Lagrange multipliers
        output (optional): print iteration info
    Returns:
        The optimum point and corresponding Lagrange multipliers."""

    nx = x0.shape[0]

    m = len(g(x0))

    if l0_m is None:
        l_m = np.zeros(m)
    else:
        l_m = l0_m.copy()

    # u0 = np.concatenate([x0, l_m])
    u0 = np.zeros(nx+m)
    u0[:nx] = x0.copy()
    u0[nx:] = l_m.copy()

    ls_solver = direct_solver.LUSolver()

    nr_solver = newton_solver.NewtonSolver(ls_solver, u0=u0, k_max=100, tol=1e-8, r=1.0)

    problem = non_linear_problem.LagrangeMultipliers(
        f=f, nx=nx, g=g, df=df, grad_f=partial_derivatives.grad_f_fd)

    u = nr_solver.solve(problem, output=output)

    x, l_m = u[:nx], u[nx:]

    return x, l_m

def augmented_lagrangian(f: Callable[[npt.NDArray[np.float64]], float],
    g_con: Callable[[npt.NDArray[np.float64]], npt.NDArray[np.float64]],
    x0: npt.NDArray[np.float64],
    df: Callable[[Callable[[float], float], float, float, float], float],
    l0_m: npt.NDArray[np.float64] = None, output: bool = False, tol: float = 1e-6,
    k_max: int = 100) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Returns the optimum point and corresponding Lagrange multipliers
    of a multi-variable function under equality constraints using the
    Augmented Lagrangian method.

    Args:
        f: the function to be optimized
        g: the constraints function
        x0: initial guess
        fd: the function to calculate the derivative of f
        l0_m (optional): initial guess for the Lagrange multipliers
        output (optional): print iteration info
        tol (optional): threshold for convergence
        k_max (optional): maximum number of iterations
    Returns:
        The optimum point and corresponding Lagrange multipliers."""

    x = x0.copy()

    m = len(g_con(x0))

    if l0_m is None:
        l = np.zeros(m)
    else:
        l = l0_m.copy()

    p = 20.0

    for k in range(1, k_max+1):

        # def f_lagrangian(x_c: npt.NDArray[np.float64]) -> float:
        #     v = g_con(x_c)
        #     return f(x_c) + np.dot(l,v) + 0.5*p*np.sum(v**2)

        aug_l = AugmentedLagrangian(f, g_con, l, p)

        x_new = gradient_methods.bfgs(f=aug_l.f_lagrangian, x0=x, df=df)

        g = g_con(x_new)
        norm_g = np.linalg.norm(g)

        l_old = l
        l = l + p*g

        dl_norm = np.linalg.norm(l-l_old)
        dx_norm = np.linalg.norm(x_new-x)

        if output:
            print(f'k = {k}, Norm g: {norm_g:.4e}, dl Norm: {dl_norm:.4e}, dx Norm: {dx_norm:.4e}')

        if norm_g < tol and dl_norm < tol and dx_norm < tol:
            print(f'k = {k}, Norm g: {norm_g:.4e}, dl Norm: {dl_norm:.4e}, dx Norm: {dx_norm:.4e}')
            return x_new, l

        if p < 1e3:
            mult = 1.2
        else:
            mult = 1.02
        p *= mult
        p = min(p, 1e8)

        x = x_new

    print(f"Warning: Maximum iterations ({k_max}) reached without converging.")

    return x, l


class AugmentedLagrangian:

    """An object wrapper that behaves like a function but holds state parameters."""

    def __init__(self, f: Callable[[npt.NDArray[np.float64]], float],
        g_con: Callable[[npt.NDArray[np.float64]], npt.NDArray[np.float64]],
        l: npt.NDArray[np.float64], p: float):

        self.f = f
        self.g_con = g_con
        self.l = l
        self.p = p

    # def __call__(self, x: npt.NDArray[np.float64]) -> float:

    #     g = self.g_con(x)
    #     return self.f(x) + np.dot(self.l, g) + 0.5*self.p*np.sum(g**2)

    def f_lagrangian(self, x: npt.NDArray[np.float64]) -> float:

        g = self.g_con(x)
        return self.f(x) + np.dot(self.l, g) + 0.5*self.p*np.sum(g**2)
