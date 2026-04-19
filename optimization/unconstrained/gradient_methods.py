
from typing import Callable, Tuple
import numpy as np
import numpy.typing as npt
from scipy.optimize import minimize

from differentiation import partial_derivatives, forward_fd as ffd
from differentiation import backward_fd as bfd, central_fd as cfd
from optimization.one_dimensional import hybrid
from optimization.unconstrained.line_search import LineSearch
from linear_systems import direct_solver
from non_linear_systems import newton_solver, non_linear_problem

def steepest_descent(f: Callable[[npt.NDarray[np.float64]], float],
    x0: npt.NDarray[np.float64] = None, fd_type: str = 'ffd', mode: str = 'min',
    tol: float = 1e-8, k_max = 1000) -> npt.NDarray[np.float64]:
    """Returns the point x(x1, x2, ..., xn) where the extreme of a
    multi-variable function, f, is found using the Steepest Descent method.
    A starting guess point, x0 should be given."""

    if mode == 'min':
        s = -1.0
    else:
        s = +1.0

    if fd_type == 'ffd':
        df = ffd.df_h
    elif fd_type == 'bfd':
        df = bfd.df_h
    else:
        df = cfd.df_h

    x = x0.copy()

    line_search = LineSearch(f)

    for k in range(1, k_max+1):

        grad_f = partial_derivatives.grad_f_fd(df,f,x)

        norm_g = np.linalg.norm(grad_f)
        if norm_g < tol:
            print('k =', k)
            break

        grad_f = s*grad_f/norm_g

        line_search.update(x, grad_f)

        h_min, h_max = h_interval(line_search.f_line, s)

        h_opt = hybrid.brent(line_search.f_line, h_min, h_max, mode)

        x = x + h_opt*grad_f

    return x

def h_interval(f_line: Callable[[float], float],
    s: float = -1.0, k_max = 100) -> Tuple[float, float]:
    """Returns the interval [h_min, h_max] for 1D Line Search."""

    h_min, h_max = 0.0, 0.1
    f_min = f_line(h_min)
    for k in range(k_max):
        f_max = f_line(h_max)
        expand = (s*f_max > s*f_min)
        if expand:
            h_min = h_max
            f_min = f_max
            h_max *= 10.0
        else:
            # print('h_interval iter k =', k)
            break
    return h_min, h_max

def conjugate_gradient(f: Callable[[npt.NDArray[np.float64]], float],
    x0: npt.NDArray[np.float64] = None, fd_type: str = 'ffd', mode: str = 'min',
    b_update: str = 'polak_ribiere', tol: float = 1e-8,
    k_max = 1000) -> npt.NDArray[np.float64]:
    """Returns the point x(x1, x2, ..., xn) where the extreme of a
    multi-variable function, f, is found using the Conjugate Gradient method.
    A starting guess point, x0 should be given."""

    # if mode == 'min':
    #     g = f
    # else:
    #     g = lambda x: -f(x)
    # res = minimize(g, x0, method='CG')
    # return res.x

    if mode == 'min':
        s = -1.0
    else:
        s = +1.0

    if fd_type == 'ffd':
        df = ffd.df_h
    elif fd_type == 'bfd':
        df = bfd.df_h
    else:
        df = cfd.df_h

    x = x0.copy()

    grad_f_old = partial_derivatives.grad_f_fd(df,f,x)

    p = s*grad_f_old

    line_search = LineSearch(f)

    for k in range(1, k_max+1):

        norm_p = np.linalg.norm(p)
        if norm_p < tol:
            print('k =', k)
            break

        # p = p/norm_p

        line_search.update(x, p)

        h_min, h_max = h_interval(line_search.f_line, s)

        h_opt = hybrid.brent(line_search.f_line, h_min, h_max, mode)

        x = x + h_opt*p

        grad_f = partial_derivatives.grad_f_fd(df,f,x)

        if b_update == 'polak_ribiere':
            num = np.dot(grad_f, grad_f-grad_f_old)
            den = np.dot(grad_f_old, grad_f_old) + 1e-12
        elif b_update == 'fletcher_reeves':
            num = np.dot(grad_f, grad_f)
            den = np.dot(grad_f_old, grad_f_old) + 1e-12
        else:
            raise ValueError('Wrong b_update choice in CG!')

        beta = num/den
        p = s*grad_f + beta*p

        grad_f_old = grad_f.copy()

    return x

def newton(f: Callable[[npt.NDArray[np.float64]], float],
    x0: npt.NDArray[np.float64], fd_type: str = 'ffd',
    mode: str = 'min') -> npt.NDArray[np.float64]:
    """Returns the point x(x1, x2, ..., xn) where the extreme of a
    multi-variable function, f, is found using Newton's method.
    A starting guess point, x0 should be given."""

    if fd_type == 'ffd':
        df = ffd.df_h
    elif fd_type == 'bfd':
        df = bfd.df_h
    else:
        df = cfd.df_h

    grad_f = partial_derivatives.grad_f_fd
    hessian_f = partial_derivatives.hessian_f_fd

    # if mode == 'min':
    #     s = -1.0
    # else:
    #     s = +1.0

    # g = lambda x: -s*f(x)
    # jac_g = lambda x: -s*grad_f(df, f, x)
    # hess_g = lambda x: -s*hessian_f(df, f, x)
    # res = minimize(g, x0, method='Newton-CG', jac=jac_g, hess=hess_g)
    # return res.x

    ls_solver = direct_solver.LUSolver()

    nr_solver = newton_solver.NewtonSolver(ls_solver=ls_solver, u0=x0, k_max=100, tol=1e-8, r=1.0)

    problem = non_linear_problem.Optimization(f=f, df=df, grad_f=grad_f, hessian_f=hessian_f)

    x = nr_solver.solve(problem, output=True)

    return x

def marquardt(f: Callable[[npt.NDArray[np.float64]], float],
    x0: npt.NDArray[np.float64], fd_type: str = 'ffd'
    ) -> npt.NDArray[np.float64]:
    """Returns the point x(x1, x2, ..., xn) where the extreme of a
    multi-variable function, f, is found using the Marquardt method.
    A starting guess point, x0 should be given."""

    ls_solver = direct_solver.LUSolver()

    nr_solver = newton_solver.MarquardtSolver(ls_solver=ls_solver, u0=x0,
        k_max=100, tol=1e-8, r=1.0, l0=1e-4, scale=2.0)

    if fd_type == 'ffd':
        df = ffd.df_h
    elif fd_type == 'bfd':
        df = bfd.df_h
    else:
        df = cfd.df_h

    grad_f = partial_derivatives.grad_f_fd
    hessian_f = partial_derivatives.hessian_f_fd

    problem = non_linear_problem.Optimization(f=f, df=df, grad_f=grad_f, hessian_f=hessian_f)

    x = nr_solver.solve(problem, output=True)
    return x

def bfgs(f: Callable[[npt.NDArray[np.float64]], float],
    x0: npt.NDArray[np.float64] = None, fd_type: str = 'ffd', mode: str = 'min',
    output: bool = False, tol: float = 1e-8, k_max = 1000) -> npt.NDArray[np.float64]:
    """Returns the point x(x1, x2, ..., xn) where the extreme of a
    multi-variable function, f, is found using the BFGS method.
    A starting guess point, x0 should be given."""

    # if mode == 'min':
    #     g = f
    # else:
    #     g = lambda x: -f(x)
    # res = minimize(g, x0, method='BFGS')
    # return res.x

    n = len(x0)
    x = x0.copy()

    if mode == 'min':
        s = -1.0
    else:
        s = +1.0

    if fd_type == 'ffd':
        df = ffd.df_h
    elif fd_type == 'bfd':
        df = bfd.df_h
    else:
        df = cfd.df_h

    H = np.eye(n)

    grad_f_old = partial_derivatives.grad_f_fd(df, f, x)

    p = grad_f_old

    line_search = LineSearch(f)

    for k in range(k_max):

        p = s*np.dot(H, grad_f_old)

        line_search.update(x, p)

        h_min, h_max = h_interval(line_search.f_line, s)

        h_opt = hybrid.brent(line_search.f_line, h_min, h_max, mode)

        dx = h_opt*p
        x = x + dx

        grad_f_new = partial_derivatives.grad_f_fd(df,f,x)
        dy = grad_f_new - grad_f_old

        norm_g = np.linalg.norm(grad_f_new)
        cor_norm = np.linalg.norm(dx)

        if output:
            print(f'k = {k}, Res Norm: {norm_g:.4e}, Cor Norm: {cor_norm:.4e}')

        converged = norm_g < tol or cor_norm < tol
        if converged:
            break

        H = update_hessian_inverse(H, dx, dy)

        grad_f_old = grad_f_new

    return x

def update_hessian_inverse(H: npt.NDArray[np.float64], dx: npt.NDArray[np.float64],
    dy: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Applies the BFGS update formula to the inverse Hessian matrix."""

    dx_dy = np.dot(dx,dy)

    n = len(dx)
    rho = 1.0/(dx_dy+1e-12)

    I = np.eye(n)

    A1 = I - rho*np.outer(dx, dy)
    A2 = I - rho*np.outer(dy, dx)

    H_new = np.dot(A1, np.dot(H, A2)) + (rho*np.outer(dx, dx))

    return H_new
