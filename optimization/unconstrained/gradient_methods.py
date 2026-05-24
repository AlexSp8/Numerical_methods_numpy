
from typing import Callable
import numpy as np
import numpy.typing as npt
from scipy.optimize import minimize

from differentiation import partial_derivatives
from optimization.one_dimensional import hybrid
from optimization.unconstrained.line_search import LineSearch
from linear_systems import direct_solver
from non_linear_systems import newton_solver, non_linear_problem

def steepest_descent(f: Callable[[np.ndarray[tuple[int], np.dtype[np.float64]]], float],
    x0: np.ndarray[tuple[int], np.dtype[np.float64]],
    df: Callable[[Callable[[float], float], float, float, float], float],
    tol: float = 1e-8, k_max = 100) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """Returns the point of minimum of a multi-variable function using the
    Steepest Descent method.

    Args:
        f: the function to be optimized
        x0: initial guess
        df: the function to calculate the derivative of f
        tol (optional): threshold for convergence
        k_max (optional): maximum iterations
    Returns:
        The point of minimum of f."""

    x = x0.copy()

    line_search = LineSearch(f)

    for k in range(1, k_max+1):

        grad_f = partial_derivatives.grad_f_fd(df, f, x)

        norm_g = np.linalg.norm(grad_f)
        if norm_g < tol:
            print(f'Converged at k = {k}')
            return x

        d = -grad_f/norm_g

        line_search.update(x, d)

        h_min, h_max = line_search.h_interval(k_max=10)

        h_opt = hybrid.brent(line_search.f_line, h_min, h_max)

        x = x + h_opt*d

    print(f"Warning: Maximum iterations ({k_max}) reached without converging.")

    return x

def conjugate_gradient(f: Callable[[np.ndarray[tuple[int], np.dtype[np.float64]]], float],
    x0: np.ndarray[tuple[int], np.dtype[np.float64]],
    df: Callable[[Callable[[float], float], float, float, float], float],
    b_update: str = 'polak_ribiere', tol: float = 1e-8,
    k_max = 1000) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """Returns the point of minimum of a multi-variable function using the
    Conjugate Gradient method.

    Args:
        f: the function to be optimized
        x0: initial guess
        df: the function to calculate the derivative of f
        tol (optional): threshold for convergence
        k_max (optional): maximum iterations
    Returns:
        The point of minimum of f."""

    # res = minimize(f, x0, method='CG')
    # return res.x

    x = x0.copy()

    grad_f = partial_derivatives.grad_f_fd(df,f,x)

    p = -grad_f

    line_search = LineSearch(f)

    for k in range(1, k_max+1):

        norm_p = np.linalg.norm(p)
        if norm_p < tol:
            print(f'Converged at k = {k}')
            return x

        # p = p/norm_p

        line_search.update(x, p)

        h_min, h_max = line_search.h_interval(k_max=10)

        h_opt = hybrid.brent(line_search.f_line, h_min, h_max)

        x = x + h_opt*p

        grad_f_new = partial_derivatives.grad_f_fd(df,f,x)

        if b_update == 'polak_ribiere':
            num = np.dot(grad_f_new, grad_f_new-grad_f)
            den = np.dot(grad_f, grad_f)
        elif b_update == 'fletcher_reeves':
            num = np.dot(grad_f_new, grad_f_new)
            den = np.dot(grad_f, grad_f)
        else:
            raise ValueError('Wrong b_update choice in CG!')

        if den < 1e-12:
            beta = 0.0
        else:
            beta = max(0.0, num/den)

        p = -grad_f_new + beta*p

        grad_f = grad_f_new.copy()

    print(f"Warning: Maximum iterations ({k_max}) reached without converging.")

    return x

def newton(f: Callable[[npt.NDArray[np.float64]], float],
    x0: npt.NDArray[np.float64],
    df: Callable[[Callable[[float], float], float, float, float], float]
    ) -> npt.NDArray[np.float64]:
    """Returns the point of minimum of a multi-variable function using the
    Newton's method.

    Args:
        f: the function to be optimized
        x0: initial guess
        df: the function to calculate the derivative of f
    Returns:
        The point of minimum of f."""

    grad_f = partial_derivatives.grad_f_fd
    hessian_f = partial_derivatives.hessian_f_fd

    # jac_f = lambda x: grad_f(df, f, x)
    # hess_f = lambda x: hessian_f(df, f, x)
    # res = minimize(f, x0, method='Newton-CG', jac=jac_f, hess=hess_f)
    # return res.x

    ls_solver = direct_solver.LUSolver()

    nr_solver = newton_solver.NewtonSolver(ls_solver, u0=x0, k_max=100, tol=1e-8, r=1.0)

    problem = non_linear_problem.Optimization(f=f, df=df, grad_f=grad_f, hessian_f=hessian_f)

    x = nr_solver.solve(problem, output=True)

    return x

def marquardt(f: Callable[[npt.NDArray[np.float64]], float],
    x0: npt.NDArray[np.float64],
    df: Callable[[Callable[[float], float], float, float, float], float]
    ) -> npt.NDArray[np.float64]:
    """Returns the point of minimum of a multi-variable function using the
    Marquardt's method.

    Args:
        f: the function to be optimized
        x0: initial guess
        df: the function to calculate the derivative of f
    Returns:
        The point of minimum of f."""

    ls_solver = direct_solver.LUSolver()

    nr_solver = newton_solver.LevenbergMarquardtSolver(ls_solver=ls_solver, u0=x0,
                                    k_max=100, tol=1e-8, r=1.0, l0=1e-4, scale=10.0)

    grad_f = partial_derivatives.grad_f_fd
    hessian_f = partial_derivatives.hessian_f_fd

    problem = non_linear_problem.Optimization(f=f, df=df, grad_f=grad_f, hessian_f=hessian_f)

    x = nr_solver.solve(problem, output=True)

    return x

def bfgs(f: Callable[[npt.NDArray[np.float64]], float],
    x0: npt.NDArray[np.float64],
    df: Callable[[Callable[[float], float], float, float, float], float],
    output: bool = False, tol: float = 1e-8, k_max = 1000) -> npt.NDArray[np.float64]:
    """Returns the point of minimum of a multi-variable function using the
    BFGS method.

    Args:
        f: the function to be optimized
        x0: initial guess
        df: the function to calculate the derivative of f
        output (optional): print iteration info
        tol (optional): threshold for convergence
        k_max (optional): maximum iterations
    Returns:
        The point of minimum of f."""

    # res = minimize(f, x0, method='BFGS')
    # return res.x

    n = len(x0)
    x = x0.copy()

    invH = np.eye(n)

    grad_f = partial_derivatives.grad_f_fd(df, f, x)

    line_search = LineSearch(f)

    for k in range(1, k_max+1):

        p = -np.dot(invH, grad_f)

        line_search.update(x, p)

        h_min, h_max = line_search.h_interval(k_max=10)

        h_opt = hybrid.brent(line_search.f_line, h_min, h_max)

        dx = h_opt*p

        x = x + dx

        grad_f_new = partial_derivatives.grad_f_fd(df, f, x)

        dy = grad_f_new - grad_f

        norm_g = np.linalg.norm(grad_f_new)
        cor_norm = np.linalg.norm(dx)

        if output:
            print(f'k = {k}, Res Norm: {norm_g:.4e}, Cor Norm: {cor_norm:.4e}')

        if norm_g < tol or cor_norm < tol:
            return x

        invH = update_hessian_inverse(invH, dx, dy)

        grad_f = grad_f_new

    print(f"Warning: Maximum iterations ({k_max}) reached without converging.")

    return x

def update_hessian_inverse(invH: npt.NDArray[np.float64], dx: npt.NDArray[np.float64],
    dy: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Applies the BFGS update formula to the inverse Hessian matrix.

    Args:
        invH: the Hessian inverse
        dx: the solution difference
        dy: the gradient difference
    Returns:
        The new Hessian inverse."""

    dx_dy = np.dot(dx, dy)

    if dx_dy <= 1e-10:
        return invH

    n = len(dx)

    rho = 1.0/dx_dy

    I = np.eye(n)

    A1 = I - rho*np.outer(dx, dy)
    A2 = I - rho*np.outer(dy, dx)

    H_new = np.dot(A1, np.dot(invH, A2)) + (rho*np.outer(dx, dx))

    return H_new

def l_bfgs(f: Callable[[npt.NDArray[np.float64]], float],
    x0: npt.NDArray[np.float64],
    df: Callable[[Callable[[float], float], float, float, float], float],
    output: bool = False, m_max: int = 10, tol: float = 1e-6,
    k_max: int = 1000) -> npt.NDArray[np.float64]:
    """Returns the point of minimum of a multi-variable function using the
    L-BFGS method.

    Args:
        f: the function to be optimized
        x0: initial guess
        df: the function to calculate the derivative of f
        output (optional): print iteration info
        m (optional): the maximum number of updates that are stored
        tol (optional): threshold for convergence
        k_max (optional): maximum iterations
    Returns:
        The point of minimum of f."""

    x = x0.copy()

    grad_f = partial_derivatives.grad_f_fd(df, f, x)

    s = []
    y = []
    rho = []

    for k in range(1, k_max+1):

        m = len(s)

        q = grad_f.copy()

        a = []

        for i in range(m-1, -1, -1):
            alpha = rho[i]*np.dot(s[i], q)
            q = q - alpha*y[i]
            a.append(alpha)

        gamma = 1.0
        if m > 0:
            s_m = s[-1]
            y_m = y[-1]
            gamma = np.dot(s_m, y_m)/np.dot(y_m, y_m)

        z = gamma*q

        a.reverse()

        for i in range(m):
            beta = rho[i]*np.dot(y[i], z)
            z = z + s[i]*(a[i] - beta)

        p = -z

        h_opt = strong_wolfe_line_search(f, df, x, p)

        dx = h_opt*p

        x = x + dx

        grad_f_new = partial_derivatives.grad_f_fd(df, f, x)

        s_k = dx
        y_k = grad_f_new - grad_f
        sy = np.dot(s_k, y_k)

        norm_g = np.linalg.norm(grad_f_new)
        cor_norm = np.linalg.norm(dx)

        if output:
            print(f'k = {k}, Res Norm: {norm_g:.4e}, Cor Norm: {cor_norm:.4e}')

        if norm_g < tol and cor_norm < tol:
            return x

        if sy > 1e-10:

            s.append(s_k)
            y.append(y_k)
            rho.append(1.0/sy)

            if len(s) > m_max:
                # s[:-1] = s[1:]
                # s[-1] = s_k
                # y[:-1] = y[1:]
                # y[-1] = y_k
                # rho[:-1] = rho[1:]
                # rho[-1] = 1.0/sy
                s.pop(0)
                y.pop(0)
                rho.pop(0)

        grad_f = grad_f_new

    print(f"Warning: Maximum iterations ({k_max}) reached without converging.")

    return x

def strong_wolfe_line_search(f: Callable[[npt.NDArray[np.float64]], float],
    df: Callable[[Callable[[float], float], float, float, float], float],
    x: npt.NDArray[np.float64], p: npt.NDArray[np.float64],
    c1: float = 1e-4, c2: float = 0.9,
    h_max: float = 10.0, k_max: int = 10) -> float:
    """Returns the optimal step to minimize a function along a direction
    with a line search that enforces positive curvature through the
    Strong Wolfe conditions (Armijo and Strong curvature).

    Args:
        f: the function to be minimized
        df: the function to calculate the derivative of f
        x: the current point
        p: the line search direction
        c1 (optional): the constant for the Armijo condition
        c2 (optional): the constant for the Strong curvature condition
        h_max (optional): the maximum step size
        k_max (optional): the maximum number of iterations
    Returns:
        The optimal step for minimization of f along the direction, d."""

    h_old = 0.0
    h = 1.0

    grad_f0 = partial_derivatives.grad_f_fd(df, f, x)

    grad_f0_p = np.dot(grad_f0, p)

    # Abort if p is not a descent direction (L-BFGS lost positive definiteness)
    if grad_f0_p >= 0:
        raise ValueError("Search direction is not a descent direction!")

    f0 = f(x)

    f_x = f0

    for k in range(1, k_max+1):

        x_new = x + h*p
        f_x_new = f(x_new)

        # Armijo condition fails or function increases: the minimum is bracketed
        if (f_x_new > f0 + c1*h*grad_f0_p) or (k > 1 and f_x_new >= f_x):
            return zoom(f, df, x, p, h_old, h, f0, grad_f0_p, c1, c2)

        grad_f = partial_derivatives.grad_f_fd(df, f, x_new)
        grad_f_p = np.dot(grad_f, p)

        # Strong Curvature condition
        if abs(grad_f_p) <= -c2*grad_f0_p:
            return h

        # grad flipped
        if grad_f_p >= 0:
            return zoom(f, df, x, p, h, h_old, f0, grad_f0_p, c1, c2)

        h_old = h
        h = min(h*2.0, h_max)

        # if h == h_old:
        #     break

        f_x = f_x_new

    return h

def zoom(f: Callable[[npt.NDArray[np.float64]], float],
    df: Callable[[Callable[[float], float], float, float, float], float],
    x: npt.NDArray[np.float64], p: npt.NDArray[np.float64],
    hl: float, hu: float, f0: float, grad_f0_p: float,
    c1: float = 1e-4, c2: float = 0.9, k_max: int = 20) -> float:
    """Returns the optimal step value to minimize a function along a direction
    with a line search that enforces positive curvature through the
    Strong Wolfe conditions (Armijo and Strong curvature).
    The value is bracketed in an interval which shrinks with bisection.

    Args:
        f: the function to be minimized
        df: the function to calculate the derivative of f
        x: the current point
        p: the line search direction
        hl: the lower limit of the interval
        hu: the upper limit of the interval
        f0: the initial function value
        grad_f0_p: the initial grad value
        c1 (optional): the constant for the Armijo condition
        c2 (optional): the constant for the Strong curvature condition
        h_max (optional): the maximum step size
        k_max (optional): the maximum number of iterations
    Returns:
        The optimal step for minimization of f along the direction, d."""

    for k in range(1, k_max+1):

        h = (hl + hu)/2.0

        x_new = x + h*p
        f_x_new = f(x_new)

        x_l = x + hl*p
        f_x_l = f(x_l)

        if (f_x_new > f0 + c1*h*grad_f0_p) or (f_x_new >= f_x_l):

            # Armijo condition fails, or f increases from the lower bound,
            # the minimum is to the left
            hu = h

        else:

            grad_f = partial_derivatives.grad_f_fd(df, f, x_new)
            grad_f_p = np.dot(grad_f, p)

            # Strong Curvature condition
            if abs(grad_f_p) <= -c2*grad_f0_p:
                return h

            # grad flipped, the minimum is outside the interval (to the left)
            if grad_f_p*(hu - hl) >= 0:
                hu = hl

            # Otherwise, the minimum is to the right.
            hl = h

    # if max_iter is reached
    return (hl + hu) / 2.0
