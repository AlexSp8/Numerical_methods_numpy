
from typing import Callable, List

def grad_f_fd(df: Callable[[Callable[[float], float], float, float], float],
    f: Callable[[List[float]], float], x: List[float],
    h: float = 1e-8) -> List[float]:
    """Returns the gradient of an n-D function, f(x),
    at a point, x, given the finite difference method."""
    n, f_x = len(x), f(x)
    grad_f = [0.0]*n
    for i in range(n):
        fi = f_1d_scalar(f, x, i)
        grad_f[i] = df(fi, x[i], f_x, h)
    return grad_f

def f_1d_scalar(f_nd_scalar: Callable[[List[float]], float],
    x: List[float], i: int) -> Callable[[float], float]:
    """Returns a 1-D scalar function, f(xi), from an n-D scalar function f(x)."""

    def f_1d(xi: float) -> float:
        x_new = x[:]
        x_new[i] = xi
        return f_nd_scalar(x_new)

    return f_1d

def hessian_f_fd(df: Callable[[Callable[[float], float], float, float], float],
    f: Callable[[List[float]], float], x: List[float],
    h: float = 1e-4) -> List[List[float]]:
    """Returns the hessian of an n-D function, f(x),
    at a point, x, given the finite difference method."""

    n = len(x)
    hessian = [None]*n

    def full_gradient(x_vec):
        return grad_f_fd(df, f, x_vec, h)

    for i in range(n):
        grad_fi = f_nd_scalar(full_gradient, i)
        row = grad_f_fd(df, grad_fi, x, h)
        hessian[i] = row

    return hessian

def f_nd_scalar(f_vec: Callable[[List[float]], List[float]],
    i: int) -> Callable[[List[float]], float]:
    """Returns an n-D scalar function, f(x)[i], from an n-D vector function, f_vec(x)."""
    def f_scalar(x: List[float]) -> float:
        return f_vec(x)[i]
    return f_scalar
