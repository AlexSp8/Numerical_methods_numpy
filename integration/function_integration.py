
from typing import Callable
import numpy as np
import numpy.typing as npt

from integration import gauss_legendre_data
from non_linear_systems import newton_solver, non_linear_problem
from linear_systems import direct_solver

def trapezoidal(f: Callable[[float], float],
    a: float, b: float, n: int = 1) -> float:
    """Returns the value of the integral of function, f, in [a,b]
    using the trapezoidal rule across n segments."""

    if n <= 0:
        raise ValueError("n segments must be positive.")

    h = (b-a)/n

    s = 0.0
    for i in range(1,n):
        x = a+i*h
        s += f(x)

    s *= 2

    return 0.5*h*(f(a)+s+f(b))

def simpson_13(f: Callable[[float], float],
    a: float, b: float, n: int = 2) -> float:
    """Returns the value of the integral of function, f, in [a,b]
    using the Simpson 1/3 rule across n segments."""

    if n <= 0:
        raise ValueError("n segments must be positive.")

    if n%2 != 0:
        raise ValueError("Even number of segments required for Simpson 1/3 integration")

    h = (b-a)/n

    s = 0.0
    for i in range(0,n,2):
        x0 = a+i*h
        s += f(x0) + 4*f(x0+h) + f(x0+2*h)

    return (h/3)*s

def simpson_38(f: Callable[[float], float], a: float, b: float) -> float:
    """Returns the value of the integral of function, f, in [a,b]
    using the Simpson 3/8 rule across n segments."""

    n = 3

    h = (b-a)/n
    return (b-a)*( f(a) + 3*(f(a+h) + f(b-h)) + f(b))/8

def simpson_mixed(f: Callable[[float], float],
    a: float, b: float, n: int = 1) -> float:
    """Returns the value of the integral of function, f, in [a,b]
    using mixed Simpson rules across n segments.
    For even segments Simpson's 1/3 is used.
    For odd segments Simpson's 3/8 is used for the last 3 segments only."""

    if n <= 0:
        raise ValueError("n segments must be positive.")

    if n == 1:
        return trapezoidal(f, a, b, n)

    if n == 2:
        return simpson_13(f, a, b, n)

    h = (b-a)/n
    m = n
    s = 0.0
    if m%2 == 1:
        s += simpson_38(f, b-3*h, b)
        m = n-3

    if m > 0:
        s += simpson_13(f, a, a+m*h, m)

    return s

def romberg(f: Callable[[float], float], a: float, b: float,
    k_max: int = 10, eps: float = 1e-8) -> float:
    """Returns the value of the integral of function, f, in [a,b]
    using the Romberg's formula."""

    s = np.zeros((k_max+1,k_max+1))

    s[0,0] = trapezoidal(f, a, b, n=1)

    for k in range(1, k_max+1):

        n = 2**k
        s[k,0] = trapezoidal(f, a, b, n)

        for j in range(1,k+1):
            w = 4**j
            s[k,j] = ( w*s[k,j-1] - s[k-1,j-1] )/(w-1)

        err = abs( (s[k,k]-s[k-1,k-1])/s[k,k] )
        if err < eps:
            # print(f"Romberg converged at k = {k}")
            break

    return s[k,k]

def gauss_legendre_1D(f: Callable[[float], float],
    a: float, b: float, n: int = 1) -> float:
    """Returns the value of the integral of function, f, in [a,b]
    using n-point Gauss-Legendre quadrature. The quadrature accurately
    fits a polynomial of order 2n-1."""

    if n <= 0:
        raise ValueError("n points must be positive.")

    c1 = (b-a)/2
    d1 = (b+a)/2

    c, w = gauss_legendre_coefficients(n)
    # c = gauss_legendre_data.GAUSS_POINTS[n]
    # w = gauss_legendre_data.GAUSS_WEIGHTS[n]

    s = 0.0
    for i in range(n):
        xi = c1*c[i] + d1
        s += w[i]*f(xi)

    return c1*s

def gauss_legendre_2D(f: Callable[[npt.NDArray], float],
    a: npt.NDArray, b: npt.NDArray, n: npt.NDArray = np.array([1, 1])) -> float:
    """Returns the value of the integral of a 2D function, f(x,y),
    using n-point Gauss-Legendre quadrature in each direction.
    The limits a[i], b[i] and the number of points n[i] in each
    direction must be given."""

    if any(n) <= 0:
        raise ValueError("n points must be positive.")

    cx, wx = gauss_legendre_coefficients(n[0])
    cy, wy = gauss_legendre_coefficients(n[1])

    ndim = len(a)
    jac, d = np.zeros(ndim), np.zeros(ndim)
    jac_pr = 1.0
    for i in range(ndim):
        jac[i] = (b[i]-a[i])/2
        d[i] = (b[i]+a[i])/2
        jac_pr *= jac[i]

    s = 0.0
    for i in range(n[0]):
        xi = jac[0]*cx[i] + d[0]
        for j in range(n[1]):
            yj = jac[1]*cy[j] + d[1]
            x = [xi, yj]
            s += wx[i]*wy[j]*f(x)

    return s*jac_pr

def gauss_legendre_nD_recursive(f: Callable[[npt.NDArray], float],
    a: npt.NDArray, b: npt.NDArray, n: npt.NDArray, dim: int = 0,
    xp: npt.NDArray = None) -> float:
    """Returns the value of the integral of a nD function, f(x,y),
    using n-point Gauss-Legendre quadrature in each direction.
    The limits a[i], b[i] and the number of points n[i] in each
    direction must be given."""

    ndim = len(a)

    if xp is None:
        xp = np.zeros(ndim)

    if dim == ndim:
        return f(xp)

    c, w = gauss_legendre_coefficients(n[dim])

    jac = (b[dim] - a[dim])/2.0
    d = (b[dim] + a[dim])/2.0

    s = 0.0
    for i in range(n[dim]):
        xp[dim] = jac*c[i] + d
        val = gauss_legendre_nD_recursive(f, a, b, n, dim+1, xp)
        s += w[i]*val

    return s*jac

def gauss_legendre_nD(f: Callable[[npt.NDArray], float],
    a: npt.NDArray, b: npt.NDArray, n:npt.NDArray) -> float:
    """Returns the value of the integral of a nD function, f(x,y),
    using n-point Gauss-Legendre quadrature in each direction.
    The limits a[i], b[i] and the number of points n[i] in each
    direction must be given."""

    import itertools

    ndim = len(a)

    weights = [0.0]*ndim
    points = [0.0]*ndim
    jac = np.zeros(ndim) # Jacobians
    d = np.zeros(ndim) # Offsets

    jac_pr = 1.0
    for i in range(ndim):
        c, w = gauss_legendre_coefficients(n[i])
        weights[i] = w
        points[i] = c
        jac[i] = (b[i] - a[i])/2.0
        d[i] = (b[i] + a[i])/2.0
        jac_pr *= jac[i]

    s = 0.0
    point = np.zeros(ndim)
    for indices in itertools.product(*(range(ni) for ni in n)):
        w = 1.0
        for dim, idx in enumerate(indices):
            point[dim] = jac[dim]*points[dim][idx] + d[dim]
            w *= weights[dim][idx]

        s += w*f(point)

    return s*jac_pr

def legendre_polynomial(x: float, n: int) -> tuple[float, float]:
    """Evaluates the n-th Legendre polynomial and its derivative at x
    using the recursive relation."""

    if n == 0:
        return 1.0, 0.0

    p_prev = 0.0
    p = 1.0

    for i in range(n):
        p_next = ((2.0*i + 1.0)*x*p - i*p_prev)/(i + 1.0)
        p_prev = p
        p = p_next

    dp = n*(x*p - p_prev)/(x**2 - 1)

    return p, dp

def gauss_legendre_coefficients(n: int = 1) -> tuple[npt.NDArray, npt.NDArray]:
    """Solves for Gauss-Legendre points and weights of degree n."""

    points = np.zeros(n)
    weights = np.zeros(n)

    # We only need to find half the roots due to symmetry
    m = (n + 1)//2

    for i in range(m):
        # Initial guess for the i-th root (Chebyshev-like spacing)
        x = np.cos(np.pi*(i + 0.75)/(n + 0.5))

        for _ in range(100):
            p, dp = legendre_polynomial(x, n)
            dx = p/dp
            x -= dx
            if abs(dx) < 1e-12:
                break

        points[i] = -x
        points[n-1-i] = x

        weight = 2.0/((1.0 - x**2)*(dp**2))
        weights[i] = weight
        weights[n-1-i] = weight

    return points, weights

def gauss_legendre_residuals(x: npt.NDArray) -> npt.NDArray:
    """Returns the residuals of the Gauss-Legendre quadrature equations
    for n-point integration."""

    n = int(len(x)/2)

    res = np.zeros(2*n)
    for i in range(2*n):
        s = 0.0
        for j in range(n):
            w, c = x[2*j], x[2*j+1]
            s += w*(c**i)
        v = (1 - (-1)**(i+1))/(i+1)
        res[i] = s - v
    return res

def solve_gauss_legendre_system(n: int = 1) -> tuple[npt.NDArray, npt.NDArray]:
    """Returns the points and weights of Gauss-Legendre quadrature
    for n-point integration."""

    if n <= 0:
        raise ValueError("n points must be positive.")

    x = np.zeros(2*n)
    for i in range(n):
        x[2*i] = (2.0/n) # Guess for weight
        x[2*i+1] = (-1.0 + (2.0 * i + 1.0) / n) # Guess for point

    ls_solver = direct_solver.LUSolver()
    nr_solver = newton_solver.NewtonSolver(ls_solver,x)
    problem = non_linear_problem.Regular(gauss_legendre_residuals)
    x = nr_solver.solve(problem)

    w, c = np.zeros(n), np.zeros(n)
    for i in range(n):
        w[i], c[i] = x[2*i], x[2*i+1]

    return c, w
