
from typing import List

from utilities import indexing
from linear_systems import direct_solvers

def vandermonde(xi: List[float] = None, yi: List[float] = None) -> List[float]:
    """Returns the coefficients of the Vandermonde polynomial
    of order n-1, given n data points."""

    n = len(xi)

    A = [ [0.0 for _ in range(n)] for _ in range(n) ]
    for i in range(n):
        for j in range(n):
            A[i][j] = xi[i]**j

    return direct_solvers.lu_decomposition_solve(A,yi)

def newton_gregory_polynomials(xp: List[float], xi: List[float] = None,
    yi: List[float] = None) -> List[List[float]]:
    """Returns the interpolation value at a list of points, xp, for every
    Newton-Gregory polynomials of order 1 to n-1, given n data points."""

    n = len(yi)
    np = len(xp)

    fd = [ [0.0 for _ in range(n)] for _ in range(n) ]
    # 0th order derivatives at data points (xi,yi)
    for i in range(n):
        fd[i][0] = yi[i]
    # derivatives of order 1 to n-1
    for j in range(1,n):
        for i in range(n-j):
            fd[i][j] = ( fd[i+1][j-1] - fd[i][j-1] )/( xi[i+j] - xi[i] )

    y_int = [ [0.0 for _ in range(np)] for _ in range(n) ]
    y_int[0] = [ fd[0][0] for _ in range(np) ]
    x_term = [1.0]*np
    for j in range(np):
        for i in range(1,n):

            x_term[j] *= (xp[j] - xi[i-1])
            y_int[i][j] = y_int[i-1][j] + fd[0][i]*x_term[j]
            # err = y_int[i][j] - y_int[i-1][j]
            # print(err)

    return y_int

def lagrange_polynomial(xp: List[float], xi: List[float] = None,
    yi: List[float] = None, m: int = 1) -> List[float]:
    """Returns the interpolation value at a list of points, xp, using
    Lagrange polynomials of given order m (default 1), given n data points."""

    np = len(xp)

    n = len(xi)
    if n < m + 1:
        m = n - 1

    y_int = [0.0]*np
    for k in range(np):
        i_start = indexing.starting_index(xp[k], xi, m)

        xi_p, yi_p = xi[i_start:i_start+m+1], yi[i_start:i_start+m+1]

        for i in range(m+1):
            p = yi_p[i]
            for j in range(m+1):
                if i != j:
                    p *= (xp[k]-xi_p[j])/(xi_p[i]-xi_p[j])
            y_int[k] += p

    return y_int

def find_spline_interval(xp: float, xi: List[float] = None) -> int:
    """Returns the linear spline interval inside which xp lies."""

    if xp < xi[0]:
        return 0

    n = len(xi)
    if xp > xi[-1]:
        return n-2

    for i in range(n-1):
        if xi[i] <= xp and xp <= xi[i+1]:
            return i

    return 0

def linear_splines(xp: List[float], xi: List[float] = None,
    yi: List[float] = None) -> List[float]:
    """Returns the interpolation value at a list of points, xp, using
    linear splines, given n data points."""

    np = len(xp)
    y_int = [0.0]*np
    for i in range(np):
        idx = find_spline_interval(xp[i], xi)

        m_slope = ( yi[idx+1]-yi[idx] )/( xi[idx+1]-xi[idx] )

        y_int[i] = yi[idx] + m_slope*(xp[i]-xi[idx])
    return y_int

def quadratic_splines(xp: List[float], xi: List[float] = None,
    yi: List[float] = None, bc: str = 'linear', v: float = 0.0) -> List[float]:
    """Returns the interpolation value at a list of points, xp, using
    quadratic splines, given n data points. A boundary condition
    for one endpoint can also be given (default is linear at the first point)"""

    n = len(xi)
    n_int = n-1
    l, d, u, b_vec = [1.0]*n_int, [1.0]*n, [0.0]*n_int, [0.0]*n

    if bc == "linear":
        u[0] = -1.0
    elif bc == "clamped":
        b_vec[0] = v
    elif bc == "average":
        b_vec[0] = (yi[1] - yi[0]) / (xi[1] - xi[0])

    for i in range(n_int):
        h = xi[i+1] - xi[i]
        b_vec[i+1] = 2*( yi[i+1]-yi[i] )/h

    b = direct_solvers.thomas_tri_diagonal(l,d,u,b_vec)
    # b = direct_solvers.lu_decomposition_solve(A,b_vec)

    # c = [0.0]*n_int
    # for i in range(n_int):
    #     h = xi[i+1] - xi[i]
    #     c[i] = ( b[i+1]-b[i] )/(2*h)

    np = len(xp)
    y_int = [0.0]*np
    for i in range(np):

        idx = find_spline_interval(xp[i], xi)

        # ci = c[idx]
        bi = b[idx]
        ci = ( b[idx+1]-bi )/( 2*(xi[idx+1]-xi[idx]) )

        dx = xp[i] - xi[idx]
        y_int[i] = yi[idx] + bi*dx + ci*(dx**2)

    return y_int

def cubic_splines(xp: List[float], xi: List[float] = None,
    yi: List[float] = None, bc: str = 'natural', v: float = 0.0) -> List[float]:
    """Returns the interpolation value at a list of points, xp, using
    cubic splines, given n data points. A boundary condition
    can also be given (default is natural at the endpoints)"""

    n = len(xi)
    n_int = n-1
    l, d, u, b_vec = [0.0]*n_int, [0.0]*n, [0.0]*n_int, [0.0]*n

    for i in range(1, n_int):
        hi = xi[i+1] - xi[i]
        hi_1 = xi[i] - xi[i-1]
        l[i-1] = hi_1
        d[i] = 2*(hi+hi_1)
        u[i] = hi
        b_vec[i] = 3*( (yi[i+1]-yi[i])/hi - (yi[i]-yi[i-1])/hi_1 )

    if bc == "natural":
        d[0], d[-1] = 1.0, 1.0
        u[0] = 0.0
        b_vec[0], b_vec[-1] = 0.0, 0.0
    elif bc == 'clamped':
        h1 = xi[1]-xi[0]
        d[0], u[0] = 2.0, 1.0
        b_vec[0] = (3/h1)*( (yi[1]-yi[0])/h1 - v )

        hn = xi[-1]-xi[-2]
        d[-1], l[-1] = 2.0, 1.0
        b_vec[-1] = (3/hn)*( v -(yi[-1]-yi[-2])/hn )

    c = direct_solvers.thomas_tri_diagonal(l,d,u,b_vec)
    # c = direct_solvers.lu_decomposition_solve(A,b_vec)

    np = len(xp)
    y_int = [0.0]*np
    for i in range(np):
        idx = find_spline_interval(xp[i], xi)

        ci = c[idx]
        hi = xi[idx+1]-xi[idx]
        bi = ( (yi[idx+1]-yi[idx])/hi - hi*(c[idx+1]+2*ci)/3 )
        di = (c[idx+1]-ci)/(3*hi)

        dx = xp[i] - xi[idx]
        y_int[i] = yi[idx] + bi*dx + ci*(dx**2) + di*(dx**3)
    return y_int

def multi_lagrange(xp: List[float], xi: List[List[float]] = None,
    yi: List[List[float]] = None, m: int = 1) -> float:
    """Returns the interpolation value at a list of points, xp, in multiple
    dimensions using Lagrange polynomials of given order m, given n data points."""

    n = len(yi)

    if len(xp) == 1:
        return lagrange_polynomial(xp, xi[0], yi, m)[0]

    temp_y = [0.0]*n
    for i in range(n):
        temp_y[i] = multi_lagrange(xp[1:], xi[1:], yi[i], m)

    return lagrange_polynomial([xp[0]], xi[0], temp_y, m)[0]

def multi_cubic_spline(xp: List[float], xi: List[List[float]] = None,
    yi: List[List[float]] = None, bc: str = 'natural', v: float = 0.0) -> float:
    """Returns the interpolation value at a list of points, xp, in multiple
    dimensions using Cubic splines, given n data points."""

    n = len(yi)

    if len(xp) == 1:
        return cubic_splines(xp, xi[0], yi, bc, v)[0]

    temp_y = [0.0]*n
    for i in range(n):
        temp_y[i] = multi_cubic_spline(xp[1:], xi[1:], yi[i], bc, v)

    return cubic_splines([xp[0]], xi[0], temp_y, bc, v)[0]

def multi_interpolation(xp: List[float], xi: List[List[float]], yi: List[List[float]],
    method_1d, **kwargs) -> float:
    """Returns the interpolation value at a list of points, xp, in multiple
    dimensions using given method_1d, given n data points."""

    if len(xp) == 1:
        return method_1d([xp[0]], xi[0], yi, **kwargs)[0]

    intermediate_yi = []
    for sub_slice in yi:
        res = multi_interpolation(xp[1:], xi[1:], sub_slice, method_1d, **kwargs)
        intermediate_yi.append(res)

    return method_1d([xp[0]], xi[0], intermediate_yi, **kwargs)[0]
