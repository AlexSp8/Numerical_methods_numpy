
import numpy as np
import numpy.typing as npt

from utilities import indexing
from linear_systems import direct_solver

def vandermonde(xi: npt.NDArray, yi: npt.NDArray) -> npt.NDArray:
    """Returns the coefficients of the Vandermonde polynomial
    of order n-1, given n data points."""

    n = xi.shape[0]

    # A = np.vander(xi, n, increasing=True)

    # powers = np.arange(n)
    # A = xi[:, np.newaxis] ** powers

    A = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            A[i][j] = xi[i]**j

    lu = direct_solver.LUSolver()
    return lu.solve(A,yi)

def newton_gregory_polynomials(xp: npt.NDArray, xi: npt.NDArray,
    yi: npt.NDArray) -> npt.NDArray:
    """Returns the interpolation value at a list of points, xp, for every
    Newton-Gregory polynomial of order 1 to n-1, given n data points."""

    n = yi.shape[0]
    npts = xp.shape[0]

    fd = np.zeros((n,n))
    # 0th order derivatives at data points (xi,yi)
    fd[:,0] = yi[:]

    # derivatives of order 1 to n-1
    for j in range(1,n):
        fd[:n-j,j] = (fd[1:n-j+1,j-1] - fd[:n-j,j-1])/( xi[j:] - xi[:n-j] )

    y_int = np.zeros((n,npts))
    y_int[0,:] = fd[0,0]

    x_term = np.ones(npts)
    for i in range(1,n):
        x_term *= (xp - xi[i-1])
        y_int[i,:] = y_int[i-1,:] + fd[0,i]*x_term
        # err = y_int[i,:] - y_int[i-1,:]
        # print(err)

    return y_int

def lagrange_polynomial(xp: npt.NDArray, xi: npt.NDArray,
    yi: npt.NDArray, m: int = 1) -> npt.NDArray:
    """Returns the interpolation value at a list of points, xp, using
    Lagrange polynomials of given order m (default 1), given n data points."""

    npts = xp.shape[0]

    n = len(xi)
    if n < m + 1:
        m = n - 1

    y_int = np.zeros(npts)
    for k in range(npts):
        i_start = indexing.starting_index(xp[k], xi, m)

        xi_p, yi_p = xi[i_start : i_start+m+1], yi[i_start : i_start+m+1]

        L = np.ones(m+1)
        for i in range(m+1):
            xj = np.delete(xi_p, i)
            L[i] = np.prod( (xp[k]-xj)/(xi_p[i]-xj) )

        y_int[k] = np.dot(yi_p, L)

    return y_int

def find_spline_interval(xp: float, xi: npt.NDArray) -> int:
    """Returns the linear spline interval inside which xp lies."""

    if xp < xi[0]:
        return 0

    n = len(xi)
    if xp > xi[-1]:
        return n-2

    return np.searchsorted(xi, xp) - 1

def linear_splines(xp: npt.NDArray, xi: npt.NDArray,
    yi: npt.NDArray) -> npt.NDArray:
    """Returns the interpolation value at a list of points, xp, using
    linear splines, given n data points."""

    # npts = xp.shape[0]

    # y_int = np.zeros(npts)
    # for i in range(npts):
    #     idx = find_spline_interval(xp[i], xi)

    #     m_slope = ( yi[idx+1]-yi[idx] )/( xi[idx+1]-xi[idx] )

    #     y_int[i] = yi[idx] + m_slope*(xp[i]-xi[idx])
    # return y_int

    idx = np.searchsorted(xi, xp) - 1
    idx = np.clip(idx, 0, len(xi) - 2)

    x0, x1 = xi[idx], xi[idx+1]
    y0, y1 = yi[idx], yi[idx+1]

    m_slope = (y1 - y0)/(x1 - x0)
    y_int = y0 + m_slope*(xp - x0)

    return y_int

def quadratic_splines(xp: npt.NDArray, xi: npt.NDArray,
    yi: npt.NDArray, bc: str = 'linear', v: float = 0.0) -> npt.NDArray:
    """Returns the interpolation value at a list of points, xp, using
    quadratic splines, given n data points. A boundary condition
    for one endpoint can also be given (default is linear at the first point)"""

    n = len(xi)
    n_int = n-1

    l = np.ones(n_int)
    d = np.ones(n)
    u = np.zeros(n_int)
    b_vec = np.zeros(n)

    if bc == "linear":
        u[0] = -1.0
    elif bc == "clamped":
        b_vec[0] = v
    elif bc == "average":
        b_vec[0] = (yi[1] - yi[0])/(xi[1] - xi[0])

    hi = np.diff(xi)
    b_vec[1:] = 2*np.diff(yi)/hi

    b = direct_solver.thomas_tri_diagonal(l,d,u,b_vec)
    # lu = direct_solver.LUSolver()
    # b = lu.solve(A,b_vec)

    c = np.diff(b)/(2*hi)

    idx = np.searchsorted(xi, xp) - 1
    idx = np.clip(idx, 0, n_int - 1)

    bi = b[idx]
    ci = c[idx]
    xi_val = xi[idx]
    yi_val = yi[idx]

    dx = xp - xi_val
    y_int = yi_val + bi*dx + ci*(dx**2)

    return y_int

def cubic_splines(xp: npt.NDArray, xi: npt.NDArray,
    yi: npt.NDArray, bc: str = 'natural', v: float = 0.0) -> npt.NDArray:
    """Returns the interpolation value at a list of points, xp, using
    cubic splines, given n data points. A boundary condition
    can also be given (default is natural at the endpoints)"""

    n = len(xi)
    n_int = n-1

    l = np.ones(n_int)
    d = np.ones(n)
    u = np.zeros(n_int)
    b_vec = np.zeros(n)

    hi = np.diff(xi)
    l[:n_int-1] = hi[:-1]
    u[1:] = hi[1:]
    d[1:n_int] = 2*(hi[:-1]+hi[1:])

    dydx = np.diff(yi)/hi

    b_vec[1:n_int] = 3*(dydx[1:]-dydx[:-1])

    if bc == "natural":
        d[0], d[-1] = 1.0, 1.0
        u[0] = 0.0
        b_vec[0], b_vec[-1] = 0.0, 0.0
    elif bc == 'clamped':
        d[0], u[0] = 2.0, 1.0
        b_vec[0] = (3/hi[0])*( dydx[0] - v )

        d[-1], l[-1] = 2.0, 1.0
        b_vec[-1] = (3/hi[-1])*( v -dydx[-1] )

    ci = direct_solver.thomas_tri_diagonal(l,d,u,b_vec)
    # lu = direct_solver.LUSolver()
    # c = lu.solve(A,b_vec)

    bi = dydx - hi*(ci[1:] + 2*ci[:-1])/3
    di = (ci[1:] - ci[:-1])/(3*hi)

    idx = np.searchsorted(xi, xp) - 1
    idx = np.clip(idx, 0, n_int - 1)

    dx = xp - xi[idx]

    y_int = yi[idx] + bi[idx]*dx + ci[idx]*(dx**2) + di[idx]*(dx**3)

    return y_int
