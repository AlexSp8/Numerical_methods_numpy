
from typing import Callable

import numpy as np

from linear_systems import direct_solver
from non_linear_systems import newton_solver, non_linear_problem

def linear_model_predictions(xi: np.ndarray[tuple[int, int], np.dtype[np.float64]],
    a: np.ndarray[tuple[int], np.dtype[np.float64]]
    ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """Returns the predictions of a polynomial at the data points
    given the coefficients of the polynomial.

    Args:
        xi[n x m]: the independent variable data points
        a: the coefficients of the polynomial
    Returns:
        The dependent variable predictions of the polynomial at the data points."""

    # X = np.column_stack((np.ones(xi.shape[0]), xi))

    n = xi.shape[0]
    m = xi.shape[1]
    X = np.zeros((n,m+1))
    X[:,0] = np.ones(n)
    X[:,1:] = xi.copy()

    return np.matmul(X,a)

def statistical_quantities(yi: np.ndarray[tuple[int], np.dtype[np.float64]],
    y_model: np.ndarray[tuple[int], np.dtype[np.float64]], m: int
    ) -> tuple[float, float, float]:
    """Returns the statistical quantities of polynomial regression
    given the data points
    and the coefficients of the polynomial regression.

    Args:
        yi: the dependent variable data points
        y_model: the dependent variable polynomial predictions at the data points
        m: the number of independent variables (dimensions)
    Returns:
        The statistical quantities of polynomial regression
        (standard error, coefficient of determination)."""

    n = yi.shape[0]

    y_mean = np.mean(yi)

    Sr, St = np.sum((yi-y_model)**2), np.sum((yi-y_mean)**2)

    std_error = (Sr/(n-(m+1)))**0.5 # standard error: spread around regression line

    r2_coef = ((St-Sr)/St) # coefficient of determination: improvement of error

    ssr = Sr

    return std_error, r2_coef, ssr

def transform_data(xi: np.ndarray[tuple[int], np.dtype[np.float64]],
    f_basis: list[Callable[[np.ndarray[tuple[int], np.dtype[np.float64]]],
                           np.ndarray[tuple[int], np.dtype[np.float64]]]]
    ) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
    """Transforms univariate data to custom basis functions
    for multivariate regression.

    Args:
        xi: independent variable data
        f_basis: a list of basis functions to be called
    Returns:
        The coefficients of linear regression."""

    # columns = [f(xi) for f in f_basis]
    # return np.column_stack(columns)

    nrows = xi.shape[0]
    ncols = len(f_basis)
    xi_mat = np.zeros((nrows, ncols))
    for j in range(ncols):
        xi_mat[:,j] = f_basis[j](xi)

    return xi_mat

def linear_regression(xi: np.ndarray[tuple[int], np.dtype[np.float64]],
    yi: np.ndarray[tuple[int], np.dtype[np.float64]]
    ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """Returns the coefficients of linear regression such that
    the linear function f(x) = a0+a1*x best describes the data.

    Args:
        xi: independent variable data
        yi: dependent variable data
    Returns:
        The coefficients of linear regression."""

    n = xi.shape[0]

    sum_x, sum_y = np.sum(xi), np.sum(yi)
    sum_xy, sum_x2 = np.dot(xi,yi), np.dot(xi,xi)

    x_mean, y_mean = sum_x/n, sum_y/n

    a = np.zeros(2)
    a[1] = (n*sum_xy - sum_x*sum_y)/( n*sum_x2 - sum_x**2 )
    a[0] = y_mean - a[1]*x_mean

    return a

def polynomial_regression(xi: np.ndarray[tuple[int], np.dtype[np.float64]],
    yi: np.ndarray[tuple[int], np.dtype[np.float64]], m: int = 1
    ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """Returns the coefficients of polynomial regression such that
    the linear function f(x) = a0+a1*x+...+am*x best describes the data.

    Args:
        xi: independent variable data
        yi: dependent variable data
        m: the order of the polynomial
    Returns:
        The coefficients of polynomial regression."""

    lu = direct_solver.LUSolver()

    n = xi.shape[0]
    if (n < m+1):
        print('Not enough data points for mth-order polynomial regression!')
        return None

    b = np.zeros(m+1)
    for p in range(m+1):
        b[p] = np.sum(yi*(xi**p))

    sum_x = np.zeros(2*m+1)
    for p in range(2*m+1):
        sum_x[p] = np.sum(xi**p)

    A = np.zeros((m+1,m+1))
    for i in range(m+1):
        for j in range(m+1):
            A[i,j] = sum_x[i+j]

    a_coef = lu.solve(A, b)

    return a_coef

def multi_linear_regression(xi: np.ndarray[tuple[int, int], np.dtype[np.float64]],
    yi: np.ndarray[tuple[int], np.dtype[np.float64]]
    ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """Returns the coefficients of multi-linear regression such that
    the multi-linear function f(x1, x2, ..., xm) = a0+a1*x1+...+am*xm best describes the data.

    Args:
        xi[n x m]: independent variable data
        yi[n]: dependent variable data
    Returns:
        The coefficients of multi-linear regression."""

    n, ndim = xi.shape[0], xi.shape[1]
    if (n < ndim+1):
        print('Not enough data points for multi-linear regression!')
        return None

    # X = np.column_stack((np.ones(n), xi))

    X = np.zeros((n,ndim+1))
    X[:,0] = np.ones(n)
    X[:,1:] = xi.copy()

    A = np.dot(X.T, X)

    b = np.dot(X.T, yi)

    lu = direct_solver.LUSolver()
    a_coef = lu.solve(A, b)

    return a_coef

def non_linear_regression(
    f: Callable[[np.ndarray[tuple[int], np.dtype[np.float64]], np.ndarray[tuple[int], np.dtype[np.float64]]],
                np.ndarray[tuple[int], np.dtype[np.float64]]],
    xi: np.ndarray[tuple[int, int], np.dtype[np.float64]], yi: np.ndarray[tuple[int], np.dtype[np.float64]],
    u0: np.ndarray[tuple[int], np.dtype[np.float64]], output: bool = False, tol: float = 1e-8,
    k_max: int = 1000, r: float = 1.0) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """Returns the coefficients of non-linear regression such that
    the non linear function f(x1, x2, ..., xm) best describes the data.

    Args:
        f: the non-linear function (model)
        xi[n x m]: the independent variables data
        yi[n]: the dependent variable data
        u0: the initial guess for the coefficients of the model
        output (optional): print iteration info
        tol (optional): threshold for convergence
        k_max (optional): maximum number of iterations
        r (optional): relaxation factor for updating the coefficients
    Returns:
        The coefficients of the non-linear model that best describe the data."""

    ls_solver = direct_solver.LUSolver()

    nr_solver = newton_solver.LevenbergMarquardtSolver(ls_solver, u0, k_max, tol, r, l0=1e-4, scale=2.0)

    problem = non_linear_problem.Regression(f=f, xi=xi, yi=yi)

    u = nr_solver.solve(problem, output=output)

    return u
