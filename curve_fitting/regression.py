
from typing import List, Callable

import numpy as np
import numpy.typing as npt

from linear_systems import direct_solver
from non_linear_systems import newton_solver, non_linear_problem

def linear_model_predictions(xi: npt.NDArray[np.float64],
    a: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Returns the predictions of polynomial at the data points
    given the coefficients, a, of the polynomial."""

    X = np.column_stack((np.ones(xi.shape[0]), xi))

    return np.matmul(X,a)

def statistical_quantities(n: int, m: int,
    yi: npt.NDArray[np.float64], y_model: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
    """Returns the statistical quantities of polynomial regression
    (standard error, coefficient of determination) given the data points
    and the coefficients, a, of the polynomial regression."""

    y_mean = np.mean(yi)

    Sr, St = np.sum((yi-y_model)**2), np.sum((yi-y_mean)**2)

    std_error = (Sr/(n-(m+1)))**0.5 # standard error: spread around regression line

    r2_coef = ((St-Sr)/St) # coefficient of determination: improvement of error

    ssr = Sr

    return std_error, r2_coef, ssr

def transform_data(xi: npt.NDArray[np.float64],
    f_basis: List[Callable]) -> npt.NDArray[np.float64]:
    """Transforms univariate data to custom basis functions
    for multivariate regression."""
    columns = [f(xi) for f in f_basis]
    return np.column_stack(columns)

def linear_regression(xi: npt.NDArray[np.float64],
    yi: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Returns the coefficients of linear regression (a0, a1) such that
    the linear function f(x) = a0+a1*x best describes the data (xi, yi)."""

    n = xi.shape[0]

    sum_x, sum_y = np.sum(xi), np.sum(yi)
    sum_xy, sum_x2 = np.dot(xi,yi), np.dot(xi,xi)

    x_mean, y_mean = sum_x/n, sum_y/n

    a = np.zeros(2)
    a[1] = (n*sum_xy - sum_x*sum_y)/( n*sum_x2 - sum_x**2 )
    a[0] = y_mean - a[1]*x_mean
    return a

def polynomial_regression(xi: npt.NDArray[np.float64],
    yi: npt.NDArray[np.float64], m: int = 1) -> npt.NDArray[np.float64]:
    """Returns the coefficients of polynomial regression (a0, a1, ..., am) such that
    the linear function f(x) = a0+a1*x+...+am*x best describes the data (xi, yi)."""

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

    a_coef = lu.solve(A,b)

    return a_coef

def multi_linear_regression(xi: npt.NDArray[np.float64],
    yi: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Returns the coefficients of multi-linear regression (a0, a1, ..., am)
    such that the multi-linear function f(x1, x2, ..., xm) = a0+a1*x1+...+am*xm
    best describes the data (xi[n x m], yi[n])."""

    lu = direct_solver.LUSolver()

    n, ndim = xi.shape[0], xi.shape[1]
    if (n < ndim+1):
        print('Not enough data points for multiple linear regression!')
        return None

    X = np.column_stack((np.ones(n), xi))

    A = np.dot(X.T, X)

    b = np.dot(X.T, yi)

    lu = direct_solver.LUSolver()
    a_coef = lu.solve(A, b)

    return a_coef

def non_linear_regression(
    f: Callable[ [npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]],
                npt.NDArray[np.float64] ],
    xi: npt.NDArray[np.float64], yi: npt.NDArray[np.float64],
    u0: npt.NDArray[np.float64], output: bool = False, tol: float = 1e-8,
    k_max: int = 1000, r: float = 1.0) -> npt.NDArray[np.float64]:
    """Returns the coefficients of non-linear regression (u0, u1, ..., um)
    such that the non linear function f(x1, x2, ..., xm)
    best describes the data (xi[n x m], yi[n])."""

    ls_solver = direct_solver.LUSolver()

    nr_solver = newton_solver.LevenbergMarquardtSolver(ls_solver=ls_solver, u0=u0,
                            k_max=k_max, tol=tol, r=r, l0=1e-4, scale=2.0)

    problem = non_linear_problem.Regression(f=f, xi=xi, yi=yi)

    u = nr_solver.solve(problem, output=output)

    return u
