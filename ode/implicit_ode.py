
from typing import Callable

import numpy as np

from linear_systems import direct_solver
from non_linear_systems import newton_solver, non_linear_problem
from ode import butcher_tableaus, explicit_ode

def implicit_solver(f: Callable[[float, np.ndarray[tuple[int]]], np.ndarray[tuple[int]]],
    t0: float, tf: float, y0: np.ndarray[tuple[int]], h: float, method: str = 'euler'
    ) -> tuple[np.ndarray[tuple[int]], np.ndarray[tuple[int, int]]]:
    """Returns the solution of a system of ODEs using a single-stage implicit method.

    Args:
        f: the right-hand side function of the ODEs
        t0: the starting value of the independent variable
        tf: the final value of the independent variable
        y0: the initial conditions
        h: the step of integration
        method (optional): the method of integration
    Returns:
        A tuple of arrays of the independent variable and the solutions of the ODEs."""

    if method not in butcher_tableaus.IMPLICIT_METHODS:
        raise ValueError(f"Method {method} is not supported.")

    method_dict = butcher_tableaus.IMPLICIT_METHODS[method]

    n = int(round((tf-t0)/h))
    neq = y0.shape[0]

    t = np.array([t0+i*h for i in range(n+1)])

    y = np.zeros((n+1, neq))

    y[0,:] = y0.copy()

    ns = method_dict['ns']

    ls_solver = direct_solver.LUSolver()

    u0 = np.zeros((ns*neq))
    for j in range(ns):
        u0[j*neq:(j+1)*neq] = y0.copy()

    nr_solver = newton_solver.NewtonSolver(ls_solver, u0=u0, k_max=100, tol=1e-8, r=1.0)

    problem = non_linear_problem.ImplicitODE(f=f, yi=y0, ti=t[0], h=h, method_dict=method_dict)

    for i in range(n):

        ti, yi = t[i], y[i]

        y[i+1] = solve_step(f, ti, yi, h, method_dict, nr_solver, problem)

    return t, y

def solve_step(f: Callable[[float, np.ndarray[tuple[int]]], np.ndarray[tuple[int]]],
    ti: float, yi: np.ndarray[tuple[int]], h: float, method_dict: dict,
    nr_solver: newton_solver.NewtonSolver, problem: non_linear_problem.NonlinearProblem
    ) -> np.ndarray[tuple[int]]:
    """Returns the solution to a single step of the implicit ODE solver.

    Args:
        f: the right-hand side function of the ODEs
        ti: the value of the independent variable at the previous step
        yi: the values of the independent variables at the previous step
        h: the step of integration
        method_dict: the dictionary that contains the method's info
        nr_solver: the Newton solver object
        problem: the non-linear problem object
    Returns:
        The solution array of the ODEs at the current step."""

    ns = method_dict['ns']
    b = method_dict['b']
    p = method_dict['p']

    neq = yi.shape[0]

    problem.update(yi=yi, ti=ti, h=h)

    u0 = np.zeros((ns*neq))
    for j in range(ns):
        u0[j*neq:(j+1)*neq] = yi.copy()
    nr_solver.update_guess(u0)

    u = nr_solver.solve(problem, output=False)

    yn = yi.copy()
    for j in range(ns):

        Y_j = u[j*neq:(j+1)*neq]
        t_j = ti + p[j]*h
        k_j = f(t_j, Y_j)

        yn += h*b[j]*k_j

    return yn

def implicit_adaptive_solver(f: Callable[[float, np.ndarray[tuple[int]]], np.ndarray[tuple[int]]],
    t0: float, tf: float, y0: np.ndarray[tuple[int]], h_min: float = 1e-3, h_max: float = 1.0,
    tol: float = 1e-4, method: str = 'euler') -> tuple[np.ndarray[tuple[int]], np.ndarray[tuple[int, int]]]:
    """Returns the solution of a system of ODEs using a single-stage implicit method
    with adaptive step.

    Args:
        f: the right-hand side function of the ODEs
        t0: the starting value of the independent variable
        tf: the final value of the independent variable
        y0: the initial conditions
        h_min (optional): the minimum step of integration
        h_max (optional): the maximum step of integration
        tol (optional): the tolerance for adapting the step size
        method (optional): the method of integration
    Returns:
        A tuple of arrays of the independent variable and the solutions of the ODEs."""

    if method not in butcher_tableaus.IMPLICIT_METHODS:
        raise ValueError(f"Method {method} is not supported.")

    method_dict = butcher_tableaus.IMPLICIT_METHODS[method]

    ns = method_dict['ns']
    order = method_dict['order']

    neq = y0.shape[0]

    t = [t0]
    y = [y0.copy()]

    h = h_max

    ti = t0
    yi = y0.copy()

    ls_solver = direct_solver.LUSolver()

    u0 = np.zeros((ns*neq))
    for j in range(ns):
        u0[j*neq:(j+1)*neq] = y0.copy()

    nr_solver = newton_solver.NewtonSolver(ls_solver, u0=u0, k_max=100, tol=1e-8, r=1.0)

    problem = non_linear_problem.ImplicitODE(f=f, yi=y0, ti=t[0], h=h, method_dict=method_dict)

    while ti < tf:

        if ti+h > tf:
            h = tf - ti

        # solve with step h
        y_low = solve_step(f, ti, yi, h, method_dict, nr_solver, problem)

        # solve with step h/2 twice
        y_mid = solve_step(f, ti, yi, h/2, method_dict, nr_solver, problem)
        y_high = solve_step(f, ti+h/2, y_mid, h/2, method_dict, nr_solver, problem)

        # new step
        error_diff = np.max(np.abs(y_high - y_low))
        max_err = error_diff/((2.0**order) - 1.0)

        exponent = 1.0/(order + 1.0)
        h_new = 0.9*h*((tol/(max_err + 1e-20))**exponent)

        if max_err <= tol or h <= h_min:
            ti = ti + h
            yi = y_high
            t.append(ti)
            y.append(yi)

        h = min( max(h_min, h_new), h_max )

    return np.array(t), np.array(y)

def adams_moulton(f: Callable[[float, np.ndarray[tuple[int]]], np.ndarray[tuple[int]]],
    t0: float, tf: float, y0: np.ndarray[tuple[int]], h: float,
    order: int = 4) -> tuple[np.ndarray[tuple[int]], np.ndarray[tuple[int, int]]]:
    """Returns the solution of a system of ODEs using the Adams-Moulton method.

    Args:
        f: the right-hand side function of the ODEs
        t0: the starting value of the independent variable
        tf: the final value of the independent variable
        y0: the initial conditions
        h: the step of integration
        order (optional): the order of the Adams-Moulton method
    Returns:
        A tuple of arrays of the independent variable and the solutions of the ODEs."""

    AM_COEFFS = {
        2: [0.5, 0.5],
        3: [5.0/12.0, 8.0/12.0, -1.0/12.0],
        4: [9.0/24.0, 19.0/24.0, -5.0/24.0, 1.0/24.0],
        5: [251.0/720.0, 646.0/720.0, -264.0/720.0,
            106.0/720.0, -19.0/720.0]
    }

    if order not in AM_COEFFS:
        raise ValueError(f"Order {order} is not supported.")

    b = AM_COEFFS[order]

    n = int(round((tf-t0)/h))
    if n < order:
        raise ValueError(f"Order {order} larger than n.")

    neq = y0.shape[0]

    t = np.array([t0+i*h for i in range(n+1)])

    y = np.zeros((n+1, neq))

    y[0,:] = y0.copy()

    for i in range(order-1):
        ti = t[i]
        _, y_rk = explicit_ode.rk4(f, ti, ti+h, y[i], h)
        y[i+1] = y_rk[-1]

    f_past = np.array([f(t[i], y[i]) for i in range(order)])

    ls_solver = direct_solver.LUSolver()

    u0 = y0.copy()

    nr_solver = newton_solver.NewtonSolver(ls_solver, u0=u0, k_max=100, tol=1e-8, r=1.0)

    problem = non_linear_problem.AdamsMoultonODE(f=f, yi=y0, ti=t[0], h=h, b0=b[0])

    for i in range(order-1, n):

        y_next = y[i,:]

        # for k in range(order-1):
        for k in range(1, order):
            # f_j = f_past[-(k+1)]
            f_j = f_past[-k]
            for ieq in range(neq):
                y_next[ieq] += h*b[k]*f_j[ieq]

        problem.update(yi=y_next, ti=ti)

        y[i+1] = nr_solver.solve(problem, output=False)

        if i+1 < n:
            f_past = f_past[1:] + [f(t[i+1], y[i+1])]

    return t, y
