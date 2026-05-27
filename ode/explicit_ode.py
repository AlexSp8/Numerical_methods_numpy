
from typing import Callable

import numpy as np

from linear_systems import direct_solver

def heun(f: Callable[[float, np.ndarray[tuple[int]]], np.ndarray[tuple[int]]],
    t0: float, tf: float, y0: np.ndarray[tuple[int]], h: float, k_max: int = 1
    ) -> tuple[np.ndarray[tuple[int]], np.ndarray[tuple[int, int]]]:
    """Returns the solution of a system of ODEs using the Heun's method.

    Args:
        f: the right-hand side function of the ODEs
        t0: the starting value of the independent variable
        tf: the final value of the independent variable
        y0: the initial conditions
        h: the step of integration
        k_max (optional): the inner correction maximum iterations
    Returns:
        A tuple of arrays of the independent variable and the solutions of the ODEs."""

    n = int(round((tf-t0)/h))
    neq = y0.shape[0]

    t = np.array([t0+i*h for i in range(n+1)])

    y = np.zeros((n+1, neq))

    y[0,:] = y0.copy()

    for i in range(n):

        fi = f(t[i], y[i])

        y_next = y[i,:] + h*fi[:]

        for k in range(k_max):

            y_old = y_next[:]

            fi_next =  f(t[i]+h, y_next)

            phi = (fi[:] + fi_next[:])/2.0
            y_next[:] = y[i,:] + h*phi[:]

            err = abs( (y_next[:]-y_old[:])/(y_next[:]+1e-8) )

            if max(err) < 1e-8:
                # print(k)
                break

        y[i+1] = y_next[:]
    return t, y

def rk2(f: Callable[[float, np.ndarray[tuple[int]]], np.ndarray[tuple[int]]],
    t0: float, tf: float, y0: np.ndarray[tuple[int]], h: float, a2: float = 0.5
    ) -> tuple[np.ndarray[tuple[int]], np.ndarray[tuple[int, int]]]:
    """Returns the solution of a system of ODEs using the RK2 method with specified a2.

    Args:
        f: the right-hand side function of the ODEs
        t0: the starting value of the independent variable
        tf: the final value of the independent variable
        y0: the initial conditions
        h: the step of integration
        a2 (optional): the value to specify the RK2 method
    Returns:
        A tuple of arrays of the independent variable and the solutions of the ODEs."""

    a1 = 1.0-a2
    p1 = 0.5/(a2+1e-12)
    q11 = p1

    n = int(round((tf-t0)/h))
    neq = y0.shape[0]

    t = np.array([t0+i*h for i in range(n+1)])

    y = np.zeros((n+1, neq))

    y[0,:] = y0.copy()

    for i in range(n):

        ti, yi = t[i], y[i]

        k1 = f(ti,yi)

        y2 = yi[:] + q11*h*k1[:]

        k2 = f(ti+p1*h, y2)

        phi = a1*k1[:] + a2*k2[:]
        y[i+1,:] = yi[:] + h*phi[:]

    return t, y

def rk3_coefficients(p1: float = 0.5, p2: float = 1.0) -> np.ndarray[tuple[int]]:
    """Returns the coefficients for the RK3 method with specified p1, p2.

    Args:
        p1, p2 (optional): the values to specify the RK3 method
    Returns:
        The remaining RK3 coefficients."""

    if p1 < 0 or p2 < 0:
        raise ValueError("For RK3 choose p1 > 0 and p2 > 0.")

    if p1 == p2:
        raise ValueError("For RK3 choose p1 != p2.")

    if abs(p1 - 2/3.0) < 1e-8:
        raise ValueError("For RK3 choose p1 != 2/3 (a3 = 0 leads to division by zero for q22).")

    A = np.array([ [p1, p2], [p1**2, p2**2] ])
    b = np.array([0.5, 1.0/3.0])

    lu = direct_solver.LUSolver()
    x = lu.solve(A,b)

    a2, a3 = x[0], x[1]

    a1 = 1.0 - a2 - a3
    q22 = 1.0/(6.0*p1*a3)

    q11 = p1
    q21 = p2 - q22

    # A = np.array([ [a2, a3], [a2*p1, a3*p2] ])
    # b = np.array([0.5 - a3*q22, (1.0/3.0) - a3*p2*q22])

    # x = lu.solve(A,b)
    # q11, q21 = x[0], x[1]

    return np.array([a1, a2, a3, q11, q21, q22])

