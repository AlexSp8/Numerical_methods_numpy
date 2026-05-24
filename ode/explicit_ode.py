
from typing import Callable

import numpy as np
import numpy.typing as npt

from linear_systems import direct_solver

def explicit_euler(f: Callable[[float, npt.NDArray], npt.NDArray], x0: float,
    xf: float, h: float, y0: npt.NDArray) -> tuple[npt.NDArray, npt.NDArray]:
    """Returns the solution of a system of ODEs using the explicit Euler method."""

    n = int(round((xf-x0)/h))
    m = len(y0)

    x = np.array([x0+i*h for i in range(n+1)])
    # x = np.linspace(x0, xf, n+1)

    y = np.zeros((n+1, m))

    y[0,:] = y0[:]

    for i in range(n):
        dy = f(x[i], y[i])
        y[i+1,:] = y[i,:] + h*dy[:]

    return x, y

def heun_method(f: Callable[[float, npt.NDArray], npt.NDArray], x0: float, xf: float,
    h: float, y0: npt.NDArray, k_max: int = 1) -> tuple[npt.NDArray, npt.NDArray]:
    """Returns the solution of a system of ODEs using Heun's method."""

    n = int(round((xf-x0)/h))
    m = len(y0)

    x = np.array([x0+i*h for i in range(n+1)])

    y = np.zeros((n+1, m))

    y[0] = y0

    for i in range(n):

        fi = f(x[i], y[i])

        y_next = y[i,:] + h*fi[:]

        for k in range(k_max):

            y_old = y_next[:]

            fi_next =  f(x[i]+h, y_next)

            phi = (fi[:] + fi_next[:])/2.0
            y_next[:] = y[i,:] + h*phi[:]

            err = abs( (y_next[:]-y_old[:])/(y_next[:]+1e-8) )

            if max(err) < 1e-8:
                # print(k)
                break

        y[i+1] = y_next[:]
    return x, y

def midpoint_method(f: Callable[[float, npt.NDArray], npt.NDArray], x0: float,
    xf: float, h: float, y0: npt.NDArray) -> tuple[npt.NDArray, npt.NDArray]:
    """Returns the solution of a system of ODEs using the midpoint method."""

    n = int(round((xf-x0)/h))
    m = len(y0)

    x = np.array([x0+i*h for i in range(n+1)])

    y = np.zeros((n+1, m))

    y[0] = y0[:]

    for i in range(n):

        xi, yi = x[i], y[i]

        fi = f(xi, yi)

        xm = xi + (h/2)

        ym = yi[:] + fi[:]*h/2

        phi = f(xm, ym)

        y[i+1,:] = yi[:] + h*phi[:]

    return x, y

def rk2(f: Callable[[float, npt.NDArray], npt.NDArray], x0: float, xf: float,
    h: float, y0: npt.NDArray, a2: float = 0.5) -> tuple[npt.NDArray, npt.NDArray]:
    """Returns the solution of a system of ODEs using the RK2 method with specified a2."""

    a1 = 1.0-a2
    p1 = 0.5/(a2+1e-12)
    q11 = p1

    n = int(round((xf-x0)/h))
    m = len(y0)

    x = np.array([x0+i*h for i in range(n+1)])

    y = np.zeros((n+1, m))

    y[0] = y0[:]

    for i in range(n):

        xi, yi = x[i], y[i]

        k1 = f(xi,yi)

        y2 = yi[:] + q11*h*k1[:]

        k2 = f(xi+p1*h, y2)

        phi = a1*k1[:] + a2*k2[:]
        y[i+1,:] = yi[:] + h*phi[:]

    return x, y

def rk3(f: Callable[[float, npt.NDArray], npt.NDArray], x0: float, xf: float,
    h: float, y0: npt.NDArray, p1: float = 0.5, p2: float = 1.0) -> tuple[npt.NDArray, npt.NDArray]:
    """Returns the solution of a system of ODEs using the RK3 method with specified p1, p2.
    Standard RK3 is: p1 = 0.5, p2 = 1.0."""

    c = rk3_coefficients(p1, p2)
    a1, a2, a3, q11, q21, q22 = c[0], c[1], c[2], c[3], c[4], c[5]

    n = int(round((xf-x0)/h))
    m = len(y0)

    x = np.array([x0+i*h for i in range(n+1)])

    y = np.zeros((n+1, m))

    y[0] = y0[:]
    for i in range(n):

        xi, yi = x[i], y[i]

        k1 = f(xi,yi)

        y2 = yi[:] + q11*h*k1[:]
        k2 = f(xi+p1*h, y2)

        y3 = yi[:] + h*(q21*k1[:]+q22*k2[:])
        k3 = f(xi+p2*h, y3)

        phi = a1*k1[:] + a2*k2[:] + a3*k3[:]
        y[i+1,:] = yi[:] + h*phi[:]

    return x, y

def rk3_coefficients(p1: float = 0.5, p2: float = 1.0) -> npt.NDArray:
    """Returns the coefficients for the RK3 method with specified p1, p2."""

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

def rk4(f: Callable[[float, npt.NDArray], npt.NDArray], x0: float, xf: float,
    h: float, y0: npt.NDArray) -> tuple[npt.NDArray, npt.NDArray]:
    """Returns the solution of a system of ODEs using the standard RK4 method."""

    n = int(round((xf-x0)/h))
    m = len(y0)

    x = np.array([x0+i*h for i in range(n+1)])

    y = np.zeros((n+1, m))

    y[0] = y0[:]

    for i in range(n):

        xi, yi = x[i], y[i]

        k1 = f(xi,yi)

        y2 = yi[:] + 0.5*k1[:]*h
        k2 = f(xi+0.5*h, y2)

        y3 = yi[:] + 0.5*k2[:]*h
        k3 = f(xi+0.5*h, y3)

        y4 = yi[:] + k3[:]*h
        k4 = f(xi+h, y4)

        phi = (k1[:] + 2.0*k2[:] + 2.0*k3[:] + k4[:])/6.0
        y[i+1,:] = yi[:] + h*phi[:]

    return x, y

def rk5(f: Callable[[float, npt.NDArray], npt.NDArray], x0: float, xf: float,
    h: float, y0: npt.NDArray) -> tuple[npt.NDArray, npt.NDArray]:
    """Returns the solution of a system of ODEs using
    the standard RK5 (Butcher's) method."""

    n = int(round((xf-x0)/h))
    m = len(y0)

    x = np.array([x0+i*h for i in range(n+1)])

    y = np.zeros((n+1, m))

    y[0] = y0[:]

    for i in range(n):

        xi, yi = x[i], y[i]

        k1 = f(xi, yi)

        y2 = yi[:] + 0.25*k1[:]*h
        k2 = f(xi+0.25*h, y2)

        y3 = yi[:] + 0.125*(k1[:]+k2[:])*h
        k3 = f(xi+0.25*h, y3)

        y4 = yi[:] + (-0.5*k2[:]+k3[:])*h
        k4 = f(xi+0.5*h, y4)

        y5 = yi[:] + (3.0*k1[:]+9.0*k4[:])*h/16.0
        k5 = f(xi+0.75*h, y5)

        y6 = yi[:] + (-3.0*k1[:]+2.0*k2[:]+12.0*k3[:]-12.0*k4[:]+8.0*k5[:])*h/7.0
        k6 = f(xi+h, y6)

        phi = (7.0*k1[:] + 32.0*k3[:] + 12.0*k4[:] + 32.0*k5[:]+7.0*k6[:])/90.0
        y[i+1,:] = yi[:] + h*phi[:]

    return x, y

def rk45_fehlberg(f: Callable[[float, npt.NDArray], npt.NDArray], x0: float, xf: float,
    y0: npt.NDArray, h_min: float = 1e-6, h_max :float = 1.0,
    eps: float = 1e-8) -> tuple[npt.NDArray, npt.NDArray]:
    """Returns the solution of a system of ODEs using
    the adaptive RK45 Fehlberg method."""

    x = [x0]
    y = [np.array(y0)]

    h = h_max

    xi = x0
    yi = y0[:]

    while xi < xf:

        if xi+h > xf:
            h = xf - xi

        k1 = f(xi, yi)

        y2 = yi[:] + h*0.25*k1[:]
        k2 = f(xi+0.25*h, y2)

        y3 = yi[:] + h*((3.0/32.0)*k1[:] + (9.0/32.0)*k2[:])
        k3 = f(xi+0.375*h, y3)

        y4 = yi[:] + h*((1932.0/2197.0)*k1[:] - (7200.0/2197.0)*k2[:]
                        + (7296.0/2197.0)*k3[:])
        k4 = f(xi+(12.0/13.0)*h, y4)

        y5 = yi[:] + h*((439.0/216.0)*k1[:] - 8.0*k2[:] + (3680.0/513.0)*k3[:]
                        - (845.0/4104.0)*k4[:])
        k5 = f(xi+h, y5)

        y6 = yi[:] + h*((-8.0/27.0)*k1[:] + 2.0*k2[:] - (3544.0/2565.0)*k3[:]
                        + (1859.0/4104.0)*k4[:] - (11.0/40.0)*k5[:])
        k6 = f(xi+0.5*h, y6)

        y_h4 = yi[:] + h*((25.0/216.0)*k1[:] + (1408.0/2565.0)*k3[:]
                            + (2197.0/4104.0)*k4[:] - 0.2*k5[:])

        y_h5 = yi[:] + h*((16.0/135.0)*k1[:] + (6656.0/12825.0)*k3[:]
                            + (28561.0/56430.0)*k4[:] - (9.0/50.0)*k5[:]
                            + (2.0/55.0)*k6[:])

        max_err = np.max(np.abs(y_h5[:]-y_h4[:]))

        h_new = 0.9*h*((eps/(max_err+1e-20))**0.2)

        if max_err < eps or h < h_min:
            xi = xi + h
            yi = y_h5
            x.append(xi)
            y.append(yi)

        h = min( max(h_min, h_new), h_max )

    return np.array(x), np.array(y)

def rk45_cash_karp(f: Callable[[float, npt.NDArray], npt.NDArray], x0: float, xf: float,
    y0: npt.NDArray, h_min: float = 1e-6, h_max :float = 1.0,
    eps: float = 1e-8) -> tuple[npt.NDArray, npt.NDArray]:
    """Returns the solution of a system of ODEs using
    the adaptive RK45 Cash-Karp method."""

    x = [x0]
    y = [np.array(y0)]

    h = h_max

    xi = x0
    yi = y0[:]


    while xi < xf:

        if xi+h > xf:
            h = xf - xi

        k1 = f(xi, yi)

        y2 = yi[:] + h*0.2*k1[:]
        k2 = f(xi+0.2*h, y2)

        y3 = yi[:] + h*((3.0/40.0)*k1[:] + (9.0/40.0)*k2[:])
        k3 = f(xi+0.3*h, y3)

        y4 = yi[:] + h*(0.3*k1[:] - 0.9*k2[:] + 1.2*k3[:])
        k4 = f(xi+0.6*h, y4)

        y5 = yi[:] + h*(-(11.0/54.0)*k1[:] + 2.5*k2[:] - (70.0/27.0)*k3[:]
                        + (35.0/27.0)*k4[:])
        k5 = f(xi+h, y5)

        y6 = yi[:] + h*((1631.0/55296.0)*k1[:] + (175.0/512.0)*k2[:] + (575.0/13824.0)*k3[:]
                        + (44275.0/110592.0)*k4[:] + (253.0/4096.0)*k5[:])
        k6 = f(xi+(7.0*h/8.0), y6)

        y_h4 = yi[:] + h*((2825.0/27648.0)*k1[:] + (18575.0/48384.0)*k3[:]
                            + (13525.0/55296.0)*k4[:] + (277.0/14336.0)*k5[:]
                            + 0.25*k6[:])

        y_h5 = yi[:] + h*((37.0/378.0)*k1[:] + (250.0/621.0)*k3[:]
                            + (125.0/594.0)*k4[:] + (512.0/1771.0)*k6[:])

        max_err = max(abs(y_h5[:]-y_h4[:]) )

        h_new = 0.9*h*((eps/(max_err+1e-20))**0.2)

        if max_err < eps or h < h_min:
            xi = xi + h
            yi = y_h5
            x.append(xi)
            y.append(yi)

        h = min( max(h_min,h_new), h_max )

    return np.array(x), np.array(y)

def rk45_dormand_prince(f: Callable[[float, npt.NDArray], npt.NDArray], x0: float, xf: float,
    y0: npt.NDArray, h_min: float = 1e-6, h_max :float = 1.0,
    eps: float = 1e-8) -> tuple[npt.NDArray, npt.NDArray]:
    """Returns the solution of a system of ODEs using
    the adaptive RK45 Dormand-Prince method."""

    x = [x0]
    y = [np.array(y0)]

    h = h_max

    xi = x0
    yi = y0[:]

    k1 = f(x0, y0)

    while xi < xf:

        if xi+h > xf:
            h = xf - xi

        y2 = yi[:] + h*0.2*k1[:]
        k2 = f(xi+0.2*h, y2)

        y3 = yi[:] + h*((3.0/40.0)*k1[:] + (9.0/40.0)*k2[:])
        k3 = f(xi+0.3*h, y3)

        y4 = yi[:] + h*((44.0/45.0)*k1[:] - (56.0/15.0)*k2[:] + (32.0/9.0)*k3[:])
        k4 = f(xi+0.8*h, y4)

        y5 = yi[:] + h*((19372.0/6561.0)*k1[:] - (25360.0/2187.0)*k2[:] + (64448.0/6561.0)*k3[:]
                        - (212.0/729.0)*k4[:])
        k5 = f(xi+(8.0/9.0)*h, y5)

        y6 = yi[:] + h*((9017.0/3168.0)*k1[:] - (355.0/33.0)*k2[:] + (46732.0/5247.0)*k3[:]
                        + (49.0/176.0)*k4[:] - (5103.0/18656.0)*k5[:])
        k6 = f(xi+h, y6)

        y_h5 = yi[:] + h*((35.0/384.0)*k1[:] + (500.0/1113.0)*k3[:] + (125.0/192.0)*k4[:]
                        - (2187.0/6784.0)*k5[:] + (11.0/84.0)*k6[:])
        k7 = f(xi+h, y_h5)

        y_h4 = yi[:] + h*((5179.0/57600.0)*k1[:] + (7571.0/16695.0)*k3[:]
                            + (393.0/640.0)*k4[:] - (92097.0/339200.0)*k5[:]
                            + (187.0/2100.0)*k6[:] + (1.0/40.0)*k7[:])

        max_err = max( abs(y_h5[:]-y_h4[:]) )

        h_new = 0.9*h*((eps/(max_err+1e-20))**0.2)

        if max_err < eps or h < h_min:
            xi = xi + h
            yi = y_h5
            x.append(xi)
            y.append(yi)
            k1 = k7

        h = min( max(h_min,h_new), h_max )

    return np.array(x), np.array(y)

def adams_bashforth(f: Callable[[float, npt.NDArray], npt.NDArray],
    x0: float, xf: float, h: float, y0: npt.NDArray,
    order: int = 1) -> tuple[npt.NDArray, npt.NDArray]:
    """Returns the solution of a system of ODEs using an
    Adams-Bashforth explicit method of specified order."""

    AB_COEFFS = {
        1: [1.0],
        2: [1.5, -0.5],
        3: [23.0/12.0, -16.0/12.0, 5.0/12.0],
        4: [55.0/24.0, -59.0/24.0, 37.0/24.0, -9.0/24.0],
        5: [1901.0/720.0, -2774.0/720.0, 2616.0/720.0,
            -1274.0/720.0, 251.0/720.0]
    }

    if order not in AB_COEFFS:
        raise ValueError(f"Order {order} is not supported.")

    b = AB_COEFFS[order]

    n = int(round((xf-x0)/h))
    if n < order:
        raise ValueError(f"Order {order} larger than n.")
        # return explicit_ode.rk4(f, x0, xf, h, y0)

    neq = len(y0)

    x = np.array([x0+i*h for i in range(n+1)])

    y = np.zeros((n+1, neq))

    y[0] = y0[:]

    for i in range(order-1):
        xi = x[i]
        _, y_rk = rk4(f, xi, xi+h, h, y[i])
        y[i+1] = y_rk[-1]

    f_past = [f(x[i], y[i]) for i in range(order)]

    for i in range(order-1, n):
        y_sum = np.zeros(neq)
        for k in range(order):
            y_sum += b[k]*f_past[-(k+1)]

        y[i+1] = y[i] + h*y_sum
        if i+1 < n:
            f_past = f_past[1:] + [f(x[i+1], y[i+1])]

    return x, y

def abm_pece(f: Callable[[float, npt.NDArray], npt.NDArray],
    x0: float, xf: float, h: float, y0: npt.NDArray,
    order: int = 1) -> tuple[npt.NDArray, npt.NDArray]:
    """Returns the solution of a system of ODEs using an
    Adams-Bashforth-Moulton PECE (Predict-Evaluate-Correct-Evaluate)
    explicit method of specified order."""

    AB_COEFFS = {
        1: [1.0],
        2: [1.5, -0.5],
        3: [23/12, -16/12, 5/12],
        4: [55/24, -59/24, 37/24, -9/24]
    }
    if order not in AB_COEFFS:
        raise ValueError(f"Order {order} is not supported.")

    AM_COEFFS = {
        1: [0.5, 0.5],
        2: [5/12, 8/12, -1/12],
        3: [9/24, 19/24, -5/24, 1/24],
        4: [251.0/720.0, 646.0/720.0, -264.0/720.0,
            106.0/720.0, -19.0/720.0]
    }
    if order not in AM_COEFFS:
        raise ValueError(f"Order {order} is not supported.")

    # We match the AB order to the AM order to keep accuracy consistent
    # AB of order k predicts for AM of order k+1
    b_ab = AB_COEFFS[order]
    b_am = AM_COEFFS[order]

    n = int(round((xf-x0)/h))
    if n < order:
        raise ValueError(f"Order {order} larger than n.")

    neq = len(y0)

    x = np.array([x0+i*h for i in range(n+1)])

    y = np.zeros((n+1, neq))

    y[0] = y0[:]

    for i in range(order-1):
        xi = x[i]
        _, y_rk = rk4(f, xi, xi+h, h, y[i])
        y[i+1] = y_rk[-1]

    f_past = [f(x[i], y[i]) for i in range(order)]

    for i in range(order-1, n):

        y_ab = np.zeros(neq)
        for k in range(order):
            y_ab += b_ab[k]*f_past[-(k+1)]

        y_predict = y[i] + h*y_ab

        f_predict = f(x[i+1], y_predict)

        y_am = b_am[0]*f_predict
        for k in range(1, order+1):
            y_am += b_am[k]*f_past[-k]

        y[i+1] = y[i] + h*y_am

        if i + 1 < n:
            f_past = f_past[1:] + [f(x[i+1], y[i+1])]

    return x, y
