
from typing import Callable
from abc import ABC, abstractmethod
import numpy as np

from linear_systems import direct_solver
from non_linear_systems import newton_solver, non_linear_problem
from ode import butcher_tableaus


class IVPSolver(ABC):

    def __init__(self, f: Callable[[float, np.ndarray[tuple[int]]], np.ndarray[tuple[int]]],
        t0: float, tf: float, y0: np.ndarray[tuple[int]], method: str):

        self.f = f
        self.t0 = t0
        self.tf = tf
        self.y0 = y0
        self.method = method

    @abstractmethod
    def solve(self) -> tuple[np.ndarray[tuple[int]], np.ndarray[tuple[int, int]]]:
        """Returns the solution of a system of ODEs.

        Args:
            self: the IVPSolver object
        Returns:
            A tuple of arrays of the independent variable and the solutions of the ODEs."""

        pass


class ExplicitIVPSolver(IVPSolver):

    def __init__(self, *args, h: float, **kwargs):

        super().__init__(*args, **kwargs)

        method = self.method
        if method not in butcher_tableaus.EXPLICIT_METHODS:
            raise ValueError(f"Method {method} is not supported.")

        self.method_dict = butcher_tableaus.EXPLICIT_METHODS[method]

        self.h = h

    def solve(self) -> tuple[np.ndarray[tuple[int]], np.ndarray[tuple[int, int]]]:

        t0, tf, h = self.t0, self.tf, self.h
        y0 = self.y0
        f = self.f
        method_dict = self.method_dict

        ns = method_dict['ns']
        b = method_dict['b']
        p = method_dict['p']
        q = method_dict['q']

        n = int(round((tf-t0)/h))
        neq = y0.shape[0]

        t = np.array([t0+i*h for i in range(n+1)])
        # t = np.linspace(t0, tf, n+1)

        y = np.zeros((n+1, neq))

        y[0,:] = y0.copy()

        for i in range(n):

            ti, yi = t[i], y[i,:]

            k_mat = np.zeros((ns, neq))

            for j in range(ns):

                yj = yi[:] + h*np.dot(q[j, :j], k_mat[:j, :])
                k_mat[j, :] = f(ti + p[j]*h, yj)

            y[i+1,:] = yi[:] + h*np.dot(b, k_mat)

        return t, y


class ExplicitAdaptiveIVPSolver(IVPSolver):

    def __init__(self, *args, h_min: float = 1e-6, h_max: float = 1.0,
        atol: float = 1e-6, rtol: float = 1e-3, **kwargs):

        super().__init__(*args, **kwargs)

        method = self.method
        if method not in butcher_tableaus.ADAPTIVE_EXPLICIT_METHODS:
            raise ValueError(f"Method {method} is not supported.")

        self.method_dict = butcher_tableaus.ADAPTIVE_EXPLICIT_METHODS[method]

        self.h_min = h_min
        self.h_max = h_max
        self.atol = atol
        self.rtol = rtol

    def solve(self) -> tuple[np.ndarray[tuple[int]], np.ndarray[tuple[int, int]]]:

        t0, tf = self.t0, self.tf
        y0 = self.y0
        f = self.f
        h_min, h_max, atol, rtol = self.h_min, self.h_max, self.atol, self.rtol
        method_dict = self.method_dict

        ns = method_dict['ns']
        b_high = method_dict['b_high']
        b_low = method_dict['b_low']
        p = method_dict['p']
        q = method_dict['q']
        order = method_dict['order_high']

        exponent = 1.0/order

        neq = y0.shape[0]

        t = [t0]
        y = [np.array(y0)]

        h = h_max

        ti = t0
        yi = y0.copy()

        while ti < tf:

            if ti+h > tf:
                h = tf - ti

            k_mat = np.zeros((ns, neq))

            for j in range(ns):

                yj = yi[:] + h*np.dot(q[j, :j], k_mat[:j, :])
                k_mat[j, :] = f(ti + p[j]*h, yj)

            y_high = yi[:] + h*np.dot(b_high, k_mat)
            y_low = yi[:] + h*np.dot(b_low, k_mat)

            # new step
            error_diff = np.abs(y_high[:]-y_low[:])
            scale = atol + rtol*np.maximum(np.abs(yi), np.abs(y_high))
            max_err = np.max(error_diff/scale)
            # max_err /= (2.0**order)-1.0 # for step-doubling (Richardson extrapolation)

            h_new = 0.9*h*((1.0/(max_err + 1e-20))**exponent)

            if max_err <= 1.0 or h <= h_min:
                ti = ti + h
                yi = y_high
                t.append(ti)
                y.append(yi)

            h = min( max(h_min, h_new), h_max )

        return np.array(t), np.array(y)


class ImplicitIVPSolver(IVPSolver):

    def __init__(self, *args, h: float, **kwargs):

        super().__init__(*args, **kwargs)

        method = self.method
        if method not in butcher_tableaus.IMPLICIT_METHODS:
            raise ValueError(f"Method {method} is not supported.")

        self.method_dict = butcher_tableaus.IMPLICIT_METHODS[method]

        self.h = h

    def solve(self) -> tuple[np.ndarray[tuple[int]], np.ndarray[tuple[int, int]]]:

        t0, tf, h = self.t0, self.tf, self.h
        y0 = self.y0
        f = self.f
        method_dict = self.method_dict

        ns = method_dict['ns']
        b = method_dict['b']

        n = int(round((tf-t0)/h))
        neq = y0.shape[0]

        t = np.array([t0+i*h for i in range(n+1)])

        y = np.zeros((n+1, neq))

        y[0,:] = y0.copy()

        ls_solver = direct_solver.LUSolver()

        u0 = np.zeros((ns*neq))
        for j in range(ns):
            u0[j*neq:(j+1)*neq] = f(t0, y0)

        nr_solver = newton_solver.NewtonSolver(ls_solver, u0=u0, k_max=100, tol=1e-8, r=1.0)

        problem = non_linear_problem.ImplicitODE(f=f, yi=y0, ti=t0, h=h, method_dict=method_dict)

        k_mat = np.zeros((ns, neq))
        for i in range(n):

            ti, yi = t[i], y[i,:]

            problem.update(yi=yi, ti=ti)

            u = nr_solver.solve(problem, output=False)
            nr_solver.update_guess(u)

            # k_mat = u.reshape(ns, neq)
            for j in range(ns):
                k_mat[j,:] = u[j*neq:(j+1)*neq]

            y[i+1,:] = yi[:] + h*np.dot(b, k_mat)

        return t, y


class ImplicitAdaptiveIVPSolver(IVPSolver):

    def __init__(self, *args, h_min: float = 1e-3, h_max: float = 1.0,
        atol: float = 1e-6, rtol: float = 1e-3, **kwargs):

        super().__init__(*args, **kwargs)

        method = self.method
        if method not in butcher_tableaus.ADAPTIVE_IMPLICIT_METHODS:
            raise ValueError(f"Method {method} is not supported.")

        self.method_dict = butcher_tableaus.ADAPTIVE_IMPLICIT_METHODS[method]

        self.h_min = h_min
        self.h_max = h_max
        self.atol = atol
        self.rtol = rtol

    def solve(self) -> tuple[np.ndarray[tuple[int]], np.ndarray[tuple[int, int]]]:

        t0, tf = self.t0, self.tf
        y0 = self.y0
        f = self.f
        h_min, h_max, atol, rtol = self.h_min, self.h_max, self.atol, self.rtol
        method_dict = self.method_dict

        ns = method_dict['ns']
        b_high = method_dict['b_high']
        if 'b_low' in method_dict:
            b_low = method_dict['b_low']
        if 'd_weights' in method_dict:
            d_weights = method_dict['d_weights']
        order = method_dict['order_high']

        exponent = 1.0/order

        neq = y0.shape[0]

        t = [t0]
        y = [np.array(y0)]

        h = h_max

        ti = t0
        yi = y0.copy()

        ls_solver = direct_solver.LUSolver()

        u0 = np.zeros((ns*neq))
        for j in range(ns):
            u0[j*neq:(j+1)*neq] = f(t0, y0)

        nr_solver = newton_solver.NewtonSolver(ls_solver, u0=u0, k_max=100, tol=1e-8, r=1.0)

        problem = non_linear_problem.ImplicitODE(f=f, yi=y0, ti=t0, h=h, method_dict=method_dict)

        k_mat = np.zeros((ns, neq))

        while ti < tf:

            if ti+h > tf:
                h = tf - ti

            problem.update(yi=yi, ti=ti, h=h)

            u = nr_solver.solve(problem, output=False)

            # k_mat = u.reshape(ns, neq)
            for j in range(ns):
                k_mat[j,:] = u[j*neq:(j+1)*neq]

            y_high = yi[:] + h*np.dot(b_high, k_mat)

            if 'b_low' in method_dict:
                y_low = yi[:] + h*np.dot(b_low, k_mat)
                # new step
                error_diff = np.abs(y_high[:]-y_low[:])
                scale = atol + rtol*np.maximum(np.abs(yi), np.abs(y_high))
                max_err = np.max(error_diff/scale)
                # max_err /= (2.0**order)-1.0 # for step-doubling (Richardson extrapolation)
            elif 'd_weights' in method_dict:
                max_err = np.linalg.norm(np.dot(d_weights, k_mat))
            else:
                raise ValueError("Unsupported Implicit Adaptive Method!")

            h_new = 0.9*h*((1.0/(max_err + 1e-20))**exponent)

            if max_err <= 1.0 or h <= h_min:

                ti = ti + h
                yi = y_high
                t.append(ti)
                y.append(yi)

                nr_solver.update_guess(u)

            h = min( max(h_min, h_new), h_max )

        return np.array(t), np.array(y)


class ExplicitMultistepIVPSolver(IVPSolver):

    def __init__(self, *args, h: float, order: int = 1, **kwargs):

        super().__init__(*args, **kwargs)

        method = butcher_tableaus.MULTISTEP_METHODS['adams-bashforth']

        if order not in method:
            raise ValueError(f"Order {order} is not supported.")

        self.order = order
        self.ns = method[order]['ns']
        self.b = method[order]['b']
        self.h = h

    def solve(self) -> tuple[np.ndarray[tuple[int]], np.ndarray[tuple[int, int]]]:

        t0, tf, h = self.t0, self.tf, self.h
        y0 = self.y0
        f = self.f
        order, ns, b = self.order, self.ns, self.b

        n = int(round((tf-t0)/h))
        if n < order:
            raise ValueError(f"Order {order} larger than n.")
            # rk4 = ExplicitIVPSolver(f, t0, tf, y0, 'rk4', h)
            # return rk4.solve()

        neq = y0.shape[0]

        t = np.array([t0+i*h for i in range(n+1)])

        y = np.zeros((n+1, neq))

        y[0,:] = y0.copy()

        if order > 1:
            rk4 = ExplicitIVPSolver(f, t0, t[order-1], y0, 'rk4', h=h)
            _, y[:order,:] = rk4.solve()

        f_past = [f(t[i], y[i,:]) for i in range(order)]

        for i in range(order-1, n):

            y_sum = np.zeros(neq)
            for k in range(order):
                y_sum += b[k]*f_past[-(k+1)]

            y[i+1,:] = y[i,:] + h*y_sum

            if i+1 < n:
                f_past = f_past[1:] + [f(t[i+1], y[i+1])]

        return t, y


class ImplicitMultistepIVPSolver(IVPSolver):

    def __init__(self, *args, h: float, order: int = 1, **kwargs):

        super().__init__(*args, **kwargs)

        method = butcher_tableaus.MULTISTEP_METHODS['adams-moulton']

        if order not in method:
            raise ValueError(f"Order {order} is not supported.")

        self.order = order
        self.ns = method[order]['ns']
        self.b = method[order]['b']
        self.h = h

    def solve(self) -> tuple[np.ndarray[tuple[int]], np.ndarray[tuple[int, int]]]:

        t0, tf, h = self.t0, self.tf, self.h
        y0 = self.y0
        f = self.f
        order, ns, b = self.order, self.ns, self.b

        n = int(round((tf-t0)/h))
        if n < order:
            raise ValueError(f"Order {order} larger than n.")
            # rk4 = ExplicitIVPSolver(f, t0, tf, y0, 'rk4', h)
            # return rk4.solve()

        neq = y0.shape[0]

        t = np.array([t0+i*h for i in range(n+1)])

        y = np.zeros((n+1, neq))

        y[0,:] = y0.copy()

        if order > 1:
            rk4 = ExplicitIVPSolver(f, t0, t[order-1], y0, 'rk4', h=h)
            _, y[:order,:] = rk4.solve()

        f_past = [f(t[i], y[i,:]) for i in range(order)]

        ls_solver = direct_solver.LUSolver()

        u0 = y0.copy()

        nr_solver = newton_solver.NewtonSolver(ls_solver, u0=u0, k_max=100, tol=1e-8, r=1.0)

        problem = non_linear_problem.AdamsMoultonODE(f=f, yi=y0, ti=t0, h=h, b0=b[0])

        for i in range(order-1, n):

            yn = y[i,:].copy()

            for k in range(1, order):
                yn[:] += h*b[k]*f_past[-k]

            problem.update(yi=yn, ti=t[i+1])

            nr_solver.update_guess(yn)

            y[i+1,:] = nr_solver.solve(problem, output=False)

            if i+1 < n:
                f_past = f_past[1:] + [f(t[i+1], y[i+1])]

        return t, y


class PredictorCorrectorIVPSolver(IVPSolver):

    def __init__(self, *args, h: float, order: int = 1, **kwargs):

        super().__init__(*args, **kwargs)

        method_ab = butcher_tableaus.MULTISTEP_METHODS['adams-bashforth']
        if order not in method_ab:
            raise ValueError(f"Order {order} is not supported.")

        self.ns_ab = method_ab[order]['ns']
        self.b_ab = method_ab[order]['b']

        method_am = butcher_tableaus.MULTISTEP_METHODS['adams-moulton']
        if order not in method_am:
            raise ValueError(f"Order {order} is not supported.")

        self.ns_am = method_am[order]['ns']
        self.b_am = method_am[order]['b']

        self.order = order
        self.h = h

    def solve(self) -> tuple[np.ndarray[tuple[int]], np.ndarray[tuple[int, int]]]:

        t0, tf, h = self.t0, self.tf, self.h
        y0 = self.y0
        f = self.f
        order, ns_ab, b_ab = self.order, self.ns_ab, self.b_ab
        ns_am, b_am = self.ns_am, self.b_am

        n = int(round((tf-t0)/h))
        if n < order:
            raise ValueError(f"Order {order} larger than n.")

        neq = y0.shape[0]

        t = np.array([t0+i*h for i in range(n+1)])

        y = np.zeros((n+1, neq))

        y[0,:] = y0.copy()

        if order > 1:
            rk4 = ExplicitIVPSolver(f, t0, t[order-1], y0, 'rk4', h=h)
            _, y[:order,:] = rk4.solve()

        f_past = [f(t[i], y[i,:]) for i in range(order)]

        for i in range(order-1, n):

            y_ab = np.zeros(neq)
            for k in range(order):
                y_ab += b_ab[k]*f_past[-(k+1)]

            y_predict = y[i,:] + h*y_ab[:]

            f_predict = f(t[i+1], y_predict)

            y_am = b_am[0]*f_predict
            for k in range(1, order):
                y_am += b_am[k]*f_past[-k]

            y[i+1,:] = y[i,:] + h*y_am[:]

            if i + 1 < n:
                f_past = f_past[1:] + [f(t[i+1], y[i+1])]

        return t, y
