
from typing import Callable

import numpy as np

def solve_parabolic_1d(
    f: Callable[[float, np.ndarray[tuple[int]], np.ndarray[tuple[int]]], np.ndarray[tuple[int]]],
    x: np.ndarray[tuple[int]], y0: np.ndarray[tuple[int]], neq: int,
    ode_solver: Callable[[Callable, float, float, np.ndarray[tuple[int]], float],
                         tuple[np.ndarray[tuple[int]], np.ndarray[tuple[int, int]]]]
    ) -> np.ndarray[tuple[int]]:
    """Returns the solution of the non-linear system of algebraic equations
    that results from FD approximation of a system of second order 1D-BVP
    problems using the Newton-Raphson method.

    Args:
        f: the production term function
        x: the domain array discretized with Finite Differences
        y0: the initial conditions
        bc: the dictionary describing the boundary conditions
        neq: the number of equations to be solved
    Returns:
        The solution to the system of ODEs (2nd order BVP-1D) at the nodes of the domain."""

    problem = Parabolic1D(f=f, x=x, neq=neq)

    # t, y = ode_solver(problem.g, t0=0.0, tf=100.0, y0=y0, h=0.5)
    # t, y = ode_solver(problem.g, t0=0.0, tf=100.0, y0=y0, h_min=1e-6, h_max=1.0)

    # t, y = ode_solver(problem.g, t0=0.0, tf=80.0, y0=y0, h=1.0, method='lobatto-iiia-3')
    t, y = ode_solver(problem.g, t0=0.0, tf=80.0, y0=y0, method='euler')

    return t, y


class Parabolic1D():

    def __init__(self,
        f: Callable[[np.ndarray[tuple[int]]], np.ndarray[tuple[int]]],
        x: np.ndarray[tuple[int]], neq: int):

        self.f = f
        self.x = x
        self.neq = neq

    def g(self, t: float, u: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:

        nnodes, nunknowns = (self.x).shape[0], u.shape[0]

        neq = self.neq

        res = np.zeros(nunknowns)

        for i in range(1, nnodes-1):

            hf = self.x[i+1] - self.x[i]
            hb = self.x[i] - self.x[i-1]

            t2_hf_hb = 2.0/(hf+hb)

            row1 = i*neq

            ui = u[row1:row1+neq]

            # 1st, 2nd derivatives
            dui = np.zeros(neq)
            d2ui = np.zeros(neq)
            for jeq in range(neq):

                row = row1 + jeq

                dui[jeq] = (u[row+neq] - u[row])/hf
                # dui[jeq] = (u[row+neq] - u[row-neq])/(hb+hf)

                d2ui[jeq] = (u[row+neq] - u[row])/hf
                d2ui[jeq] -= (u[row] - u[row-neq])/hb
                d2ui[jeq] *= t2_hf_hb

            # source term
            # fi = self.f(self.x[i], ui, dui) # explicit form
            fi = self.f(self.x[i], ui, dui, d2ui) # implicit form

            for jeq in range(neq):
                row = row1 + jeq
                # res[row] = d2ui[jeq] + fi[jeq] # explicit form
                res[row] = fi[jeq] # implicit form

        # boundaries
        for j in range(neq):

            res[j] = 0.0

            row = (nnodes-1)*neq+j
            res[row] = 0.0

        return res
