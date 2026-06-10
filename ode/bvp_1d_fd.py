
from typing import Callable

import numpy as np

from linear_systems import direct_solver
from non_linear_systems import newton_solver
from non_linear_systems.non_linear_problem import TransientBVP1DFD
from ode.time_integration import TimeIntegration

def create_mesh(nnodes: int, xd: np.ndarray[tuple[int]],
    p: int = 1) -> np.ndarray[tuple[int]]:
    """Returns the mesh coordinates given the number of nodes,
    the domain bounds and the stretching factor for power law
    refinement. p > 1: refinement towards x0,
    0 < p < 1: refinement towards xf.

    Args:
        nnodes: the number of nodes
        xd: the domain bounds
        p (optional): the power law refinement factor
            For p > 1: refinement towards x0.
            For 0 < p < 1: refinement towards xf
    Returns:
        The array of mesh coordinates."""

    x0, xf = xd[0], xd[1]
    L = xf - x0

    x = np.zeros(nnodes)
    for i in range(nnodes):
        xi = i/(nnodes-1)
        gi = xi**p
        # gi = 0.5*(1.0-np.cos(p*np.pi*xi))
        # gi = np.tanh(p*(2*xi-1))/(2*np.tanh(p)) + 0.5
        x[i] = x0 + L*gi

    return x

def solve_bvp_1d_fd(
    f: Callable[[float, np.ndarray[tuple[int]], np.ndarray[tuple[int]], np.ndarray[tuple[int]]],
                np.ndarray[tuple[int]]], x: np.ndarray[tuple[int]], bc: dict[str, Callable],
    neq: int = 1, u_guess: np.ndarray[tuple[int]] = None) -> np.ndarray[tuple[int]]:
    """Returns the solution of the non-linear system of algebraic equations
    that results from Finite Differences (FD) approximation of a system of second order 1D-BVP
    equations using the Newton-Raphson method.

    Args:
        f: the differential equation in residual form
        x: the discretized domain array
        bc: the dictionary containing the boundary conditions functions
        neq: the number of equations to be solved
    Returns:
        The FD solution to the system of ODEs (2nd order BVP-1D) at the nodes of the domain."""

    nnodes = x.shape[0]

    nunknowns = nnodes*neq

    ls_solver = direct_solver.LUSolver()

    if u_guess is None:
        u0 = np.ones(nunknowns)
    else:
        u0 = u_guess.copy()

    nr_solver = newton_solver.NewtonSolver(ls_solver, u0=u0, k_max=100, tol=1e-8, r=1.0)

    t0, tf, dt0, dt_min, dt_max, atol, rtol = 0.0, 0.0, float('inf'), float('inf'), float('inf'), 0.0, 0.0
    time_int = TimeIntegration(t0, tf, dt0, dt_min, dt_max, atol, rtol, u0)
    problem = TransientBVP1DFD(f=f, x=x, bc=bc, time_int=time_int)

    return nr_solver.solve(problem, output=True)
