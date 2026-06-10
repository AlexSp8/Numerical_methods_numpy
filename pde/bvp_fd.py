
from typing import Callable
import numpy as np

from linear_systems.linear_solver import LinearSolver
from non_linear_systems import newton_solver, non_linear_problem
from ode.time_integration import TimeIntegration
from ode import plot


class BVPSolver():

    def __init__(self):
        pass

    def set_ls_solver(self, ls_solver: LinearSolver):

        self.ls_solver = ls_solver

    def set_nr_solver(self, u0: np.ndarray[tuple[int]],
        k_max: int = 100, tol: float = 1e-8, r: float = 1.0):

        self.nr_solver = newton_solver.NewtonSolver(self.ls_solver, u0, k_max, tol, r)

    def set_time_int(self, t0: float, tf: float, dt: float,
        dt_min: float, dt_max: float, atol: float, rtol: float,
        u0: np.ndarray[tuple[int]]):

        self.time_int = TimeIntegration(t0, tf, dt, dt_min, dt_max, atol, rtol, u0)

    def set_problem(self, f_res: Callable, bc: dict[str, Callable], theta: float = 1.0):

        self.problem = non_linear_problem.TransientBVP1DFD(
            f=f_res, x=self.x, bc=bc, time_int=self.time_int, theta=theta)

    def solve(self, dtw: float = 1.0):

        u = self.nr_solver.u0

        nunknowns, nnodes = u.shape[0], self.x.shape[0]
        neq = round(nunknowns/nnodes)

        dt = self.time_int.dt

        tw = self.time_int.t

        while self.time_int.t < self.time_int.tf:

            self.time_int.print_info()

            u = self.nr_solver.solve(self.problem, output=True)

            self.nr_solver.update_guess(u)

            dt = self.time_int.update_time_step(u)

            self.time_int.update(dt, u)

            if abs(self.time_int.t-tw) > dtw:
                tw = self.time_int.t
                for ieq in range(neq):
                    ui_num = u[ieq::neq].copy()
                    plot.plot_ode_system_evolution(self.x, [ui_num])

        return u
    
    def create_mesh(self, nnodes: int, xd: np.ndarray[tuple[int]],
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

        self.x = x
    
    def check_dt(self, D0: float, v0: float, dt0: float) -> float:
        """Checks if the time discretization is stable.

        Args:
            x: the spatial domain
            D0: the diffusion number
            v0: the convection velocity
            dt0: the initial time step to check
        Returns:
            The initial time step if it is stable. Otherwise the maximum stable dt."""

        dt = dt0

        D = np.maximum(np.abs(D0), 1e-12)
        v = np.maximum(np.abs(v0), 1e-12)

        # dx = np.diff(self.x)
        dx = self.x[1:] - self.x[0:-1]
        dx_min = np.min(dx)
        dx_max = np.max(dx)

        dt_max_diffusion = (dx_min**2)/(2*D)
        dt_max_convection = (2*D)/(v**2)
        dt_max = min(dt_max_diffusion, dt_max_convection)

        cfl_diffusion = dt/(dx_min**2)
        courant = (v*dt)/dx_min

        print("="*40)
        print("    CONVECTION-DIFFUSION DIAGNOSTICS")
        print("="*40)
        print(f"dx_min:              {dx_min:.4f}")
        print(f"dx_max:              {dx_max:.4f}\n")
        print(f"Courant (C):         {courant:.4f}\n")
        print(f"Max dt (Diffusion):  {dt_max_diffusion:.6f}")
        print(f"Max dt (Convection): {dt_max_convection:.6f}")
        print(f"Max dt Allowed:      {dt_max:.6f}")
        print(f"Time Step (dt):      {dt:.6f}")
        print("-"*40)

        if dt > dt_max:
            print(f"dt = {dt:.6f} exceeds the maximum stable dt_max = {dt_max:.6f}.")
            print(f"Setting dt = {dt_max:.6f}.")
            dt = dt_max
            courant = (v*dt)/dx_min
            print(f"Courant (C):         {courant:.4f}\n")

        print("="*40)

        return dt

    def check_dx(self, D0: float, v0: float) -> None:
        """Checks if the spatial discretization is stable.

        Args:
            x: the spatial domain
            D0: the diffusion number
            v0: the convection velocity
        Returns:
            The initial time step if it is stable. Otherwise the maximum stable dt."""

        D = np.maximum(np.abs(D0), 1e-12)
        v = np.maximum(np.abs(v0), 1e-12)

        # dx = np.diff(self.x)
        dx = self.x[1:] - self.x[0:-1]
        dx_max = np.max(dx)

        pe = (v*dx_max)/D

        dx_stable = D*2.0/v
        if pe > 2.0:
            print(f"dx = {dx_max:.4f} > dx_max = {dx_stable:.4f} (Pe = {pe:.1f} > 2.0)")
            raise ValueError('Too large dx for this system')
