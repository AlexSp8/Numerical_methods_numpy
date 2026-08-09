
from typing import Callable
import numpy as np
from numpy.typing import NDArray

from linear_systems.linear_solver import LinearSolver
from non_linear_systems import newton_solver, non_linear_problem
from ode.time_integration import TimeIntegration
from ode import plot
from integration import fe_gauss_integration
from pde.fe_stabilization import FEStabilization
from pde.mesh_discretization import MeshDiscretization
from pde.boundary_conditions import BoundaryConditions


class BVPSetup():

    def __init__(self, neq: int, mesh: MeshDiscretization, boundary: BoundaryConditions,
        f_res: Callable):

        self.neq = neq
        self.mesh = mesh
        self.boundary = boundary
        self.f_res = f_res

    def set_ls_solver(self, ls_solver: LinearSolver):

        self.ls_solver = ls_solver

    def set_nr_solver(self, u0: NDArray[np.float64],
        k_max: int = 100, tol: float = 1e-8, r: float = 1.0):

        self.nr_solver = newton_solver.NewtonSolver(self.ls_solver, u0, k_max, tol, r)

    def set_time_int(self, t0: float, tf: float, dt: float,
        dt_min: float, dt_max: float, atol: float, rtol: float,
        u0: NDArray[np.float64]):

        self.time_int = TimeIntegration(t0, tf, dt, dt_min, dt_max, atol, rtol, u0)

    def solve(self, dtw: float = 1.0):

        neq = self.neq

        u = self.nr_solver.u0

        dt = self.time_int.dt

        tw = self.time_int.t

        while True:

            self.time_int.print_info()

            u = self.nr_solver.solve(self.problem, output=True)

            if abs(self.time_int.t-tw) > dtw:
                tw = self.time_int.t
                for ieq in range(neq):
                    ui_num = u[ieq::neq].copy()
                    plot.plot_ode_system_evolution(self.mesh.x_mesh, [ui_num])

            if abs(self.time_int.t - self.time_int.tf) < 1e-12:
                break

            self.nr_solver.update_guess(u)

            dt = self.time_int.update_time_step(u)

            self.time_int.update(dt, u)

        return u


class BVPSetupFD(BVPSetup):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def set_fd_problem(self, theta: float = 1.0):

        self.problem = non_linear_problem.FiniteDifferencesBVP(
            f=self.f_res, neq=self.neq, mesh=self.mesh, boundary=self.boundary,
            time_int=self.time_int, theta=theta)


class BVPSetupFE(BVPSetup):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def set_fe_gauss_int(self, ng: int):

        nbf = self.mesh.nbf
        fe_type = self.mesh.fe_type

        if fe_type == 'line':
            self.fe_gauss_int = fe_gauss_integration.FEGaussIntegrationLine(nbf=nbf, ng=ng)
        elif fe_type == 'quadrangle':
            self.fe_gauss_int = fe_gauss_integration.FEGaussIntegrationQuadrangle(nbf=nbf, ng_1d=ng)
        elif fe_type == 'triangle':
            self.fe_gauss_int = fe_gauss_integration.FEIntegrationTriangle(nbf=nbf, ng=ng)
        elif fe_type == 'hexahedron':
            self.fe_gauss_int = fe_gauss_integration.FEGaussIntegrationHexahedron(nbf=nbf, ng_1d=ng)
        elif fe_type == 'tetrahedron':
            self.fe_gauss_int = fe_gauss_integration.FEIntegrationTetrahedron(nbf=nbf, ng=ng)
        else:
            raise ValueError('Unsupported fe_type in set_fe_gauss_int!')

    def set_fe_stabilization(self, f_quantities: Callable, f_strong: Callable,
        f_operator: Callable):

        fe_order = self.mesh.fe_order
        if fe_order == 'linear':
            order = 1
        elif fe_order == 'quadratic' or fe_order == 'serendipity':
            order = 2
        else:
            raise ValueError('Unsupported fe_order in set_fe_stabilization!')

        self.fe_stab = FEStabilization(f_quantities=f_quantities,
                        f_strong=f_strong, f_operator=f_operator, order=order)

    def set_fe_problem(self, theta: float = 1.0):

        self.problem = non_linear_problem.FiniteElementsBVP(f=self.f_res, neq=self.neq, mesh=self.mesh,
           boundary=self.boundary, time_int=self.time_int, gauss_int=self.fe_gauss_int, fe_stab=self.fe_stab,
           theta=theta)
