
from typing import Callable

import numpy as np

from linear_systems import direct_solver
from non_linear_systems import newton_solver, non_linear_problem

class FiniteDifferencesMOL1D():
    """The method of lines (MOL) is based on discretizing the spatial dimensions
    and solving a system of ODEs, dudt = f(t,x,u,u',u"), involving all the unknowns."""

    def __init__(self,
        f: Callable[[float, float, np.ndarray[tuple[int]],
                     np.ndarray[tuple[int]], np.ndarray[tuple[int]]],
                     np.ndarray[tuple[int]]],
        x: np.ndarray[tuple[int]], neq: int, bc: dict[str, Callable]):

        self.f = f
        self.x = x
        self.neq = neq
        self.bc = bc

        self.is_Dirichlet_left, self.is_Dirichlet_right = self._classify_bc()

    def _classify_bc(self):
        """Checks which boundary conditions involve derivatives of the unknowns
        and which are of Dirichlet type.

        Args:
            self: the FiniteDifferencesMOL1D class object
        Returns:
            Two boolean arrays for the BCs on the left and right boundaries, respectively.
            True indicates a Dirichlet type BC and False indicates a derivative type BC."""

        t0, x0, xf = 1.0, self.x[0], self.x[-1]
        neq = self.neq

        bc_left = self.bc['left']
        bc_right = self.bc['right']

        ub_dum = np.ones(neq)
        res_left = bc_left(t0, x0, ub_dum, np.zeros(neq))
        res_right = bc_right(t0, xf, ub_dum, np.zeros(neq))

        is_Dirichlet_left = np.zeros(neq, dtype=bool)
        is_Dirichlet_right = np.zeros(neq, dtype=bool)

        dudx_p = np.full(neq, 1e-8)

        res_p = bc_left(t0, x0, ub_dum, dudx_p)
        for i in range(neq):
            if np.abs(res_p[i] - res_left[i]) < 1e-12:
                is_Dirichlet_left[i] = True

        res_p = bc_right(t0, xf, ub_dum, dudx_p)
        for i in range(neq):
            if np.abs(res_p[i] - res_right[i]) < 1e-12:
                is_Dirichlet_right[i] = True

        # print(is_Dirichlet_left, is_Dirichlet_right)

        return is_Dirichlet_left, is_Dirichlet_right

    def g(self, t: float, u: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
        """Returns the right hand side of the parabolic problem
        dudt = f(t, x, u, dudx, d2udx2) at time t.

        Args:
            t: the time instant
            u: the vector of unknowns
        Returns:
            An array of the right hand side of the parabolic problem
            for all unknowns at time t."""

        f = self.f
        x = self.x
        neq = self.neq

        nnodes, nunknowns = x.shape[0], u.shape[0]

        dudt = np.zeros(nunknowns)

        for i in range(1, nnodes-1):

            hf = x[i+1] - x[i]
            hb = x[i] - x[i-1]
            t2_hf_hb = 2.0/(hf+hb)

            row1 = i*neq
            row_f = row1+neq

            ui_b = u[row1-neq:row1]
            ui = u[row1:row_f]
            ui_f = u[row_f:row_f+neq]

            # duidx = (ui_f - ui)/hf # Forward
            # duidx = (ui - ui_b)/hb # Backward
            duidx = (ui_f - ui_b)/(hb+hf) # Central

            d2uidx2 = ((ui_f-ui)/hf - (ui-ui_b)/hb)*t2_hf_hb

            # dudt[row1:row_f] = d2uidx2[:] +f(t, x[i], ui, duidx) # explicit form
            dudt[row1:row_f] = f(t, x[i], ui, duidx, d2uidx2) # implicit form

        # boundaries
        dudt[:neq] = self.left_bc_residual(t, u)
        dudt[-neq:] = self.right_bc_residual(t, u)

        return dudt

    def left_bc_residual(self, t: float, u: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
        """Returns the time derivatives on the left boundary of a 1D domain.

        Args:
            t: the time instant
            u: the vector of unknowns
        Returns:
            The time derivatives of all unknowns on the left boundary
            after applying the boundary conditions."""

        bc = self.bc['left']
        f = self.f
        x = self.x
        neq = self.neq
        is_Dirichlet = self.is_Dirichlet_left

        hf = x[1] - x[0]
        ub = u[0:neq]
        u_f = u[neq:2*neq]

        u_ghost = np.zeros(neq)

        def res_bc(u_ghost):
            dudx_guess = (u_f - u_ghost)/(2.0*hf)
            res = bc(t, x[0], ub, dudx_guess)
            for j in range(neq):
                if is_Dirichlet[j]:
                    res[j] = u_ghost[j] - (2.0*ub[j]-u_f[j])
            return res

        u0 = ub.copy()
        for j in range(neq):
            if is_Dirichlet[j]:
                u0[j] = (2.0*ub[j]-u_f[j])
        ls_solver = direct_solver.LUSolver()
        nr_solver = newton_solver.NewtonSolver(ls_solver, u0=u0, k_max=100, tol=1e-8, r=1.0)
        problem = non_linear_problem.Regular(res_bc)
        u_ghost = nr_solver.solve(problem, output=False)

        dudx_b = (u_f - u_ghost)/(2.0*hf)
        d2udx2_b = (u_f - 2.0*ub + u_ghost)/(hf**2)

        dudt_f = f(t, x[0], ub, dudx_b, d2udx2_b)

        dudt_b = np.zeros(neq)

        res_left = bc(t, x[0], ub, np.zeros(neq))
        for j in range(neq):

            if is_Dirichlet[j]:
                # dudt = -(dR/dt)/(dR/du)
                # dR/dt
                dt_pert = 1e-8
                res_p = bc(t+dt_pert, x[0], ub, np.zeros(neq))
                dR_dt = (res_p[j] - res_left[j])/dt_pert

                # dR/du
                du_pert = 1e-8
                ub_pert = ub.copy()
                ub_pert[j] += du_pert
                res_p = bc(t, x[0], ub_pert, np.zeros(neq))
                dR_du = (res_p[j] - res_left[j])/du_pert

                penalty = res_left[j]*10.0

                if np.abs(dR_du) > 1e-12:
                    dudt_b[j] = -(dR_dt+penalty)/dR_du

            else:
                dudt_b[j] = dudt_f[j]

        return dudt_b

    def right_bc_residual(self, t: float, u: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
        """Returns the time derivatives on the right boundary of a 1D domain.

        Args:
            t: the time instant
            u: the vector of unknowns
        Returns:
            The time derivatives of all unknowns on the right boundary
            after applying the boundary conditions."""

        bc = self.bc['right']
        f = self.f
        x = self.x
        neq = self.neq
        is_Dirichlet = self.is_Dirichlet_right

        nnodes = x.shape[0]

        hb = x[nnodes-1] - x[nnodes-2]
        ub = u[-neq:]
        u_b = u[nnodes-2*neq:nnodes-neq]

        u_ghost = np.zeros(neq)

        def res_bc(u_ghost):
            dudx_guess = (u_ghost - u_b)/(2.0*hb)
            res = bc(t, x[nnodes-1], ub, dudx_guess)
            for j in range(neq):
                if is_Dirichlet[j]:
                    res[j] = u_ghost[j] - (2.0*ub[j]-u_b[j])
            return res

        u0 = ub.copy()
        for j in range(neq):
            if is_Dirichlet[j]:
                u0[j] = (2.0*ub[j]-u_b[j])
        ls_solver = direct_solver.LUSolver()
        nr_solver = newton_solver.NewtonSolver(ls_solver, u0=u0, k_max=100, tol=1e-8, r=1.0)
        problem = non_linear_problem.Regular(res_bc)
        u_ghost = nr_solver.solve(problem, output=False)

        dudx_b = (u_ghost - u_b)/(2.0*hb)
        d2udx2_b = (u_ghost - 2.0*ub + u_b)/(hb**2)

        dudt_f = f(t, x[nnodes-1], ub, dudx_b, d2udx2_b)

        dudt_b = np.zeros(neq)

        res_right = bc(t, x[nnodes-1], ub, np.zeros(neq))
        for j in range(neq):

            if is_Dirichlet[j]:
                # dudt = -(dR/dt)/(dR/du)
                # dR/dt
                dt_pert = 1e-8
                res_p = bc(t+dt_pert, x[nnodes-1], ub, np.zeros(neq))
                dR_dt = (res_p[j] - res_right[j])/dt_pert

                # dR/du
                du_pert = 1e-8
                ub_pert = ub.copy()
                ub_pert[j] += du_pert
                res_p = bc(t, x[nnodes-1], ub_pert, np.zeros(neq))
                dR_du = (res_p[j] - res_right[j])/du_pert

                penalty = res_right[j]*10.0

                if np.abs(dR_du) > 1e-12:
                    dudt_b[j] = -(dR_dt+penalty)/dR_du

            else:
                dudt_b[j] = dudt_f[j]

        return dudt_b
