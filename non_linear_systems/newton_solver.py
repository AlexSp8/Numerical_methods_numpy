
import time
import numpy as np
import scipy

# from utilities import matrix_operations
from non_linear_systems.non_linear_problem import NonlinearProblem
from linear_systems.linear_solver import LinearSolver

class NewtonSolver:

    def __init__(self, ls_solver: LinearSolver, u0: NDArray[np.float64],
        k_max: int = 1000, tol: float = 1e-8, r: float = 1.0):

        self.ls_solver = ls_solver
        self.u0 = u0
        self.k_max = k_max
        self.tol = tol
        self.r = r
        self.is_full = True
        self.jac = None

    def update_guess(self, u0: NDArray[np.float64]):

        self.u0 = u0

    def get_norms(self, du: NDArray[np.float64], res: NDArray[np.float64]
        ) -> tuple[float, float]:
        """Returns the norms of the correction and residual vectors.

        Args:
            du: the correction vector
            res: the residual vector
        Returns:
            The correction and residual norms."""

        cor_norm = np.linalg.norm(du)
        res_norm = np.linalg.norm(res)

        return cor_norm, res_norm

    def solve(self, problem: NonlinearProblem, output: bool = False,
        ls_solver: LinearSolver = None) -> NDArray[np.float64]:
        """Returns the solution of a non-linear system of algebraic equations
        around an initial guess using the Newton-Raphson method.

        Args:
            problem: the non-linear problem to be solved (object containing the residual method)
            output (optional): flag for iteration info output
            ls_solver (optional): a linear solver object to solve the linearized system
        Returns:
            The solution vector of the non-linear system, u."""

        # u = scipy.optimize.fsolve(problem.res, self.u0)
        # return u

        if ls_solver is not None:
            solver = ls_solver
        else:
            solver = self.ls_solver

        u = self.u0.copy()

        n, neq = u.shape[0], u.shape[0]

        res = np.zeros(neq)

        is_full = self.is_full
        if self.jac is None:
            self.jac = np.zeros((neq,n))
            is_full = True

        for k in range(1, self.k_max+1):

            t_start = time.perf_counter()

            res = problem.f_res(u)

            if is_full:
                self.jac = problem.jacobian(u, res)
                # for i in range(n):
                #     print(" ".join(f"{x:.2e}" for x in jac[i, :]))

            du = solver.solve(self.jac, -res)
            # solver.x0 = du

            cor_norm, res_norm = self.get_norms(du, res)
            t_end = time.perf_counter()
            if output:
                dt = t_end - t_start
                print(f'k = {k}, Res Norm: {res_norm:.4e}, Cor Norm: {cor_norm:.4e}, Time = {dt:.2f} s')
                # print(matrix_operations.condition_number(jac))

            if k == 1:
                cor_norm0 = cor_norm

            if (res_norm < self.tol) and (cor_norm < self.tol):
                break

            if cor_norm < np.sqrt(self.tol):
                is_full = False
            else:
                is_full = True

            u = u + self.r*du

        if k >= self.k_max:
            print(f"Warning: Newton maximum iterations ({k}) reached without converging.")

        if cor_norm0 < np.sqrt(self.tol):
            self.is_full = False
    
        return u

class LevenbergMarquardtSolver(NewtonSolver):

    def __init__(self, *args, l0: float = 1e-3, scale: float = 10.0, **kwargs):

        super().__init__(*args, **kwargs)
        self.l0 = l0
        self.scale = scale

    def solve(self, problem: NonlinearProblem,
        output: bool = False, ls_solver: LinearSolver = None
        ) -> NDArray[np.float64]:
        """Returns the solution of a non-linear system of algebraic equations
        around an initial guess using the Levenberg-Marquardt method.

        Args:
            problem: the non-linear problem to be solved (object containing the residual method)
            output (optional): flag for iteration info output
            ls_solver (optional): a linear solver object to solve the linearized system
        Returns:
            The solution vector of the non-linear system, u."""

        # u = scipy.optimize.fsolve(problem.res, self.u0)
        # return u

        if ls_solver is not None:
            solver = ls_solver
        else:
            solver = self.ls_solver

        u = self.u0.copy()

        res = problem.f_res(u)

        ssr = np.sum(res**2)

        n, neq = self.u0.shape[0], res.shape[0]

        damp, s = self.l0, self.scale

        is_full = self.is_full
        is_full = True

        for k in range(1, self.k_max+1):

            if is_full:
                jac = problem.jacobian(u, res)

            # print(matrix_operations.condition_number(jac))

            # if m != n:
            jt = jac.T
            j_mat = np.matmul(jt, jac)
            b = -np.matmul(jt, res)
            # else:
            #     j_mat = jac
            #     b = -res

            diag_j = np.diag(np.diag(j_mat))

            if np.any(np.diag(diag_j) == 0):
                diag_j = np.eye(n)

            for k2 in range(1, 21):

                try:
                    jtj_damped = j_mat + damp*diag_j # np.eye(n)
                    du = solver.solve(jtj_damped, b)
                except:
                    damp *= s
                    continue

                # solver.x0 = du

                u_trial = u + self.r*du

                res_trial = problem.f_res(u_trial)

                ssr_trial = np.sum(res_trial**2)

                if ssr_trial < ssr:
                    u, ssr, res = u_trial, ssr_trial, res_trial
                    damp /= s  # Move towards Newton
                    break
                else:
                    damp *= s  # Move towards Steepest Descent

                if damp > 1e12:
                    break

            cor_norm, res_norm = self.get_norms(du, res)
            if output:
                print(f'k = {k} (Inner k2: {k2}), Res Norm: {res_norm:.4e}, Cor Norm: {cor_norm:.4e}')

            if (res_norm < self.tol) or (cor_norm < self.tol):
                return u

            if cor_norm < np.sqrt(self.tol):
                is_full = False

        print(f"Warning: Maximum iterations ({self.k_max}) reached without converging.")

        return u
