
from typing import Callable, Tuple
import numpy as np
import numpy.typing as npt
import scipy

# from utilities import matrix_operations
from non_linear_systems.non_linear_problem import NonlinearProblem
from linear_systems.linear_solver import LinearSolver

class NewtonSolver:

    def __init__(self, ls_solver: LinearSolver, u0: npt.NDArray[np.float64],
        k_max: int = 1000, tol: float = 1e-8, r: float = 1.0):

        self.ls_solver = ls_solver
        self.u0 = u0
        self.k_max = k_max
        self.tol = tol
        self.r = r

    def get_norms(self, du: npt.NDArray[np.float64], res: npt.NDArray[np.float64]
        ) -> Tuple[float, float]:
        """Returns the correction and residual norms."""
        cor_norm = np.linalg.norm(du)
        res_norm = np.linalg.norm(res)
        return cor_norm, res_norm

    def is_converged(self, cor_norm: float, res_norm: float, mode: str = 'and') -> bool:
        """Checks the correction and residual norm convergence criteria."""
        if mode == 'and':
            return (res_norm < self.tol) and (cor_norm < self.tol)
        else:
            return (res_norm < self.tol) or (cor_norm < self.tol)

    def solve(self, problem: NonlinearProblem,
        output: bool = False, ls_solver: LinearSolver = None) -> npt.NDArray[np.float64]:
        """Returns the solution of a non-linear system of algebraic equations, F,
        around an initial guess, u0, using the Newton-Raphson method"""

        # u = scipy.optimize.fsolve(problem.res, self.u0)
        # return u

        solver = ls_solver or self.ls_solver

        u = self.u0.copy()
        for k in range(1, self.k_max+1):

            res = problem.f_res(u)

            if hasattr(problem, 'jac'):
                jac = problem.jac(u)
            else:
                jac = self.jacobian(res, problem.res, u)

            # print(matrix_operations.condition_number(jac))

            du = solver.solve(jac, -res)
            # solver.x0 = du

            cor_norm, res_norm = self.get_norms(du, res)
            if output:
                print(f'k = {k}, Res Norm: {res_norm:.4e}, Cor Norm: {cor_norm:.4e}')

            if self.is_converged(cor_norm, res_norm, 'and'):
                break

            u = u + self.r*du

        return u

    @staticmethod
    def jacobian(res: npt.NDArray[np.float64],
        f_res: Callable[[npt.NDArray[np.float64]], npt.NDArray[np.float64]],
        u: npt.NDArray[np.float64], h: float = 1e-8) -> npt.NDArray[np.float64]:
        """Returns the Jacobian of a system of non-linear algebraic equations
        around the values, u, of the unknowns."""

        n, m = len(u), len(res)

        jac = np.zeros((m,n))
        for j in range(n):

            u_val = u[j]
            u[j] += h
            res_f = f_res(u)
            u[j] = u_val

            # u_val = u[j]
            # u[j] -= h
            # res_b = f_res(u)
            # u[j] = u_val

            jac[:,j] = (res_f-res)/h
            # jac[:,j] = (res_f-res_b)/(2*h)

        return jac

class MarquardtSolver(NewtonSolver):

    def __init__(self, *args, l0: float = 1e-3, scale: float = 10.0, **kwargs):

        super().__init__(*args, **kwargs)
        self.l0 = l0
        self.scale = scale

    def solve(self, problem: NonlinearProblem,
        output: bool = False, ls_solver: LinearSolver = None) -> npt.NDArray[np.float64]:
        """Returns the solution of a non-linear system of algebraic equations, F,
        around an initial guess, u0, using the Marquardt method"""

        # u = scipy.optimize.fsolve(problem.res, self.u0)
        # return u

        solver = ls_solver or self.ls_solver

        s = self.scale

        n = len(self.u0)

        u = self.u0.copy()

        for k in range(1, self.k_max+1):

            res = problem.f_res(u)

            if hasattr(problem, 'jac'):
                jac = problem.jac(u)
            else:
                jac = self.jacobian(res, problem.res, u)

            # print(matrix_operations.condition_number(jac))

            l = self.l0
            for k2 in range(1, 21):

                try:
                    jac_damped = jac + l*np.eye(n)
                    du = solver.solve(jac_damped, -res)
                except:
                    l *= s
                    continue

                # solver.x0 = du

                u_trial = u + self.r*du

                res_trial = problem.f_res(u_trial)

                is_better = (np.linalg.norm(res_trial) < np.linalg.norm(res))
                if is_better:
                    u = u_trial
                    res = res_trial
                    l /= s  # Move closer to pure Newton
                    break
                else:
                    l *= s  # Move closer to Steepest Descent

                if l > 1e12:
                    break

            cor_norm, res_norm = self.get_norms(du, res)
            if output:
                print(f'k = {k} (k2: {k2}), Res Norm: {res_norm:.4e}, Cor Norm: {cor_norm:.4e}')

            if self.is_converged(cor_norm, res_norm, 'and'):
                break

        return u
