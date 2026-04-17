
import math
import numpy as np
import numpy.typing as npt

from non_linear_systems import newton_solver, non_linear_problem
from linear_systems import direct_solver

def F(x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Test system of non-linear algebraic equations in residual form."""
    return np.array([
        x[0]**2 + x[1]**2 - 4,
        math.exp(x[0]) + x[1] - 1
    ])

def non_linear_systems():

    print('Non-Linear Systems')

    u0 = np.array([1.0, -1.0])
    # u0 = np.array([-2.0, 1.0])

    print('\nNewton-Raphson')
    ls_solver = direct_solver.LUSolver()
    nr_solver = newton_solver.NewtonSolver(ls_solver, u0=u0, k_max=100, tol=1e-8, r=1.0)
    problem = non_linear_problem.Regular(F)
    print('x =', nr_solver.solve(problem, output=True))

if __name__ == '__main__':
    non_linear_systems()
