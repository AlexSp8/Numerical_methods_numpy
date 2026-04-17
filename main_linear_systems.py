
import numpy as np

from utilities import matrix_operations
from linear_systems import direct_solver, iterative_solver

def main_linear_systems():

    print('Linear Systems')

    print('\nDirect solvers')

    A = np.array([
        [0.3, 0.52, 1.0],
        [0.5, 1.0, 1.9],
        [0.1, 0.3, 0.5]
    ])
    b = np.array([-0.01, 0.67, -0.44])

    print('\nCramer rule 3x3: ', end='')
    cramer = direct_solver.CramerSolver()
    print(cramer.solve(A, b))

    print('\nGauss Elimination: ', end='')
    gauss = direct_solver.GaussSolver()
    print(gauss.solve(A, b))

    print('\nLU decomposition: ', end='')
    lu = direct_solver.LUSolver()
    print(lu.solve(A, b))

    print('\nNumPy solve ', end='')
    print(np.linalg.solve(A,b))

    print('\nThomas: ', end='')
    A = np.array([
        [2.04, -1.0, 0.0, 0.0],
        [-1.0, 2.04, -1.0, 0.0],
        [0.0, -1.0, 2.04, -1.0],
        [0.0, 0.0, -1.0, 2.04]
    ])
    b = np.array([40.8, 0.8, 0.8, 200.8])
    l, d, u = matrix_operations.tri_diagonal(A)
    print(direct_solver.thomas_tri_diagonal(l, d, u, b))

    print('\nCholesky decomposition:')
    A = np.array([
        [6.0, 15.0, 55.0],
        [15.0, 55.0, 225.0],
        [55.0, 225.0, 979.0]
    ])
    print(matrix_operations.cholesky_decomposition(A))

    print('\nIterative solvers')
    A = np.array([
        [3.0, -0.1, -0.2],
        [0.1, 7.0, -0.3],
        [0.3, -0.2, 10.0]
    ])
    b = np.array([7.85, -19.3, 71.4])

    print('\nJacobi')
    jacobi = iterative_solver.JacobiSolver(x0=None, k_max=100, tol=1e-8)
    print(jacobi.solve(A, b, output=True))

    print('\nSOR')
    sor = iterative_solver.SORSolver(x0=None, k_max=100, tol=1e-8, w=1.5)
    print(sor.solve(A, b, output=True))

    print('\nSteepest Descent')
    sd = iterative_solver.SDSolver(x0=None, k_max=100, tol=1e-8)
    print(sd.solve(A, b, output=True))

    print('\nConjugate Gradient')
    cg = iterative_solver.CGSolver(x0=None, k_max=100, tol=1e-8)
    print(cg.solve(A, b, output=True))

if __name__ == '__main__':
    main_linear_systems()
