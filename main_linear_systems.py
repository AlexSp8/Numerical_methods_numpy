
import numpy as np

from utilities import matrix_operations
from linear_systems import direct_solvers, iterative_solvers

def main_linear_systems():

    A = np.array([
        [0.3, 0.52, 1.0],
        [0.5, 1.0, 1.9],
        [0.1, 0.3, 0.5]
    ])
    b = np.array([-0.01, 0.67, -0.44])

    print('Cramer rule 3x3: ', end='')
    print(direct_solvers.cramer(A, b))

    print('\nGauss Elimination: ', end='')
    print(direct_solvers.gauss_elimination(A, b))

    print('\nLU decomposition: ', end='')
    print(direct_solvers.lu_decomposition_solve(A, b))

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
    print(direct_solvers.thomas_tri_diagonal(l, d, u, b))

    print('\nCholesky decomposition:')
    A = np.array([
        [6.0, 15.0, 55.0],
        [15.0, 55.0, 225.0],
        [55.0, 225.0, 979.0]
    ])
    print(matrix_operations.cholesky_decomposition(A))

    # print('\n')
    # A = [
    #     [3.0, -0.1, -0.2],
    #     [0.1, 7.0, -0.3],
    #     [0.3, -0.2, 10.0]
    # ]
    # b = [7.85, -19.3, 71.4]

    # print(iterative_solvers.jacobi(A, b, x0=None, output=True), '\n')
    # print(iterative_solvers.sor_solver(A, b, x0=None, output=True), '\n')
    # print(iterative_solvers.steepest_descent(A, b, x0=None, output=True), '\n')
    # print(iterative_solvers.conjugate_gradient(A, b, x0=None, output=True))

if __name__ == '__main__':
    main_linear_systems()
