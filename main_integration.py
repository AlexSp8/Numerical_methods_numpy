
from pathlib import Path

import numpy as np
import numpy.typing as npt

from integration import data_integration, function_integration
from utilities import io_utils

# Path to the Datasets directory
DATASETS_DIR = Path(__file__).parent / 'integration/data'

def f(x: float) -> float:
    return 0.2 + 25*x - 200*(x**2) + 675*(x**3) - 900*(x**4) + 400*(x**5)
    # return math.sin(x)

def f_2D(xi: npt.NDArray) -> float:
    x, y = xi[0], xi[1]
    return 2*x*y + 2*x - x**2 -2*(y**2) + 72

def main():

    print('Integration')

    print('\nTrapezoidal')
    print(function_integration.trapezoidal(f, a=0.0, b=0.8, n=10))

    print('\nSimpson 1/3')
    print(function_integration.simpson_13(f, a=0.0, b=0.8, n=4))

    print('\nSimpson 3/8')
    print(function_integration.simpson_38(f, a=0.0, b=0.8))

    print('\nSimpson Mixed')
    print(function_integration.simpson_mixed(f, a=0.0, b=0.8, n=13))

    print('\nRomberg')
    print(function_integration.romberg(f, a=0.0, b=0.8))

    print('\nGauss-Legendre')
    # print(function_integration.gauss_legendre_coefficients(n=1))
    # print(function_integration.solve_gauss_legendre_system(n=1))
    print(function_integration.gauss_legendre_1D(f, a=0.0, b=0.8, n=3))

    print('\nTrapezoidal data')
    file_path = DATASETS_DIR / 'integration_data1.txt'
    data = np.loadtxt(file_path, delimiter=None, skiprows=0, usecols=(0,1))
    x = data[:,0]
    y = data[:,-1]
    print(data_integration.trapezoidal(x, y))

    print('\nSimpson 1/3 data')
    print(data_integration.simpson_13(x, y))

    print('\nSimpson mixed data')
    print(data_integration.simpson_mixed(x, y))

    print('\nMulti-dimensional data')
    x, y = np.array([0, 4, 8]), np.array([0, 3, 6])
    z = np.array([
        [72, 64, 24],
        [54, 70, 54],
        [0, 40, 48]
    ])
    print(data_integration.simpson_mixed_2D(x, y, z))

    print(data_integration.simpson_mixed_multi([x, y],z))

    print('\nMulti-dimensional Gauss')
    a, b = np.array([0.0, 0.0]), np.array([8.0, 6.0])
    n = np.array([2,2])
    print(function_integration.gauss_legendre_nD(f_2D, a, b, n))
    print(function_integration.gauss_legendre_nD_recursive(f_2D, a, b, n))
    print(function_integration.gauss_legendre_2D(f_2D, a, b, n))

if __name__ == '__main__':
    main()
