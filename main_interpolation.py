
from pathlib import Path

import numpy as np
import numpy.typing as npt

from utilities import indexing
from curve_fitting import interpolation

# Path to the Datasets directory
DATASETS_DIR = Path(__file__).parent / 'curve_fitting/data'

def f(xi: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    yi = np.log(xi)
    return yi

def main():

    print('\nInterpolation')

    print('\nVandermonde')
    xi = np.array([-2.0, 0.0, 1.0])
    yi = np.array([-27.0, -1.0, 0.0])
    xi_sorted, yi_sorted = indexing.bubble_sort(xi,yi)
    print(interpolation.vandermonde(xi_sorted, yi_sorted))

    print('\nNewton-Gregory')
    file_path = DATASETS_DIR / 'newton_interpolation_data.txt'
    data = np.loadtxt(file_path, delimiter=None, skiprows=0, usecols=(0,1))
    xi = data[:,0]
    # yi = data[:,-1]
    yi = f(xi)
    xi_sorted, yi_sorted = indexing.bubble_sort(xi,yi)
    xp = np.array([2.0])
    print(interpolation.newton_gregory_polynomials(xp,xi_sorted, yi_sorted))

    print('\nLagrange')
    file_path = DATASETS_DIR / 'newton_interpolation_data.txt'
    data = np.loadtxt(file_path, delimiter=None, skiprows=0, usecols=(0,1))
    xi = data[:,0]
    # yi = data[:,-1]
    yi = f(xi)
    xi_sorted, yi_sorted = indexing.bubble_sort(xi,yi)
    xp = np.array([2.0])
    print(interpolation.lagrange_polynomial(xp,xi_sorted,yi_sorted,m=3))

    print('\nLinear Splines')
    file_path = DATASETS_DIR / 'linear_spline_interpolation_data.txt'
    data = np.loadtxt(file_path, delimiter=None, skiprows=0, usecols=(0,1))
    xi = data[:,0]
    yi = data[:,-1]
    xi_sorted, yi_sorted = indexing.bubble_sort(xi,yi)
    xp = np.array([5.0])
    print(interpolation.linear_splines(xp, xi_sorted, yi_sorted))
    print('Quadratic Splines')
    print(interpolation.quadratic_splines(xp, xi_sorted, yi_sorted))
    print('Cubic Splines')
    print(interpolation.cubic_splines(xp, xi_sorted, yi_sorted))

if __name__ == '__main__':
    main()
