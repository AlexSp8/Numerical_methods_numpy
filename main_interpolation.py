
from pathlib import Path
from typing import List
import math

from utilities import io_utils, indexing
from curve_fitting import interpolation

# Path to the Datasets directory
DATASETS_DIR = Path(__file__).parent / 'curve_fitting/data'

def f(xi: List[float]) -> List[float]:
    n = len(xi)
    yi = [math.log(xi[i]) for i in range(n)]
    return yi

def main():

    print('\nInterpolation')

    print('\nVandermonde')
    xi = [-2.0, 0.0, 1.0]
    yi = [-27.0, -1.0, 0.0]
    xi_sorted, yi_sorted = indexing.bubble_sort(xi,yi)
    print(interpolation.vandermonde(xi_sorted, yi_sorted))

    print('\nNewton-Gregory')
    file_path = DATASETS_DIR / 'newton_interpolation_data.txt'
    data = io_utils.read_file_float_data(file_path)
    xi = [row[0] for row in data]
    # yi = [row[1] for row in data]
    yi = f(xi)
    xi_sorted, yi_sorted = indexing.bubble_sort(xi,yi)
    xp = [2.0]
    print(interpolation.newton_gregory_polynomials(xp,xi_sorted, yi_sorted))

    print('\nLagrange')
    file_path = DATASETS_DIR / 'newton_interpolation_data.txt'
    data = io_utils.read_file_float_data(file_path)
    xi = [row[0] for row in data]
    # yi = [row[1] for row in data]
    yi = f(xi)
    xi_sorted, yi_sorted = indexing.bubble_sort(xi,yi)
    xp = [2.0]
    print(interpolation.lagrange_polynomial(xp,xi_sorted,yi_sorted,m=3))

    print('\nLinear Splines')
    file_path = DATASETS_DIR / 'linear_spline_interpolation_data.txt'
    data = io_utils.read_file_float_data(file_path)
    xi = [row[0] for row in data]
    yi = [row[1] for row in data]
    xi_sorted, yi_sorted = indexing.bubble_sort(xi,yi)
    xp = [5.0]
    print(interpolation.linear_splines(xp, xi_sorted, yi_sorted))
    print('Quadratic Splines')
    print(interpolation.quadratic_splines(xp, xi_sorted, yi_sorted))
    print('Cubic Splines')
    print(interpolation.cubic_splines(xp, xi_sorted, yi_sorted))

    print('\nMulti-dimensional Interpolation')
    xi = [[2.0, 9.0], [1.0, 6.0]]
    yi = [[60.0, 55.0], [57.5, 70.0]]
    xp = [[5.25, 4.8]]
    print('Multi-dimensional Lagrange')
    for p in xp:
        print(interpolation.multi_lagrange(p, xi, yi, m=1))
    print('Multi-dimensional Cubic Splines')
    for p in xp:
        print(interpolation.multi_cubic_spline(p, xi, yi))
    print('Multi-dimensional Interpolation')
    for p in xp:
        print(interpolation.multi_interpolation(
            p, xi, yi,method_1d=interpolation.lagrange_polynomial))
    for p in xp:
        print(interpolation.multi_interpolation(
            p, xi, yi,method_1d=interpolation.cubic_splines))

if __name__ == '__main__':
    main()
