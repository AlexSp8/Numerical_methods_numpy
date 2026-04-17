
from typing import List
import math
import numpy as np
import numpy.typing as npt

from optimization.unconstrained_1D import bracketing, open_methods, hybrid
from optimization.unconstrained_multi import direct_methods, gradient_methods

def f(x: float) -> float:
    """Test algebraic function."""
    return 2*math.sin(x)-(x**2)/10

def f_2D(x: npt.NDArray[np.float64]) -> float:
    x, y = x[0], x[1]
    # return (x**2 + y - 11)**2 + (x + y**2 - 7)**2
    return 2*x*y + 2*x - x**2 - 2*(y**2)

def main():

    # print('1D unconstrained optimization')

    # print('\nBracketing')
    # print('\nGolden-section search')
    # print(bracketing.golden_section_search(f, a=0, b=4, mode='max'))
    # print('\nMulti-Golden-section search')
    # print(bracketing.multi_bracketing(f, a=0, b=10, n=2, method='golden-section', mode='max'))

    # print('\nParabolic Interpolation')
    # print(open_methods.parabolic_interpolation(f, a=0, b=4, mode='max'))
    # print('\nMulti-Parabolic Interpolation')
    # print(open_methods.multi_open(f, a=0, b=10, n=2, method='parabolic_interpolation', mode='max'))

    # print('\nSecant')
    # print(open_methods.secant(f, x0=1.0))
    # print('\nMulti-Secant')
    # print(open_methods.multi_open(f, a=0, b=10, n=2, method='secant'))

    # print('\nNewton')
    # print(open_methods.newton(f, x0=1.0))
    # print('\nMulti-Newton')
    # print(open_methods.multi_open(f, a=0, b=10, n=2, method='newton'))

    # print('\nBrent method')
    # print(hybrid.brent(f, a0=0, b0=4, mode='max', output=True))
    # print('\nMulti-Brent method')
    # print(hybrid.multi_hybrid(f, a=0, b=10, n=2, method='brent', mode='max', output=True))

    print('\nMulti unconstrained optimization')
    print('\nPowell')
    # Starting guess
    x0 = np.array([0.0, 0.0])
    # Initial direction set: Standard unit vectors for 2D [[1, 0], [0, 1]]
    d_vectors = np.array([
        [1.0, 0.0],
        [0.0, 1.0]
    ])
    # print(direct_methods.powell(f_2D, x0, d_vectors, mode='max'))

    # print('\nSteepest Descent')
    # print(gradient_methods.steepest_descent(f_2D, x0, mode='max'))
    # print('\nConjugate Gradient')
    # print(gradient_methods.conjugate_gradient(f_2D, x0, mode='max'))
    # print('\nNewton')
    print(gradient_methods.newton(f_2D, x0))
    print('\nMarquardt')
    print(gradient_methods.marquardt(f_2D, x0, mode='max'))
    print('\nBFGS')
    print(gradient_methods.bfgs(f_2D, x0, mode='max'))

if __name__ == '__main__':
    main()
