
from typing import Callable
import numpy as np
import numpy.typing as npt

from optimization.unconstrained import direct_methods, gradient_methods

def g(x: npt.NDArray[np.float64]) -> float:
    x, y = x[0], x[1]
    # return (x**2 + y - 11)**2 + (x + y**2 - 7)**2
    return 2*x*y + 2*x - x**2 - 2*(y**2)

def g_p(x: npt.NDArray[np.float64], g: Callable[[npt.NDArray[np.float64]], float],
    xp: npt.NDarray[np.float64], d: npt.NDarray[np.float64],
    mode: str = 'min', p: float = 1e1) -> float:
    """Test function with Gauss Radial penalty constraints"""

    if mode == 'min':
        s = +1.0
    else:
        s = -1.0

    penalty = 0.0
    for i in range(xp.shape[0]):
        dist_sq = np.sum(((x - xp[i])/d[i])**2)
        penalty += np.exp(-0.5*dist_sq)

    return g(x) + s*p*penalty

def main():

    print('\nUnconstrained optimization')

    f = g

    # print('\nConstrained')
    # f = lambda x: g_p(x, g, xp=np.array([[2.0, 1.0]]), d=np.array([0.1, 0.1]), mode='max', p=1e1)

    print('\nPowell')
    # Starting guess
    x0 = np.array([0.0, 0.0])
    # Initial direction set: Standard unit vectors for 2D [[1, 0], [0, 1]]
    d_vectors = np.array([
        [1.0, 0.0],
        [0.0, 1.0]
    ])
    x = direct_methods.powell(f, x0, d_vectors, mode='max')
    print(f"x = {x}, f = {g(x)}")

    print('\nSteepest Descent')
    x = gradient_methods.steepest_descent(f, x0, fd_type='ffd', mode='max')
    print(f"x = {x}, f = {g(x)}")

    print('\nConjugate Gradient')
    x = gradient_methods.conjugate_gradient(f, x0, fd_type='ffd', mode='max')
    print(f"x = {x}, f = {g(x)}")

    print('\nNewton')
    x = gradient_methods.newton(f, x0, fd_type='ffd')
    print(f"x = {x}, f = {g(x)}")

    print('\nMarquardt')
    x = gradient_methods.marquardt(f, x0, fd_type='ffd')
    print(f"x = {x}, f = {g(x)}")

    print('\nBFGS')
    x = gradient_methods.bfgs(f, x0, fd_type='ffd', mode='max')
    print(f"x = {x}, f = {g(x)}")

if __name__ == '__main__':
    main()
