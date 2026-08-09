
import sys
from pathlib import Path

from typing import Callable
import numpy as np
from numpy.typing import NDArray

current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent

if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from optimization.unconstrained import direct_methods, gradient_methods
from differentiation import forward_fd as ffd

def g(x: NDArray[np.float64]) -> float:
    x, y = x[0], x[1]
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

def g_max(x: NDArray[np.float64]) -> float:
    x, y = x[0], x[1]
    return -(2*x*y + 2*x - x**2 - 2*(y**2))

def g_p(x: NDArray[np.float64],
    g: Callable[[NDArray[np.float64]], float],
    xp: NDArray[np.float64],
    d: NDArray[np.float64],
    p: float = 1e1) -> float:
    """Test function with Gauss Radial penalty constraints"""

    penalty = 0.0
    for i in range(xp.shape[0]):
        dist_sq = np.sum(((x - xp[i])/d[i])**2)
        penalty += np.exp(-0.5*dist_sq)
        # penalty -= np.exp(-0.5*dist_sq)

    return g(x) + p*penalty

def main():

    print('\nUnconstrained optimization')

    # f = g
    f = g_max

    # print('\nConstrained')
    # f = lambda x: g_p(x, g, xp=np.array([[2.0, 1.0]]), d=np.array([0.1, 0.1]), mode='max', p=1e1)
    # f = lambda x: g_p(x, g_max, xp=np.array([[2.0, 1.0]]), d=np.array([0.1, 0.1]), mode='max', p=1e1)

    print('\nPowell')
    # Starting guess
    x0 = np.array([0.0, 0.0])
    # Initial direction set: Standard unit vectors for 2D [[1, 0], [0, 1]]
    d_vectors = np.array([
        [1.0, 0.0],
        [0.0, 1.0]
    ])
    x = direct_methods.powell(f, x0, d_vectors)
    print(f"x = {x}, f = {f(x)}")

    print('\nSteepest Descent')
    df = ffd.df_h
    x = gradient_methods.steepest_descent(f, x0, df=df)
    print(f"x = {x}, f = {f(x)}")

    print('\nConjugate Gradient')
    df = ffd.df_h
    x = gradient_methods.conjugate_gradient(f, x0, df=df)
    print(f"x = {x}, f = {f(x)}")

    print('\nNewton')
    df = ffd.df_h
    x = gradient_methods.newton(f, x0, df=df)
    print(f"x = {x}, f = {f(x)}")

    print('\nMarquardt')
    df = ffd.df_h
    x = gradient_methods.marquardt(f, x0, df=df)
    print(f"x = {x}, f = {f(x)}")

    print('\nBFGS')
    df = ffd.df_h
    x = gradient_methods.bfgs(f, x0, df=df, output=True)
    print(f"x = {x}, f = {f(x)}")

    print('\nL-BFGS')
    df = ffd.df_h
    x = gradient_methods.l_bfgs(f, x0, df=df, output=True)
    print(f"x = {x}, f = {f(x)}")

if __name__ == '__main__':
    main()
