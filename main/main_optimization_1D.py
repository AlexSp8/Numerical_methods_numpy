
from typing import Callable
import sys
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent

if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from optimization.one_dimensional import bracketing, open_methods, hybrid
from differentiation import forward_fd as ffd

def f(x: float) -> float:
    """Test algebraic function."""
    return 2*np.sin(x)-(x**2)/10

def f_max(x: float) -> float:
    return -f(x)

def f_p(x: float, f: Callable[[float], float],
    xp: NDArray[np.float64],
    d: NDArray[np.float64],
    p: float = 1e1) -> float:
    """Test function with Gauss Radial penalty constraints"""

    penalty = 0
    for i in range(len(xp)):
        penalty += np.exp(-0.5*((x - xp[i])/d[i])**2)
        # penalty -= np.exp(-0.5*((x - xp[i])/d[i])**2)

    return f(x) + p*penalty

def main():

    print('1D optimization')

    g = f
    # g = f_max

    # print('\nConstrained')
    # g = lambda x: f_p(x, f, xp=np.array([0.5, 1.5]), d=np.array([0.1, 0.1]), p=1e1)
    # g = lambda x: f_p(x, f_max, xp=np.array([0.5, 1.5]), d=np.array([0.1, 0.1]), p=1e1)

    print('\nBracketing')
    print('\nGolden-section search')
    x = bracketing.golden_section_search(g, a=0, b=4)
    print(f"x = {x}, f = {f(x)}")

    print('\nMulti-Golden-section search')
    xi = bracketing.multi_bracketing(g, a=0, b=10, n=2, method='golden-section')
    for x in xi:
        print(f"x = {x}, f = {f(x)}")

    print('\nParabolic Interpolation')
    x = open_methods.parabolic_interpolation(g, a=0, b=4)
    print(f"x = {x}, f = {f(x)}")
    print('\nMulti-Parabolic Interpolation')
    xi = bracketing.multi_bracketing(g, a=0, b=10, n=2, method='parabolic_interpolation')
    for x in xi:
        print(f"x = {x}, f = {f(x)}")

    print('\nSecant')
    df = ffd.df_h
    x = open_methods.secant(g, x0=1.0, df=df)
    print(f"x = {x}, f = {f(x)}")
    print('\nMulti-Secant')
    xi = bracketing.multi_bracketing(g, a=0, b=10, n=2, method='secant', df=df)
    for x in xi:
        print(f"x = {x}, f = {f(x)}")

    print('\nNewton')
    df = ffd.df_h
    d2f = ffd.d2f_h
    x = open_methods.newton(g, x0=1.0, df=df, d2f=d2f)
    print(f"x = {x}, f = {f(x)}")
    print('\nMulti-Newton')
    xi = bracketing.multi_bracketing(g, a=0, b=10, n=2, method='newton', df=df, d2f=d2f)
    for x in xi:
        print(f"x = {x}, f = {f(x)}")

    print('\nBrent method')
    x = hybrid.brent(g, a0=0, b0=4, output=True)
    print(f"x = {x}, f = {f(x)}")
    print('\nMulti-Brent method')
    xi = bracketing.multi_bracketing(g, a=0, b=10, n=2, method='brent')
    for x in xi:
        print(f"x = {x}, f = {f(x)}")

if __name__ == '__main__':
    main()
