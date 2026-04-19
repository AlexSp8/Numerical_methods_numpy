
from typing import Callable, List
import numpy as np
import numpy.typing as npt

from optimization.one_dimensional import bracketing, open_methods, hybrid

def f(x: float) -> float:
    """Test algebraic function."""
    return 2*np.sin(x)-(x**2)/10

def f_p(x: float, f: Callable[[float], float],
    xp: npt.NDarray[np.float64], d: npt.NDarray[np.float64],
    mode: str = 'min', p: float = 1e1) -> float:
    """Test function with Gauss Radial penalty constraints"""

    if mode == 'min':
        s = +1.0
    else:
        s = -1.0

    penalty = 0
    for i in range(len(xp)):
        penalty += s*np.exp(-0.5*((x - xp[i])/d[i])**2)

    return f(x) + p*penalty

def main():

    print('1D optimization')

    g = f

    # print('\nConstrained')
    # g = lambda x: f_p(x, f, xp=np.array([0.5, 1.5]), d=np.array([0.1, 0.1]), mode='max', p=1e1)

    print('\nBracketing')
    print('\nGolden-section search')
    x = bracketing.golden_section_search(g, a=0, b=4, mode='max')
    print(f"x = {x}, f = {f(x)}")

    print('\nMulti-Golden-section search')
    xi = bracketing.multi_bracketing(g, a=0, b=10, n=2, method='golden-section', mode='max')
    for x in xi:
        print(f"x = {x}, f = {f(x)}")

    print('\nParabolic Interpolation')
    x = open_methods.parabolic_interpolation(g, a=0, b=4, mode='max')
    print(f"x = {x}, f = {f(x)}")
    print('\nMulti-Parabolic Interpolation')
    xi = open_methods.multi_open(g, a=0, b=10, n=2, method='parabolic_interpolation', mode='max')
    for x in xi:
        print(f"x = {x}, f = {f(x)}")

    print('\nSecant')
    x = open_methods.secant(g, x0=1.0)
    print(f"x = {x}, f = {f(x)}")
    print('\nMulti-Secant')
    xi = open_methods.multi_open(g, a=0, b=10, n=2, method='secant')
    for x in xi:
        print(f"x = {x}, f = {f(x)}")

    print('\nNewton')
    x = open_methods.newton(g, x0=1.0)
    print(f"x = {x}, f = {f(x)}")
    print('\nMulti-Newton')
    xi = open_methods.multi_open(g, a=0, b=10, n=2, method='newton')
    for x in xi:
        print(f"x = {x}, f = {f(x)}")

    print('\nBrent method')
    x = hybrid.brent(g, a0=0, b0=4, mode='max', output=True)
    print(f"x = {x}, f = {f(x)}")
    print('\nMulti-Brent method')
    xi = hybrid.multi_hybrid(g, a=0, b=10, n=2, method='brent', mode='max', output=True)
    for x in xi:
        print(f"x = {x}, f = {f(x)}")

if __name__ == '__main__':
    main()
