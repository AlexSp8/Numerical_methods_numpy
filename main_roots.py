
import math

from roots import bracketing, open_methods, hybrid

def f(x: float) -> float:
    """Test algebraic function."""
    return math.exp(-x) - x
    # return math.sin(x)*math.exp(x)
    # return x**10 - 1.0

def main_roots():

    print('Roots')

    print('\nBisection')
    print('xr =', bracketing.bisection(f, a=0, b=10))
    print('xr =', bracketing.multi_bracketing(f, a=0, b=10, n=10, method='bisection'))

    print('\nFalse position')
    print('xr =', bracketing.multi_bracketing(f, a=0, b=10, n=10, method='false position'))

    print('\nFixed Point Iteration')
    print('xr =', open_methods.fixed_point(f, x0=0))

    print('\nNewton-Raphson')
    print('xr =', open_methods.newton_raphson(f, x0=0))

    print('\nRalston-Rabinowitz')
    print('xr =', open_methods.ralston_rabinowitz(f, x0=0))

    print('\nBrent')
    print('xr =', hybrid.brent(f, a=0, b=1))

if __name__ == '__main__':
    main_roots()
