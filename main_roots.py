
import numpy as np

from roots import bracketing, open_methods, hybrid

def f(x: float) -> float:
    """Test algebraic function."""
    # return np.exp(-x) - x
    return np.sin(x)*np.exp(x)
    # return x**10 - 1.0

def main_roots():

    print('Roots')

    print('\nBracketing')

    print('\nBisection')
    print('xr =', bracketing.bisection(f, a=0, b=10))
    print('\nMulti-Bisection')
    print('xr =', bracketing.multi_bracketing(f, a=0, b=10, n=10, method='bisection'))

    print('\nFalse Position')
    print('xr =', bracketing.false_position(f, a=0, b=10))
    print('\nMulti-False position')
    print('xr =', bracketing.multi_bracketing(f, a=0, b=10, n=10, method='false_position'))

    print('\nOpen')

    print('\nFixed Point Iteration')
    print('xr =', open_methods.fixed_point(f, x0=1.0))
    # print('\nMulti-Fixed Point')
    # print('xr =', open_methods.multi_open(f, a=0, b=10, n=10, method='fixed_point'))

    print('\nSecant')
    print('xr =', open_methods.secant(f, x0=1.0, x1=2.0))
    print('\nMulti-Secant')
    print('xr =', open_methods.multi_open(f, a=0, b=10, n=10, method='secant'))

    print('\nIQI')
    print('xr =', open_methods.iqi(f, x0=1.0, x1=2.0, x2=3.0))
    print('\nMulti-IQI')
    print('xr =', open_methods.multi_open(f, a=0, b=10, n=10, method='iqi'))

    print('\nNewton-Raphson')
    print('xr =', open_methods.newton_raphson(f, x0=1.0))
    print('\nMulti-Newton-Raphson')
    print('xr =', open_methods.multi_open(f, a=0, b=10, n=10, method='newton_raphson'))

    print('\nRalston-Rabinowitz')
    print('xr =', open_methods.ralston_rabinowitz(f, x0=1.0))
    print('\nMulti-Ralston-Rabinowitz')
    print('xr =', open_methods.multi_open(f, a=0, b=10, n=10, method='ralston_rabinowitz'))

    print('\nHybrid')

    print('\nBrent')
    print('xr =', hybrid.brent(f, a0=0, b0=10))
    print('\nMulti-Brent')
    print('xr =', hybrid.multi_hybrid(f, a=0, b=10, n=10, method='brent'))

    print('\nRidders')
    print('xr =', hybrid.ridders(f, a=0, b=10))
    print('\nMulti-Ridders')
    print('xr =', hybrid.multi_hybrid(f, a=0, b=10, n=10, method='ridders'))

    print('\nChandrupatla')
    print('xr =', hybrid.chandrupatla(f, a0=0, b0=10))
    print('\nMulti-Chandrupatla')
    print('xr =', hybrid.multi_hybrid(f, a=0, b=10, n=10, method='chandrupatla'))

if __name__ == '__main__':
    main_roots()
