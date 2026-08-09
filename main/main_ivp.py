
import sys
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent

if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from ode import plot, ivp

def f1(t: float, u: NDArray[np.float64]) -> NDArray[np.float64]:
    dudt = np.array([-2*(t**3) + 12*(t**2) - 20*t + 8.5])
    return dudt # Explicit Euler, RK2, RK3, RK4 test

def f2(t: float, u: NDArray[np.float64]) -> NDArray[np.float64]:
    u1 = u[0]
    dudt = np.array([4*np.exp(0.8*t) -0.5*u1])
    return dudt # Heun's, Midpoint RK4 test

def f_sys(t: float, u: NDArray[np.float64]) -> NDArray[np.float64]:
    u1, u2 = u[0], u[1]
    dudt = np.array([-0.5*u1, 4.0-0.3*u2-0.1*u1])
    return dudt # Explicit Euler, RK4 System

def f_adapt(t: float, u: NDArray[np.float64]) -> NDArray[np.float64]:
    u1 = u[0]
    dudt = np.array([10*np.exp(-((t-2)**2)/(2*(0.075**2)))-0.6*u1])
    return dudt # Adaptive test

def f_im(t: float, u: NDArray[np.float64]) -> NDArray[np.float64]:
    u1 = u[0]
    dudt = np.array([-1000.0*u1 + 3000.0 - 2000.0*np.exp(-t)])
    return dudt # Implicit test

def main():

    print('ODEs')

    print('\nExplicit')
    f, t0, tf, u0, method, h = f1, 0.0, 4.0, np.array([1.0]), 'euler', 0.5
    ivp_solver = ivp.ExplicitIVPSolver(f, t0, tf, u0, method, h=h)
    t, u = ivp_solver.solve()
    f, t0, tf, u0, method, h = f_sys, 0.0, 2.0, np.array([4.0, 6.0]), 'euler', 0.5
    ivp_solver = ivp.ExplicitIVPSolver(f, t0, tf, u0, method, h=h)
    t, u = ivp_solver.solve()
    f, t0, tf, u0, method, h = f2, 0.0, 4.0, np.array([2.0]), 'euler', 1.0
    ivp_solver = ivp.ExplicitIVPSolver(f, t0, tf, u0, method, h=h)
    t, u = ivp_solver.solve()

    f, t0, tf, u0, method, h = f2, 0.0, 4.0, np.array([2.0]), 'heun', 1.0
    ivp_solver = ivp.ExplicitIVPSolver(f, t0, tf, u0, method, h=h)
    t, u = ivp_solver.solve()

    f, t0, tf, u0, method, h = f2, 0.0, 4.0, np.array([2.0]), 'midpoint', 1.0
    ivp_solver = ivp.ExplicitIVPSolver(f, t0, tf, u0, method, h=h)
    t, u = ivp_solver.solve()

    f, t0, tf, u0, method, h = f1, 0.0, 4.0, np.array([1.0]), 'ralston', 0.5
    ivp_solver = ivp.ExplicitIVPSolver(f, t0, tf, u0, method, h=h)
    t, u = ivp_solver.solve()

    f, t0, tf, u0, method, h = f1, 0.0, 4.0, np.array([1.0]), 'rk3', 0.5
    ivp_solver = ivp.ExplicitIVPSolver(f, t0, tf, u0, method, h=h)
    t, u = ivp_solver.solve()

    f, t0, tf, u0, method, h = f1, 0.0, 4.0, np.array([1.0]), 'rk4', 0.5
    ivp_solver = ivp.ExplicitIVPSolver(f, t0, tf, u0, method, h=h)
    t, u = ivp_solver.solve()
    f, t0, tf, u0, method, h = f_sys, 0.0, 2.0, np.array([4.0, 6.0]), 'rk4', 0.5
    ivp_solver = ivp.ExplicitIVPSolver(f, t0, tf, u0, method, h=h)
    t, u = ivp_solver.solve()
    f, t0, tf, u0, method, h = f2, 0.0, 4.0, np.array([2.0]), 'rk4', 1.0
    ivp_solver = ivp.ExplicitIVPSolver(f, t0, tf, u0, method, h=h)
    t, u = ivp_solver.solve()

    f, t0, tf, u0, method, h = f1, 0.0, 4.0, np.array([1.0]), 'rk5', 0.5
    ivp_solver = ivp.ExplicitIVPSolver(f, t0, tf, u0, method, h=h)
    f, t0, tf, u0, method, h = f2, 0.0, 4.0, np.array([2.0]), 'rk5', 1.0
    ivp_solver = ivp.ExplicitIVPSolver(f, t0, tf, u0, method, h=h)
    t, u = ivp_solver.solve()


    print('\nExplicit Adaptive')
    f, t0, tf, u0, method = f_adapt, 0.0, 4.0, np.array([2.0]), 'fehlberg'
    h_min, h_max, atol, rtol = 1e-6, 2.0, 1e-8, 1e-6
    ivp_solver = ivp.ExplicitAdaptiveIVPSolver(f, t0, tf, u0, method,
        h_min=h_min, h_max=h_max, atol=atol, rtol=rtol)
    t, u = ivp_solver.solve()

    f, t0, tf, u0, method = f_adapt, 0.0, 4.0, np.array([2.0]), 'cash-karp'
    h_min, h_max, atol, rtol = 1e-6, 2.0, 1e-8, 1e-6
    ivp_solver = ivp.ExplicitAdaptiveIVPSolver(f, t0, tf, u0, method,
        h_min=h_min, h_max=h_max, atol=atol, rtol=rtol)
    t, u = ivp_solver.solve()

    f, t0, tf, u0, method = f_adapt, 0.0, 4.0, np.array([2.0]), 'dormand-prince'
    h_min, h_max, atol, rtol = 1e-6, 2.0, 1e-8, 1e-6
    ivp_solver = ivp.ExplicitAdaptiveIVPSolver(f, t0, tf, u0, method,
        h_min=h_min, h_max=h_max, atol=atol, rtol=rtol)
    t, u = ivp_solver.solve()


    print('\nImplicit')
    f, t0, tf, u0, method, h = f_im, 0.0, 5.0, np.array([0.0]), 'euler', 0.05
    ivp_solver = ivp.ImplicitIVPSolver(f, t0, tf, u0, method, h=h)
    t, u = ivp_solver.solve()
    f, t0, tf, u0, method, h = f_im, 0.0, 5.0, np.array([0.0]), 'midpoint', 0.01
    ivp_solver = ivp.ImplicitIVPSolver(f, t0, tf, u0, method, h=h)
    t, u = ivp_solver.solve()
    f, t0, tf, u0, method, h = f_im, 0.0, 5.0, np.array([0.0]), 'crank-nicolson', 0.01
    ivp_solver = ivp.ImplicitIVPSolver(f, t0, tf, u0, method, h=h)
    t, u = ivp_solver.solve()
    f, t0, tf, u0, method, h = f_im, 0.0, 5.0, np.array([0.0]), 'gauss-legendre-2', 0.05
    ivp_solver = ivp.ImplicitIVPSolver(f, t0, tf, u0, method, h=h)
    t, u = ivp_solver.solve()
    f, t0, tf, u0, method, h = f_im, 0.0, 5.0, np.array([0.0]), 'gauss-legendre-3', 0.05
    ivp_solver = ivp.ImplicitIVPSolver(f, t0, tf, u0, method, h=h)
    t, u = ivp_solver.solve()
    f, t0, tf, u0, method, h = f_im, 0.0, 5.0, np.array([0.0]), 'radau-iia-2', 0.05
    ivp_solver = ivp.ImplicitIVPSolver(f, t0, tf, u0, method, h=h)
    t, u = ivp_solver.solve()
    f, t0, tf, u0, method, h = f_im, 0.0, 5.0, np.array([0.0]), 'radau-iia-3', 0.05
    ivp_solver = ivp.ImplicitIVPSolver(f, t0, tf, u0, method, h=h)
    t, u = ivp_solver.solve()
    f, t0, tf, u0, method, h = f_im, 0.0, 5.0, np.array([0.0]), 'lobatto-iiic-2', 0.05
    ivp_solver = ivp.ImplicitIVPSolver(f, t0, tf, u0, method, h=h)
    t, u = ivp_solver.solve()
    f, t0, tf, u0, method, h = f_im, 0.0, 5.0, np.array([0.0]), 'lobatto-iiia-3', 0.05
    ivp_solver = ivp.ImplicitIVPSolver(f, t0, tf, u0, method, h=h)
    t, u = ivp_solver.solve()


    print('\nImplicit Adaptive')
    f, t0, tf, u0, method = f_im, 0.0, 5.0, np.array([0.0]), 'sdirk-3a'
    h_min, h_max, atol, rtol = 1e-3, 2.0, 1e-6, 1e-3
    ivp_solver = ivp.ImplicitAdaptiveIVPSolver(f, t0, tf, u0, method,
        h_min=h_min, h_max=h_max, atol=atol, rtol=rtol)
    t, u = ivp_solver.solve()

    f, t0, tf, u0, method = f_im, 0.0, 5.0, np.array([0.0]), 'cash-sdirk-4'
    h_min, h_max, atol, rtol = 1e-4, 2.0, 1e-8, 1e-4
    ivp_solver = ivp.ImplicitAdaptiveIVPSolver(f, t0, tf, u0, method,
        h_min=h_min, h_max=h_max, atol=atol, rtol=rtol)
    t, u = ivp_solver.solve()

    f, t0, tf, u0, method = f_im, 0.0, 5.0, np.array([0.0]), 'radau-iia-5-adaptive'
    h_min, h_max, atol, rtol = 1e-3, 2.0, 1e-6, 1e-2
    ivp_solver = ivp.ImplicitAdaptiveIVPSolver(f, t0, tf, u0, method,
        h_min=h_min, h_max=h_max, atol=atol, rtol=rtol)
    t, u = ivp_solver.solve()


    print('\nAdams-Bashforth')
    f, t0, tf, u0, method, h = f2, 0.0, 4.0, np.array([2.0]), 'rk4', 0.5
    ivp_solver = ivp.ExplicitMultistepIVPSolver(f, t0, tf, u0, method, h=h, order=5)
    t, u = ivp_solver.solve()


    print('\nAdams Moulton')
    f, t0, tf, u0, method, h = f2, 0.0, 4.0, np.array([2.0]), 'rk4', 0.5
    ivp_solver = ivp.ImplicitMultistepIVPSolver(f, t0, tf, u0, method, h=h, order=5)
    t, u = ivp_solver.solve()


    print('\nABM PECE')
    f, t0, tf, u0, method, h = f2, 0.0, 4.0, np.array([2.0]), 'rk4', 0.5
    ivp_solver = ivp.PredictorCorrectorIVPSolver(f, t0, tf, u0, method, h=h, order=5)
    t, u = ivp_solver.solve()

    for i in range(t.shape[0]):
        print(t[i], u[i])
    plot.plot_ode_evolution(t, u)

if __name__ == '__main__':
    main()
