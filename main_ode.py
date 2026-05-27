
import numpy as np

from ode import plot
from ode import ivp

def f1(t: float, y: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    return np.array([-2*(t**3) + 12*(t**2) - 20*t + 8.5]) # Explicit Euler, RK2, RK3, RK4 test

def f2(t: float, y: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    y1 = y[0]
    return np.array([4*np.exp(0.8*t) -0.5*y1]) # Heun's, Midpoint RK4 test

def f_sys(t: float, y: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    y1, y2 = y[0], y[1]
    return np.array([-0.5*y1, 4.0-0.3*y2-0.1*y1]) # Explicit Euler, RK4 System

def f_adapt(t: float, y: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    y1 = y[0]
    return np.array([10*np.exp(-((t-2)**2)/(2*(0.075**2)))-0.6*y1]) # Adaptive test

def f_im(t: float, y: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    y1 = y[0]
    return np.array([-1000.0*y1 + 3000.0 - 2000.0*np.exp(-t)]) # Implicit test

def main():

    print('ODEs')

    # print('\nExplicit')
    # f, t0, tf, y0, method, h = f1, 0.0, 4.0, np.array([1.0]), 'euler', 0.5
    # f, t0, tf, y0, method, h = f_sys, 0.0, 2.0, np.array([4.0, 6.0]), 'euler', 0.5
    # f, t0, tf, y0, method, h = f2, 0.0, 4.0, np.array([2.0]), 'euler', 1.0

    # f, t0, tf, y0, method, h = f2, 0.0, 4.0, np.array([2.0]), 'heun', 1.0

    # f, t0, tf, y0, method, h = f2, 0.0, 4.0, np.array([2.0]), 'midpoint', 1.0

    # f, t0, tf, y0, method, h = f1, 0.0, 4.0, np.array([1.0]), 'ralston', 0.5

    # f, t0, tf, y0, method, h = f1, 0.0, 4.0, np.array([1.0]), 'rk3', 0.5

    # f, t0, tf, y0, method, h = f1, 0.0, 4.0, np.array([1.0]), 'rk4', 0.5
    # f, t0, tf, y0, method, h = f_sys, 0.0, 2.0, np.array([4.0, 6.0]), 'rk4', 0.5
    # f, t0, tf, y0, method, h = f2, 0.0, 4.0, np.array([2.0]), 'rk4', 1.0

    # f, t0, tf, y0, method, h = f1, 0.0, 4.0, np.array([1.0]), 'rk5', 0.5
    # f, t0, tf, y0, method, h = f2, 0.0, 4.0, np.array([2.0]), 'rk5', 1.0

    # ivp_solver = ivp.ExplicitIVPSolver(f, t0, tf, y0, method, h=h)


    # print('\nExplicit Adaptive')
    # f, t0, tf, y0, method = f_adapt, 0.0, 4.0, np.array([2.0]), 'fehlberg'
    # h_min, h_max, atol, rtol = 1e-6, 2.0, 1e-8, 1e-6

    # f, t0, tf, y0, method = f_adapt, 0.0, 4.0, np.array([2.0]), 'cash-karp'
    # h_min, h_max, atol, rtol = 1e-6, 2.0, 1e-8, 1e-6

    # f, t0, tf, y0, method = f_adapt, 0.0, 4.0, np.array([2.0]), 'dormand-prince'
    # h_min, h_max, atol, rtol = 1e-6, 2.0, 1e-8, 1e-6

    # ivp_solver = ivp.ExplicitAdaptiveIVPSolver(f, t0, tf, y0, method,
    #     h_min=h_min, h_max=h_max, atol=atol, rtol=rtol)


    # print('\nImplicit')
    # f, t0, tf, y0, method, h = f_im, 0.0, 5.0, np.array([0.0]), 'euler', 0.05
    # f, t0, tf, y0, method, h = f_im, 0.0, 5.0, np.array([0.0]), 'midpoint', 0.01
    # f, t0, tf, y0, method, h = f_im, 0.0, 5.0, np.array([0.0]), 'crank-nicolson', 0.01
    # f, t0, tf, y0, method, h = f_im, 0.0, 5.0, np.array([0.0]), 'gauss-legendre-2', 0.05
    # f, t0, tf, y0, method, h = f_im, 0.0, 5.0, np.array([0.0]), 'gauss-legendre-3', 0.05
    # f, t0, tf, y0, method, h = f_im, 0.0, 5.0, np.array([0.0]), 'radau-iia-2', 0.05
    # f, t0, tf, y0, method, h = f_im, 0.0, 5.0, np.array([0.0]), 'radau-iia-3', 0.05
    # f, t0, tf, y0, method, h = f_im, 0.0, 5.0, np.array([0.0]), 'lobatto-iiic-2', 0.05
    # f, t0, tf, y0, method, h = f_im, 0.0, 5.0, np.array([0.0]), 'lobatto-iiia-3', 0.05

    # ivp_solver = ivp.ImplicitIVPSolver(f, t0, tf, y0, method, h=h)


    # print('\nImplicit Adaptive')
    # f, t0, tf, y0, method = f_im, 0.0, 5.0, np.array([0.0]), 'sdirk-3a'
    # h_min, h_max, atol, rtol = 1e-3, 2.0, 1e-6, 1e-3

    # f, t0, tf, y0, method = f_im, 0.0, 5.0, np.array([0.0]), 'cash-sdirk-4'
    # h_min, h_max, atol, rtol = 1e-4, 2.0, 1e-8, 1e-4

    # f, t0, tf, y0, method = f_im, 0.0, 5.0, np.array([0.0]), 'radau-iia-5-adaptive'
    # h_min, h_max, atol, rtol = 1e-3, 2.0, 1e-6, 1e-2

    # ivp_solver = ivp.ImplicitAdaptiveIVPSolver(f, t0, tf, y0, method,
    #     h_min=h_min, h_max=h_max, atol=atol, rtol=rtol)


    # print('\nAdams-Bashforth')
    # f, t0, tf, y0, method, h = f2, 0.0, 4.0, np.array([2.0]), 'rk4', 0.5
    # ivp_solver = ivp.ExplicitMultistepIVPSolver(f, t0, tf, y0, method, h=h, order=5)


    # print('\nAdams Moulton')
    # f, t0, tf, y0, method, h = f2, 0.0, 4.0, np.array([2.0]), 'rk4', 0.5
    # ivp_solver = ivp.ImplicitMultistepIVPSolver(f, t0, tf, y0, method, h=h, order=5)


    print('\nABM PECE')
    f, t0, tf, y0, method, h = f2, 0.0, 4.0, np.array([2.0]), 'rk4', 0.5
    ivp_solver = ivp.PredictorCorrectorIVPSolver(f, t0, tf, y0, method, h=h, order=5)


    t, y = ivp_solver.solve()

    for i in range(t.shape[0]):
        print(t[i], y[i])
    plot.plot_ode_evolution(t, y)

if __name__ == '__main__':
    main()
