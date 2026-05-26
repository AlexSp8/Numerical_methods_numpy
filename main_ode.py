
import numpy as np

from ode import explicit_ode, implicit_ode, plot

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

    print('\nExplicit Euler')
    y0 = np.array([1.0])
    t, y = explicit_ode.explicit_euler(f1, t0=0.0, tf=4.0, y0=y0, h=0.5)
    y0 = np.array([4.0, 6.0])
    t, y = explicit_ode.explicit_euler(f_sys, t0=0.0, tf=2.0, y0=y0, h=0.5)
    y0 = np.array([2.0])
    t, y = explicit_ode.explicit_euler(f2, t0=0.0, tf=4.0, y0=y0, h=1.0)

    print('\nHeun method')
    y0 = np.array([2.0])
    t, y = explicit_ode.heun(f2, t0=0.0, tf=4.0, y0=y0, h=1.0, k_max=1)

    print('\nMidpoint method')
    y0 = np.array([2.0])
    t, y = explicit_ode.midpoint(f2, t0=0.0, tf=4.0, y0=y0, h=1.0)

    print('\nRK2')
    y0 = np.array([1.0])
    t, y = explicit_ode.rk2(f1, t0=0.0, tf=4.0, y0=y0, h=0.5, a2=0.5)

    print('\nRK3')
    y0 = np.array([2.0])
    t, y = explicit_ode.rk3(f1, t0=0.0, tf=4.0, y0=y0, h=0.5, p1=0.5, p2=1.0)

    print('\nRK4')
    y0 = np.array([1.0])
    t, y = explicit_ode.rk4(f1, t0=0.0, tf=4.0, y0=y0, h=0.5)
    y0 = np.array([2.0])
    t, y = explicit_ode.rk4(f2, t0=0.0, tf=4.0, y0=y0, h=1.0)
    y0 = np.array([4.0, 6.0])
    t, y = explicit_ode.rk4(f_sys, t0=0.0, tf=2.0, y0=y0, h=0.5)

    print('\nRK5')
    y0 = np.array([1.0])
    t, y = explicit_ode.rk5(f1, t0=0.0, tf=4.0, y0=y0, h=0.5)
    y0 = np.array([2.0])
    t, y = explicit_ode.rk5(f2, t0=0.0, tf=4.0, y0=y0, h=1.0)

    print('\nRK45 Fehlberg')
    y0 = np.array([0.5])
    t, y = explicit_ode.rk45_fehlberg(f_adapt, t0=0.0, tf=4.0, y0=y0, h_min=1e-6, h_max=2.0)

    print('\nRK45 Cash-Karp')
    y0 = np.array([0.5])
    t, y = explicit_ode.rk45_cash_karp(f_adapt, t0=0.0, tf=4.0, y0=y0, h_min=1e-6, h_max=2.0)

    print('\nRK45 Dormand-Prince')
    y0 = np.array([0.5])
    t, y = explicit_ode.rk45_dormand_prince(f_adapt, t0=0.0, tf=4.0, y0=y0, h_min=1e-6, h_max=2.0)

    print('\nAdams-Bashforth')
    y0 = np.array([2.0])
    t, y = explicit_ode.adams_bashforth(f2, t0=0.0, tf=4.0, y0=y0, h=1.0, order=4)

    print('\nABM PECE')
    y0 = np.array([2.0])
    t, y = explicit_ode.abm_pece(f2, t0=0.0, tf=4.0, y0=y0, h=1.0, order=4)

    print('\nImplicit')
    y0 = np.array([0.0])
    t, y = implicit_ode.implicit_solver(f_im, t0=0.0, tf=5.0, y0=y0, h=0.005, method='euler')

    print('\nImplicit Adaptive')
    y0 = np.array([0.0])
    t, y = implicit_ode.implicit_adaptive_solver(f_im, t0=0.0, tf=5.0, y0=y0, method='euler')

    print('\nAdams-Moulton')
    y0 = np.array([2.0])
    t, y = implicit_ode.adams_moulton(f2, t0=0.0, tf=4.0, y0=y0, h=1.0, order=4)

    for i in range(t.shape[0]):
        print(t[i], y[i])
    plot.plot_ode_evolution(t, y)

if __name__ == '__main__':
    main()
