
from pathlib import Path

import numpy as np
import numpy.typing as npt

from ode import explicit_ode, implicit_ode, plot

# Path to the Datasets directory
DATASETS_DIR = Path(__file__).parent / 'ode/data'

def f1(x: float, y: npt.NDArray) -> npt.NDArray:
    return np.array([-2*(x**3) + 12*(x**2) - 20*x + 8.5]) # Explicit Euler, RK2, RK3, RK4 test

def f2(x: float, y: npt.NDArray) -> npt.NDArray:
    y1 = y[0]
    return np.array([4*np.exp(0.8*x) -0.5*y1]) # Heun's, Midpoint RK4 test

def f_sys(x: float, y: npt.NDArray) -> npt.NDArray:
    y1, y2 = y[0], y[1]
    return np.array([-0.5*y1, 4.0-0.3*y2-0.1*y1]) # Explicit Euler, RK4 System

def f_adapt(x: float, y: npt.NDArray) -> npt.NDArray:
    y1 = y[0]
    return np.array([10*np.exp(-((x-2)**2)/(2*(0.075**2)))-0.6*y1]) # Adaptive test

def f_im(x: float, y: npt.NDArray) -> npt.NDArray:
    y1 = y[0]
    return np.array([-1000.0*y1 + 3000.0 - 2000.0*np.exp(-x)]) # Implicit test

def main():

    print('ODEs')

    print('\nExplicit Euler')
    y0 = np.array([1.0])
    x, y = explicit_ode.explicit_euler(f1, x0=0.0, xf=4.0, h=0.5, y0=y0)
    y0 = np.array([4.0, 6.0])
    x, y = explicit_ode.explicit_euler(f_sys, x0=0.0, xf=2.0, h=0.5, y0=y0)
    y0 = np.array([2.0])
    x, y = explicit_ode.explicit_euler(f2, x0=0.0, xf=4.0, h=1.0, y0=y0)

    print('\nHeun method')
    y0 = np.array([2.0])
    x, y = explicit_ode.heun_method(f2, x0=0.0, xf=4.0, h=1.0, y0=y0, k_max=1)

    print('\nMidpoint method')
    y0 = np.array([2.0])
    x, y = explicit_ode.midpoint_method(f2, x0=0.0, xf=4.0, h=1.0, y0=y0)

    print('\nRK2')
    y0 = np.array([1.0])
    x, y = explicit_ode.rk2(f1, x0=0.0, xf=4.0, h=0.5, y0=y0, a2=0.5)

    print('\nRK3')
    y0 = np.array([2.0])
    x, y = explicit_ode.rk3(f1, x0=0.0, xf=4.0, h=0.5, y0=y0, p1=0.5, p2=1.0)

    print('\nRK4')
    y0 = np.array([1.0])
    x, y = explicit_ode.rk4(f1, x0=0.0, xf=4.0, h=0.5, y0=y0)
    y0 = np.array([2.0])
    x, y = explicit_ode.rk4(f2, x0=0.0, xf=4.0, h=1.0, y0=y0)
    y0 = np.array([4.0, 6.0])
    x, y = explicit_ode.rk4(f_sys, x0=0.0, xf=2.0, h=0.5, y0=y0)

    print('\nRK5')
    y0 = np.array([1.0])
    x, y = explicit_ode.rk5(f1, x0=0.0, xf=4.0, h=0.5, y0=y0)
    y0 = np.array([2.0])
    x, y = explicit_ode.rk5(f2, x0=0.0, xf=4.0, h=1.0, y0=y0)

    print('\nRK45 Fehlberg')
    y0 = np.array([0.5])
    x, y = explicit_ode.rk45_fehlberg(f_adapt, x0=0.0, xf=4.0, y0=y0, h_min=1e-6, h_max=2.0)

    print('\nRK45 Cash-Karp')
    y0 = np.array([0.5])
    x, y = explicit_ode.rk45_cash_karp(f_adapt, x0=0.0, xf=4.0, y0=y0, h_min=1e-6, h_max=2.0)

    print('\nRK45 Dormand-Prince')
    y0 = np.array([0.5])
    x, y = explicit_ode.rk45_dormand_prince(f_adapt, x0=0.0, xf=4.0, y0=y0, h_min=1e-6, h_max=2.0)

    # print('\nImplicit')
    # x, y = implicit_ode.implicit_solver(f_im, x0=0.0, xf=0.4, h=0.05, y0=[0.0],
    #                                    method='euler')

    print('\nAdams-Bashforth')
    y0 = np.array([2.0])
    x, y = explicit_ode.adams_bashforth(f2, x0=0.0, xf=4.0, h=1.0, y0=y0, order=4)

    # print('\nAdams-Moulton')
    # x, y = implicit_ode.adams_moulton(f2, x0=0.0, xf=4.0, h=1.0, y0=[2.0],
    #                                    order=4)

    print('\nABM PECE')
    y0 = np.array([2.0])
    x, y = explicit_ode.abm_pece(f2, x0=0.0, xf=4.0, h=1.0, y0=y0, order=4)

    for i in range(len(x)):
        print(x[i], y[i])
    plot.plot_ode_evolution(x, y)

if __name__ == '__main__':
    main()
