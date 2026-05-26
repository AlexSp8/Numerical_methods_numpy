
import numpy as np

from ode import bvp_1d, plot, explicit_ode, implicit_ode
from pde import parabolic_1d

def f(x: float, y: np.ndarray[tuple[int]], dy: np.ndarray[tuple[int]],
    d2y: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    h, Ta = 0.01, 20.0
    y1 = y[0]
    return h*(Ta-y1) + d2y

def main():

    print('PDEs')

    xd = np.array([0.0, 10.0])
    x = bvp_1d.create_mesh(nnodes=11, xd=xd, p=1.0)

    neq = 1

    y0 = np.zeros((x.shape[0]*neq))
    y0[0] = 40.0
    y0[-1] = 200.0

    ode_solver = implicit_ode.implicit_solver

    t, y = parabolic_1d.solve_parabolic_1d(f, x, y0, neq=neq, ode_solver=ode_solver)

    for i in range(0, t.shape[0], 20):
        y_str = "[" + ", ".join([f"{val:.2e}" for val in y[i,:]]) + "]"
        print(f't = {t[i]:.2e}, y: {y_str}')
        plot.plot_ode_evolution(x, y[i,:])
    # plot.plot_ode_evolution(t, y)

if __name__ == '__main__':
    main()
