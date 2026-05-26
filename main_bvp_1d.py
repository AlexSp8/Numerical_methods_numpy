
import numpy as np

from ode import bvp_1d, plot

def f(x: float, y: np.ndarray[tuple[int]], dy: np.ndarray[tuple[int]],
    d2y: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    h, Ta = 0.01, 20.0
    y1 = y[0]
    return h*(Ta-y1) + d2y

def f_harm(x: float, y: np.ndarray[tuple[int]], dy: np.ndarray[tuple[int]],
    d2y: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    return 4.0*y[0] + d2y

def f_bratu(x: float, y: np.ndarray[tuple[int]], dy: np.ndarray[tuple[int]],
    d2y: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    return np.exp(y[0]) + d2y

# def f_system(x: float, y: np.ndarray[tuple[int]], dy: np.ndarray[tuple[int]],
#     d2y: np.ndarray[tuple[int]]) -> tuple[np.ndarray[tuple[int]]]:
#     return y[1], y[0]
def f_system(x: float, y: np.ndarray[tuple[int]], dy: np.ndarray[tuple[int]],
    d2y: np.ndarray[tuple[int]]) -> tuple[np.ndarray[tuple[int]]]:
    return y[1] + d2y[0], y[0] + d2y[1]

def main():

    print('BVP 1D')

    bc = {
        'left': {
            'type': 'dirichlet', 'value': [40.0], 'a_robin': None,
        },
        'right': {
            'type': 'dirichlet', 'value': [200.0], 'a_robin': None,
        }
    }
    xd = np.array([0.0, 10.0])
    x = bvp_1d.create_mesh(nnodes=11, xd=xd, p=1.0)
    y = bvp_1d.solve_bvp_1d(f, x, bc, neq=1)

    # #harmonic oscillator
    # bc = {
    # "left":  {"type": "dirichlet", "value": [1.0]},
    # "right": {"type": "dirichlet", "value": [0.0]}
    # }
    # x = bvp_1d.create_mesh(nnodes=11, xd=[0.0, np.pi])
    # y = bvp_1d.solve_bvp_1d(f_harm, x, bc, neq=1)

    # #Bratu
    # x = bvp_1d.create_mesh(nnodes=11, xd=[0.0, 1.0])
    # bc = {
    #     "left":  {"type": "dirichlet", "value": [0.0]},
    #     "right": {"type": "dirichlet", "value": [0.0]}
    # }
    # y = bvp_1d.solve_bvp_1d(f, x, bc, neq=1)

    for i in range(x.shape[0]):
        print(x[i], y[i])
    plot.plot_ode_evolution(x, y)

    # #Heat & Mass Transfer
    # x = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    # bc = {
    #     "left":  {"type": "dirichlet", "value": [100.0, 0.0]},
    #     "right": {"type": "neumann",   "value": [0.0, 1.0]}
    # }
    # y = bvp_1d.solve_bvp_1d(f_system, x, bc, neq=2)

    # neq = 2
    # nnodes = x.shape[0]

    # y_lists = [ [y[i*neq+j] for i in range(nnodes)] for j in range(neq) ]
    # print(f"{'Node':>5} | {'x':>10} | {'y1':>12} | {'y2':>12}")
    # print("-"*45)
    # for i in range(nnodes):
    #     y1 = y_lists[0][i]
    #     y2 = y_lists[1][i]
    #     print(f"{i:5d} | {x[i]:10.4f} | {y1:12.6e} | {y2:12.6e}")
    # plot.plot_ode_system_evolution(x, y_lists)

if __name__ == '__main__':
    main()
