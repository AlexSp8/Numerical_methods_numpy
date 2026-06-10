
import sys
from pathlib import Path

import numpy as np

current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent

if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from ode import bvp_1d_fd, plot

def bc_left_0(t: float, x: float, u_b: np.ndarray[tuple[int]],
    dudx_b: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    res = np.array([u_b[0] - 40.0])
    return res

def bc_right_0(t: float, x: float, u_b: np.ndarray[tuple[int]],
    dudx_b: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    res = np.array([u_b[-1] - 200.0])
    return res

def f(t: float, x: float, u: np.ndarray[tuple[int]], dudt: np.ndarray[tuple[int]],
    dudx: np.ndarray[tuple[int]], d2udx2: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    h, Ta = 0.01, 20.0
    y1 = u[0]
    return h*(Ta-y1) + d2udx2

# Harmonic oscillator
def bc_left_harm(t: float, x: float, u_b: np.ndarray[tuple[int]],
    dudx_b: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    res = np.array([u_b[0] - 1.0])
    return res

def bc_right_harm(t: float, x: float, u_b: np.ndarray[tuple[int]],
    dudx_b: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    res = np.array([u_b[-1] - 0.0])
    return res

def f_harm(t: float, x: float, u: np.ndarray[tuple[int]], dudt: np.ndarray[tuple[int]],
    dudx: np.ndarray[tuple[int]], d2udx2: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    return 4.0*u[0] + d2udx2

# Bratu
def bc_left_bratu(t: float, x: float, u_b: np.ndarray[tuple[int]],
    dudx_b: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    res = np.array([u_b[0] - 0.0])
    return res

def bc_right_bratu(t: float, x: float, u_b: np.ndarray[tuple[int]],
    dudx_b: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    res = np.array([u_b[-1] - 0.0])
    return res

def f_bratu(t: float, x: float, u: np.ndarray[tuple[int]], dudt: np.ndarray[tuple[int]],
    dudx: np.ndarray[tuple[int]], d2udx2: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    return np.exp(u[0]) + d2udx2

# Heat and Mass transfer coupled
def bc_left_hm(t: float, x: float, u_b: np.ndarray[tuple[int]],
    dudx_b: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    res = np.array([u_b[0] - 100.0, u_b[1] - 0.0])
    return res

def bc_right_hm(t: float, x: float, u_b: np.ndarray[tuple[int]],
    dudx_b: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    res = np.array([dudx_b[-2] - 0.0, dudx_b[-1]-10.0])
    return res

def f_hm(t: float, x: float, u: np.ndarray[tuple[int]], dudt: np.ndarray[tuple[int]],
    dudx: np.ndarray[tuple[int]], d2udx2: np.ndarray[tuple[int]]) -> tuple[np.ndarray[tuple[int]]]:
    return u[1] + d2udx2[0], u[0] + d2udx2[1]

# Coupled system
def f_coupled(t: float, x: float, u: np.ndarray[tuple[int]], dudt: np.ndarray[tuple[int]],
    dudx: np.ndarray[tuple[int]], d2udx2: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    res_0 = d2udx2[0] + u[0]*dudx[1] - 2.0 - np.pi*(x**2)*np.cos(np.pi*x)
    res_1 = d2udx2[1] + (dudx[0]**2)*u[1] - (4.0*x**2 - np.pi**2)*np.sin(np.pi*x)
    return np.array([res_0, res_1])

def left_bc_coupled(t: float, x: float, u_b: np.ndarray[tuple[int]], dudx_b: np.ndarray[tuple[int]]
    ) -> np.ndarray[tuple[int]]:
    return np.array([u_b[0] - 0.0, dudx_b[1] - np.pi*np.cos(u_b[0])])

def right_bc_coupled(t: float, x: float, u_b: np.ndarray[tuple[int]], dudx_b: np.ndarray[tuple[int]]
    ) -> np.ndarray[tuple[int]]:
    return np.array([u_b[0] - 1.0, u_b[1] - 0.0])

# Beam
def left_clamped_boundary(t: float, x: float, u_b: np.ndarray[tuple[int]],
    dudx_b: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    res_displacement = u_b[0] - 0.0  # u0 = 0
    res_slope        = dudx_b[0] - 0.0 # du0/dx = 0
    return np.array([res_displacement, res_slope])

def right_free_boundary(t: float, x: float, u_b: np.ndarray[tuple[int]],
    dudx_b: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    res_moment = u_b[1] - 0.0  # u1 = 0
    res_shear  = dudx_b[1] - 0.0 # du1/dx = 0
    return np.array([res_moment, res_shear])

def right_clamped_boundary(t: float, x: float, u_b: np.ndarray[tuple[int]],
    dudx_b: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    res_displacement = u_b[0] - 0.0   # u0 = 0
    res_slope        = dudx_b[0] - 0.0  # du0/dx = 0
    return np.array([res_displacement, res_slope])

def f_beam(t: float, x: float, u: np.ndarray[tuple[int]], dudt: np.ndarray[tuple[int]],
    dudx: np.ndarray[tuple[int]], d2udx2: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    # u0'' - u1 = 0
    res_0 = d2udx2[0] - u[1]

    # u1'' + q*u0 - g = 0
    q_x = 10.0 # Example stiffness coefficient
    g_x = -1.0 # Example downward uniform load
    res_1 = d2udx2[1] + q_x * u[0] - g_x

    return np.array([res_0, res_1])

def main():

    print('BVP 1D')

    bc = {'left': bc_left_0, 'right': bc_right_0}
    xd = np.array([0.0, 10.0])
    x = bvp_1d_fd.create_mesh(nnodes=11, xd=xd, p=1.0)
    u = bvp_1d_fd.solve_bvp_1d_fd(f, x, bc, neq=1)

    for i in range(x.shape[0]):
        print(x[i], u[i])
    plot.plot_ode_evolution(x, u)

    #harmonic oscillator
    bc = {'left': bc_left_harm, 'right': bc_right_harm}
    x = bvp_1d_fd.create_mesh(nnodes=11, xd=[0.0, np.pi])
    u = bvp_1d_fd.solve_bvp_1d_fd(f_harm, x, bc, neq=1)

    for i in range(x.shape[0]):
        print(x[i], u[i])
    plot.plot_ode_evolution(x, u)

    #Bratu
    bc = {'left': bc_left_bratu, 'right': bc_right_bratu}
    x = bvp_1d_fd.create_mesh(nnodes=11, xd=[0.0, 1.0])
    u = bvp_1d_fd.solve_bvp_1d_fd(f, x, bc, neq=1)

    for i in range(x.shape[0]):
        print(x[i], u[i])
    plot.plot_ode_evolution(x, u)


    #Heat & Mass Transfer
    bc = {'left': bc_left_hm, 'right': bc_right_hm}
    x = bvp_1d_fd.create_mesh(nnodes=5, xd=[0.0, 1.0], p=1.0)
    u = bvp_1d_fd.solve_bvp_1d_fd(f_hm, x, bc, neq=2)

    neq, nnodes = 2, x.shape[0]

    y_lists = np.zeros((neq,nnodes))
    for i in range(neq):
        for j in range(nnodes):
            y_lists[i,j] = u[j*neq+i]

    print(f"{'Node':>5} | {'x':>10} | {'y1':>12} | {'y2':>12}")
    print("-"*45)
    for i in range(nnodes):
        y1 = y_lists[0,i]
        y2 = y_lists[1,i]
        print(f"{i:5d} | {x[i]:10.4f} | {y1:12.6e} | {y2:12.6e}")
    plot.plot_ode_system_evolution(x, y_lists)

    # Coupled system
    bc = {'left': left_bc_coupled, 'right': right_bc_coupled}
    xd = np.array([0.0, 1.0])
    x = bvp_1d_fd.create_mesh(nnodes=21, xd=xd, p=1.0)
    u = bvp_1d_fd.solve_bvp_1d_fd(f_coupled, x, bc, neq=2)

    neq, nnodes = 2, x.shape[0]

    y_lists = np.zeros((neq,nnodes))
    for i in range(neq):
        for j in range(nnodes):
            y_lists[i,j] = u[j*neq+i]

    print(f"{'Node':>5} | {'x':>10} | {'y1':>12} | {'y2':>12}")
    print("-"*45)
    for i in range(nnodes):
        y1 = y_lists[0,i]
        y2 = y_lists[1,i]
        print(f"{i:5d} | {x[i]:10.4f} | {y1:12.6e} | {y2:12.6e}")
    plot.plot_ode_system_evolution(x, y_lists)

    # Beam
    bc = {'left': left_clamped_boundary, 'right': right_clamped_boundary}
    x = bvp_1d_fd.create_mesh(nnodes=11, xd=[0.0, 1.0])
    u = bvp_1d_fd.solve_bvp_1d_fd(f_beam, x, bc, neq=2)

    neq, nnodes = 2, x.shape[0]

    y_lists = np.zeros((neq,nnodes))
    for i in range(neq):
        for j in range(nnodes):
            y_lists[i,j] = u[j*neq+i]

    print(f"{'Node':>5} | {'x':>10} | {'y1':>12} | {'y2':>12}")
    print("-"*45)
    for i in range(nnodes):
        y1 = y_lists[0,i]
        y2 = y_lists[1,i]
        print(f"{i:5d} | {x[i]:10.4f} | {y1:12.6e} | {y2:12.6e}")
    plot.plot_ode_system_evolution(x, y_lists)

if __name__ == '__main__':
    main()
