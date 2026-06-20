
import sys
from pathlib import Path

import numpy as np

current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent.parent

if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from linear_systems import direct_solver
from pde import bvp_setup, mesh_discretization, plot_pde, boundary_conditions

# Heat Transfer
def bc_x0(t: float, x: np.ndarray[tuple[int]], u: np.ndarray[tuple[int]],
    grad_u: np.ndarray[tuple[int, int]]) -> np.ndarray[tuple[int]]:
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - 0.0
    return res

def bc_xf(t: float, x: np.ndarray[tuple[int]], u: np.ndarray[tuple[int]],
    grad_u: np.ndarray[tuple[int, int]]) -> np.ndarray[tuple[int]]:
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - 0.0
    return res

def bc_y0(t: float, x: np.ndarray[tuple[int]], u: np.ndarray[tuple[int]],
    grad_u: np.ndarray[tuple[int, int]]) -> np.ndarray[tuple[int]]:
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - 0.0
    return res

def bc_yf(t: float, x: np.ndarray[tuple[int]], u: np.ndarray[tuple[int]],
    grad_u: np.ndarray[tuple[int, int]]) -> np.ndarray[tuple[int]]:
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - np.sin(np.pi*x[0])
    return res

def initial_condition(t0: float, x: np.ndarray[tuple[int, int]], neq: int
    ) -> np.ndarray[tuple[int]]:
    nnodes = x.shape[0]
    u0 = np.zeros(nnodes*neq)
    ieq = 0
    u0[ieq::neq] = 0.0 # x
    return u0

def f_res(t: float, x: np.ndarray[tuple[int]], u: np.ndarray[tuple[int]],
    dudt: np.ndarray[tuple[int]], grad_u: np.ndarray[tuple[int, int]],
    hess_u: np.ndarray[tuple[int, int, int]]) -> np.ndarray[tuple[int]]:
    neq = u.shape[0]
    res = np.zeros(neq)
    d2udx2, d2udy2 = hess_u[0,0,0], hess_u[1,1,0]
    res[0] = (d2udx2+ d2udy2)
    return res

def u_analytical(t: float, x: np.ndarray[tuple[int, int]]) -> np.ndarray[tuple[int]]:
    return np.sinh(np.pi*x[:,1])*np.sin(np.pi*x[:,0])/np.sinh(np.pi)

def main():
    
    neq = 1

    x0 = np.array([0.0, 0.0])
    xf = np.array([1.0, 1.0])
    nodes_dim = np.array([11, 11])

    mesh = mesh_discretization.MeshDiscretization(x0, xf, nodes_dim)

    p = np.array([1, 1])
    mesh.set_rectangular_mesh(p, major_order='row')

    # bc = [[bc_x0, bc_y0], [bc_xf, bc_yf]]
    bc = [bc_x0, bc_xf, bc_y0, bc_yf]
    boundary = boundary_conditions.BoundaryConditions(neq, bc)

    bvp_solver = bvp_setup.BVPSetupFD(neq, mesh, boundary, f_res)
    
    ls_solver = direct_solver.LUSolver()
    bvp_solver.set_ls_solver(ls_solver)

    t0, tf, dt, dt_min, dt_max, atol, rtol = 0.0, 20.0, 1e-1, 1e-1, 1.0, 1e-3, 1e-2
    dt = tf
    
    u0 = initial_condition(t0, mesh.x_mesh, neq)

    bvp_solver.set_nr_solver(u0, k_max=100, tol=1e-8, r=1.0)

    bvp_solver.set_time_int(t0, tf, dt, dt_min, dt_max, atol, rtol, u0)

    bvp_solver.set_fd_problem(theta=1.0)

    u = bvp_solver.solve(dtw=tf)

    x = mesh.x_mesh
    u_exact = u_analytical(tf, x)
    err = np.max(np.abs(u_exact-u))
    print(f'err = {err:.4e}')

    plot_pde.plot_contour(mesh, u, u_exact)

if __name__ == '__main__':
    main()
