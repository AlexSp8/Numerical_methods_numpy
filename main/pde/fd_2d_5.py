
import sys
from pathlib import Path

import numpy as np

current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent.parent

if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from linear_systems import direct_solver
from pde import bvp_setup, mesh_discretization, plot_pde, boundary_conditions

# System
def bc_x0(t: float, x: np.ndarray[tuple[int]], u: np.ndarray[tuple[int]],
    grad_u: np.ndarray[tuple[int, int]]) -> np.ndarray[tuple[int]]:
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - 0.0
    res[1] = u[1] - t*np.cos(np.pi*x[1])
    return res

def bc_xf(t: float, x: np.ndarray[tuple[int]], u: np.ndarray[tuple[int]],
    grad_u: np.ndarray[tuple[int, int]]) -> np.ndarray[tuple[int]]:
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - 0.0
    res[1] = u[1] - t*np.cos(np.pi*x[1])
    return res

def bc_y0(t: float, x: np.ndarray[tuple[int]], u: np.ndarray[tuple[int]],
    grad_u: np.ndarray[tuple[int, int]]) -> np.ndarray[tuple[int]]:
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - 0.0
    res[1] = u[1] - t*np.cos(2.0*np.pi*x[0])
    return res

def bc_yf(t: float, x: np.ndarray[tuple[int]], u: np.ndarray[tuple[int]],
    grad_u: np.ndarray[tuple[int, int]]) -> np.ndarray[tuple[int]]:
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - 0.0
    res[1] = u[1] - t*np.cos(2.0*np.pi*x[0])*(-1.0)
    return res

def initial_condition(t0: float, x: np.ndarray[tuple[int, int]], neq: int
    ) -> np.ndarray[tuple[int]]:
    nnodes = x.shape[0]
    u0 = np.zeros(nnodes*neq)
    return u0

def f_res(t: float, x: np.ndarray[tuple[int]], u: np.ndarray[tuple[int]],
    dudt: np.ndarray[tuple[int]], grad_u: np.ndarray[tuple[int, int]],
    hess_u: np.ndarray[tuple[int, int, int]]) -> np.ndarray[tuple[int]]:

    pi = np.pi
    pi2 = pi**2
    xm, ym = x[0], x[1]
    
    u_t     = dudt[0]
    d2udx2  = hess_u[0, 0, 0]
    d2udxdy = hess_u[0, 1, 0]
    
    v_t     = dudt[1]
    d2vdy2  = hess_u[1, 1, 1]
    d2vdxdy = hess_u[0, 1, 1]
    
    s1 = ((2.0*t + pi2*(t**2))*np.sin(pi*xm)*np.sin(2.0*pi*ym) 
          - 2.0*pi2*t*np.sin(2.0*pi*xm)*np.sin(pi*ym))
          
    s2 = ((1.0 + pi2*t)*np.cos(2.0*pi*xm)*np.cos(pi*ym) 
          - 2.0*pi2*(t**2)*np.cos(pi*xm)*np.cos(2.0*pi*ym))
    
    res_u = u_t - d2udx2 - d2vdxdy - s1
    res_v = v_t - d2vdy2 - d2udxdy - s2
    
    return np.array([res_u, res_v])

def u_analytical(t: float, x: np.ndarray[tuple[int, int]]) -> np.ndarray[tuple[int]]:

    neq = 2
    total_nodes = x.shape[0]
    u_exact = np.zeros(total_nodes*neq)

    xm = x[:,0]
    ym = x[:,1]
    
    u_exact[0::neq] = (t**2)*np.sin(np.pi*xm)*np.sin(2.0*np.pi*ym)
    u_exact[1::neq] = t*np.cos(2.0*np.pi*xm)*np.cos(np.pi*ym)

    return u_exact

def main():
    
    neq = 2

    x0 = np.array([0.0, 0.0])
    xf = np.array([1.0, 1.0])
    nodes_dim = np.array([21, 21])

    mesh = mesh_discretization.MeshDiscretization(x0, xf, nodes_dim)

    p = np.array([1, 1])
    mesh.set_rectangular_mesh(p, major_order='row')

    bc = [bc_x0, bc_xf, bc_y0, bc_yf]
    boundary = boundary_conditions.BoundaryConditions(neq, bc)

    bvp_solver = bvp_setup.BVPSetupFD(neq, mesh, boundary, f_res)
    
    # ls_solver = direct_solver.LUSolver()
    ls_solver = direct_solver.SparseSolver()
    bvp_solver.set_ls_solver(ls_solver)

    t0, tf, dt, dt_min, dt_max, atol, rtol = 0.0, 0.5, 1e-1, 5e-2, 1.0, 1e-2, 1e-1
    # dt = tf
    
    u0 = initial_condition(t0, mesh.x_mesh, neq)

    bvp_solver.set_nr_solver(u0, k_max=100, tol=1e-8, r=1.0)

    bvp_solver.set_time_int(t0, tf, dt, dt_min, dt_max, atol, rtol, u0)

    bvp_solver.set_fd_problem(theta=1.0)

    u = bvp_solver.solve(dtw=tf)
    
    x = mesh.x_mesh
    u_exact = u_analytical(tf, x)
    err = np.max(np.abs(u_exact-u))
    print(f'err = {err:.4e}')

    for i in range(neq):
        ue = u_exact[i::neq]
        un = u[i::neq]
        plot_pde.plot_contour(mesh, un, ue)

if __name__ == '__main__':
    main()
