
import sys
from pathlib import Path

import numpy as np

current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent.parent

if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from linear_systems import direct_solver
from ode import plot
from pde import bvp_setup, mesh_discretization, boundary_conditions

# Time-dependent BCs
def bc_x0(t: float, x: float, u: np.ndarray[tuple[int]],
    grad_u: np.ndarray[tuple[int, int]]) -> np.ndarray[tuple[int]]:
    # Exact value at x=0 is sin(t)*(0) = 0
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - 0.0
    return res

def bc_xf(t: float, x: float, u: np.ndarray[tuple[int]],
    grad_u: np.ndarray[tuple[int, int]]) -> np.ndarray[tuple[int]]:
    # Exact value at x=2 is sin(t)*(4 + 2) = 6*sin(t)
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - 6.0*np.sin(t)
    return res

def initial_condition(t0: float, x: np.ndarray[tuple[int, int]], neq: int
    ) -> np.ndarray[tuple[int]]:
    nnodes = x.shape[0]
    u0 = np.zeros(nnodes*neq)
    return u0

def f_res(t: float, x: np.ndarray[tuple[int]], u: np.ndarray[tuple[int]],
    dudt: np.ndarray[tuple[int]], grad_u: np.ndarray[tuple[int, int]],
    hess_u: np.ndarray[tuple[int, int, int]]) -> np.ndarray[tuple[int]]:
    neq = u.shape[0]
    res = np.zeros(neq)
    d2udx2 = hess_u[0,0,0]
    dudx = grad_u[0,0]
    D, v = 1.0, 2.0
    xm = x[0]
    S = np.cos(t)*(xm**2 + xm) + 4.0*xm*np.sin(t)
    res[0] = dudt[0] - (D*d2udx2 - v*dudx + S)
    return res

def u_analytical(t: float, x: np.ndarray[tuple[int, int]]) -> np.ndarray[tuple[int]]:
    xm = x[:,0]
    return np.sin(t)*(xm**2+xm)

def main():

    neq = 1

    x0 = np.array([0.0])
    xf = np.array([2.0])
    nodes_dim = np.array([31])
    
    mesh = mesh_discretization.MeshDiscretization(x0, xf, nodes_dim)

    p = np.array([1])
    mesh.set_rectangular_mesh(p, major_order='row')

    bc = [bc_x0, bc_xf]
    boundary = boundary_conditions.BoundaryConditions(neq, bc)

    bvp_solver = bvp_setup.BVPSetupFD(neq, mesh, boundary, f_res)
    
    ls_solver = direct_solver.LUSolver()
    bvp_solver.set_ls_solver(ls_solver)

    t0, tf, dt, dt_min, dt_max, atol, rtol = 0.0, 2.0, 1e-2, 1e-2, 0.5, 1e-2, 1e-1
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

    for ieq in range(neq):
        ui_num, ui_ex = u[ieq::neq].copy(), u_exact[ieq::neq].copy()
        plot.plot_ode_system_evolution(x, [ui_num, ui_ex])

if __name__ == '__main__':
    main()