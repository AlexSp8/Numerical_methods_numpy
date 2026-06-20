
import sys
from pathlib import Path

import numpy as np

current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent.parent

if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from linear_systems import direct_solver
from ode import plot
from pde import bvp_setup, boundary_conditions, mesh_discretization

# Problem 4
def bc_x0(t: float, x: float, u: np.ndarray[tuple[int]],
    grad_u: np.ndarray[tuple[int, int]]) -> np.ndarray[tuple[int]]:
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - 0.0     # Dirichlet
    res[1] = grad_u[0,1] - 0.0  # Neumann
    return res

def bc_xf(t: float, x: float, u: np.ndarray[tuple[int]],
    grad_u: np.ndarray[tuple[int, int]]) -> np.ndarray[tuple[int]]:
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = grad_u[0,0] - 0.0  # Neumann
    res[1] = u[1] - 0.0     # Dirichlet
    return res

def initial_condition(t0: float, x: np.ndarray[tuple[int, int]], neq: int
    ) -> np.ndarray[tuple[int]]:

    nnodes = x.shape[0]
    u0 = np.zeros(nnodes*neq)

    u0[0::neq] = np.sin(x[:,0])
    u0[1::neq] = np.cos(x[:,0])

    return u0

def f_res(t: float, x: np.ndarray[tuple[int]], u: np.ndarray[tuple[int]],
    dudt: np.ndarray[tuple[int]], grad_u: np.ndarray[tuple[int, int]],
    hess_u: np.ndarray[tuple[int, int, int]]) -> np.ndarray[tuple[int]]:

    neq = u.shape[0]
    res = np.zeros(neq)

    xm = x[0]

    ext = np.exp(-t)
    cos_x = np.cos(xm)
    sin_x = np.sin(xm)

    S0 = ext*(cos_x - sin_x)
    S1 = -ext*(cos_x + sin_x)
    
    d2udx2 = hess_u[0,0,:]

    res[0] = dudt[0] - (d2udx2[1] + S0)
    res[1] = dudt[1] - (-d2udx2[0] + S1)
    return res

def u_analytical(t: float, x: np.ndarray[tuple[int, int]]) -> np.ndarray[tuple[int]]:
    nnodes = x.shape[0]
    neq = 2
    u_exact = np.zeros(nnodes*neq)
    xm = x[:,0]
    u_exact[0::neq] = np.exp(-t)*np.sin(xm)
    u_exact[1::neq] = np.exp(-t)*np.cos(xm)
    return u_exact

def main():

    neq = 2

    x0 = np.array([0.0])
    xf = np.array([np.pi/2])
    nodes_dim = np.array([11])
    
    mesh = mesh_discretization.MeshDiscretization(x0, xf, nodes_dim)

    p = np.array([1])
    mesh.set_rectangular_mesh(p, major_order='row')

    bc = [bc_x0, bc_xf]
    boundary = boundary_conditions.BoundaryConditions(neq, bc)

    bvp_solver = bvp_setup.BVPSetupFD(neq, mesh, boundary, f_res)
    
    ls_solver = direct_solver.LUSolver()
    bvp_solver.set_ls_solver(ls_solver)

    t0, tf, dt, dt_min, dt_max, atol, rtol = 0.0, 0.5, 0.01, 1e-3, 0.1, 1e-3, 1e-2
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