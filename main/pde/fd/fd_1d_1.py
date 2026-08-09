
import sys
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent.parent.parent

if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from linear_systems import direct_solver
from ode import plot
from pde import bvp_setup, mesh_discretization, boundary_conditions

ALPHA = 1.0
BETA = 0.0

# Linear, Dirichlet-Neumann
def bc_x0(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    grad_u: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - ALPHA
    return res

def bc_xf(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    grad_u: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)
    dudx = grad_u[0,:]
    res[0] = dudx[0] - BETA
    return res

def initial_condition(t0: float, x: NDArray[np.float64], neq: int
    ) -> NDArray[np.float64]:
    nnodes = x.shape[0]
    u0 = np.zeros(nnodes*neq)
    return u0

def f_res(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    dudt: NDArray[np.float64], grad_u: NDArray[np.float64],
    hess_u: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)
    d2udx2 = hess_u[0,0,:]
    dudx = grad_u[0,:]
    res[0] = d2udx2[0] - 3.0*dudx[0] + 2.0*u[0] - 2*x[0]
    return res

def u_analytical(t: float, x: NDArray[np.float64], neq: int) -> NDArray[np.float64]:
    xm = x[:,0]
    c2 = (BETA-1.0-np.exp(1)*(ALPHA-1.5))/(2.0*(np.exp(1))**2-np.exp(1))
    c1 = ALPHA-1.5-c2
    nnodes = xm.shape[0]
    u = np.zeros(nnodes*neq)
    u = c1*np.exp(xm) + c2*np.exp(2.0*xm)+xm+1.5
    return u

def main():

    neq = 1

    x0 = np.array([0.0])
    xf = np.array([1.0])
    mesh = mesh_discretization.MeshDiscretizationFD(x0=x0, xf=xf)
    nodes_dim = np.array([21])
    p = np.array([0.7])
    mesh.set_rectangular_mesh_fd(nodes_dim=nodes_dim, p=p, major_order='row')

    f_bc = [bc_x0, bc_xf]
    bc_type = [['dirichlet'], ['neumann']]
    boundary = boundary_conditions.BoundaryConditions(f_bc, bc_type)

    bvp_solver = bvp_setup.BVPSetupFD(neq, mesh, boundary, f_res)

    ls_solver = direct_solver.LUSolver()
    bvp_solver.set_ls_solver(ls_solver)

    t0, tf, dt, dt_min, dt_max, atol, rtol = 0.0, 1.0, 1e-1, 1e-1, 1.0, 1e-3, 1e-2
    dt = tf

    u0 = initial_condition(t0, mesh.x_mesh, neq)

    bvp_solver.set_nr_solver(u0, k_max=100, tol=1e-8, r=1.0)

    bvp_solver.set_time_int(t0, tf, dt, dt_min, dt_max, atol, rtol, u0)

    bvp_solver.set_fd_problem(theta=1.0)

    u = bvp_solver.solve(dtw=tf)

    x = mesh.x_mesh

    u_exact = u_analytical(tf, x, neq)
    err = np.max(np.abs(u_exact-u))
    print(f'err = {err:.4e}')

    for ieq in range(neq):
        ui_num, ui_ex = u[ieq::neq].copy(), u_exact[ieq::neq].copy()
        plot.plot_ode_system_evolution(x, [ui_num, ui_ex])

if __name__ == '__main__':
    main()
