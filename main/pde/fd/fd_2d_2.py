
import sys
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent.parent.parent

if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from linear_systems import direct_solver
from pde import bvp_setup, mesh_discretization, plot_pde, boundary_conditions

# Linear: Mixed derivatives, time-dependent BCs
def bc_x0(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    grad_u: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - 0.0
    return res

def bc_xf(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    grad_u: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)
    y = x[1]
    res[0] = u[0] - y*(t**2)
    return res

def bc_y0(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    grad_u: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - 0.0
    return res

def bc_yf(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    grad_u: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - x[0]*(t**2)
    return res

def initial_condition(t0: float, x: NDArray[np.float64], neq: int
    ) -> NDArray[np.float64]:
    nnodes = x.shape[0]
    u0 = np.zeros(nnodes*neq)
    u0[0::neq] = np.sin(np.pi*x[:,0])*np.sin(np.pi*x[:,1])
    return u0

def f_res(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    dudt: NDArray[np.float64], grad_u: NDArray[np.float64],
    hess_u: NDArray[np.float64]) -> NDArray[np.float64]:

    neq = u.shape[0]
    res = np.zeros(neq)

    beta = 0.5
    pi = np.pi

    xm, ym = x[0], x[1]

    d2udx2 = hess_u[0,0,:]
    d2udy2 = hess_u[1,1,:]
    d2udxdy = hess_u[0,1,:]

    sin_x, sin_y = np.sin(pi*xm), np.sin(pi*ym)
    cos_x, cos_y = np.cos(pi*xm), np.cos(pi*ym)

    s = (2.0*pi**2*sin_x*sin_y*np.cos(pi*t) - beta*pi**2*cos_x*cos_y*np.cos(pi*t)
         - pi*sin_x*sin_y*np.sin(pi*t) + 2.0*xm*ym*t - beta*t**2)

    res[0] = dudt[0] - (d2udx2[0] + beta*d2udxdy[0] + d2udy2[0]) - s

    return res

def u_analytical(t: float, x: NDArray[np.float64], neq: int) -> NDArray[np.float64]:
    nnodes = x.shape[0]
    u = np.zeros(nnodes*neq)
    xm, ym = x[:,0], x[:,1]
    u = np.sin(np.pi*xm)*np.sin(np.pi*ym)*np.cos(np.pi*t) + xm*ym*(t**2)
    return u

def main():

    neq = 1

    x0 = np.array([0.0, 0.0])
    xf = np.array([1.0, 1.0])
    mesh = mesh_discretization.MeshDiscretizationFD(x0=x0, xf=xf)
    nodes_dim = np.array([5,5])
    p = np.array([1,1])
    mesh.set_rectangular_mesh_fd(nodes_dim=nodes_dim, p=p, major_order='row')

    f_bc = [bc_x0, bc_xf, bc_y0, bc_yf]
    bc_type = [['dirichlet'], ['dirichlet'], ['dirichlet'], ['dirichlet']]
    boundary = boundary_conditions.BoundaryConditions(f_bc, bc_type)

    bvp_solver = bvp_setup.BVPSetupFD(neq, mesh, boundary, f_res)

    ls_solver = direct_solver.LUSolver()
    bvp_solver.set_ls_solver(ls_solver)

    t0, tf, dt, dt_min, dt_max, atol, rtol = 0.0, 0.5, 1e-2, 1e-2, 0.1, 1e-2, 1e-1
    # dt = tf

    u0 = initial_condition(t0, mesh.x_mesh, neq)

    bvp_solver.set_nr_solver(u0, k_max=100, tol=1e-8, r=1.0)

    bvp_solver.set_time_int(t0, tf, dt, dt_min, dt_max, atol, rtol, u0)

    bvp_solver.set_fd_problem(theta=1.0)

    u = bvp_solver.solve(dtw=tf)

    x = mesh.x_mesh
    u_exact = u_analytical(tf, x, neq)
    err = np.max(np.abs(u_exact-u))
    print(f'err = {err:.4e}')

    plot_pde.plot_contour(mesh, u, u_exact)

if __name__ == '__main__':
    main()
