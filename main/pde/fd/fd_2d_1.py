
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

# Linear: Conflicting corner BCs (Dirichlet+Dirichlet, Dirichlet+Neumann)
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
    res[0] = u[0] - 0.0
    return res

def bc_y0(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    grad_u: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)
    dudy = grad_u[1, :]
    res[0] = dudy[0] - 0.0
    return res

def bc_yf(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    grad_u: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - 1.0
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
    d2udx2, d2udy2 = hess_u[0,0,:], hess_u[1,1,:]
    res[0] = d2udx2[0] + d2udy2[0]
    return res

def u_analytical(t: float, x: NDArray[np.float64], neq: int) -> NDArray[np.float64]:
    xi = x[:,0]
    yi = x[:,1]
    nnodes = xi.shape[0]
    u_exact = np.zeros(nnodes*neq)

    terms = 100
    for k in range(1, terms*2, 2):
        u_k = (1.0/k)*(np.cosh(k*np.pi*yi)/np.cosh(k*np.pi))*np.sin(k*np.pi*xi)
        u_exact += u_k

    u_exact *= 4.0/np.pi

    for i in range(nnodes):
        if np.abs(xi[i] - 0.0) < 1e-8 and np.abs(yi[i]-1.0) < 1e-8:
            r = i*neq
            u_exact[r:r+neq] = 0.5
        if np.abs(xi[i] - 1.0) < 1e-8 and np.abs(yi[i]-1.0) < 1e-8:
            r = i*neq
            u_exact[r:r+neq] = 0.5

    return u_exact

def main():

    neq = 1

    x0 = np.array([0.0, 0.0])
    xf = np.array([1.0, 1.0])
    mesh = mesh_discretization.MeshDiscretizationFD(x0=x0, xf=xf)
    nodes_dim = np.array([11, 11])
    p = np.array([1,1])
    mesh.set_rectangular_mesh_fd(nodes_dim=nodes_dim, p=p, major_order='row')

    # bc = [[bc_x0, bc_y0], [bc_xf, bc_yf]]
    f_bc = [bc_x0, bc_xf, bc_y0, bc_yf]
    bc_type = [['dirichlet'], ['dirichlet'], ['neumann'], ['dirichlet']]
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

    plot_pde.plot_contour(mesh, u, u_exact)

if __name__ == '__main__':
    main()
