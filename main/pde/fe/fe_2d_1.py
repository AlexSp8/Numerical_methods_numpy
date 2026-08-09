
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
    res[0] = - 0.0 # dudy = 0
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

def f_weak(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    dudt: NDArray[np.float64], grad_u: NDArray[np.float64], w: float,
    grad_w: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)
    dudx, dudy = grad_u[0,:], grad_u[1,:]
    dwdx, dwdy = grad_w[0], grad_w[1]
    res[0] = dudx[0]*dwdx + dudy[0]*dwdy
    return res

def f_strong(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    dudt: NDArray[np.float64], grad_u: NDArray[np.float64],
    hess_u: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)
    d2udx2, d2udy2 = hess_u[0,0,:], hess_u[1,1,:]
    res[0] = d2udx2[0] + d2udy2[0]
    return res

def physical_quantities(t: float, x: NDArray[np.float64], u: NDArray[np.float64]
    ) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    neq = u.shape[0]
    D = np.ones(neq)
    v = np.zeros(neq)
    s = np.zeros(neq)
    return v, D, s

def stabilization_operator(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    w: float, dwdt: float, grad_w: NDArray[np.float64],
    hess_w: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    p_mat = np.zeros((neq, neq))
    return p_mat

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
    nbf = 4
    fe_type = 'quadrangle'
    # nbf = 3
    # fe_type = 'triangle'
    x0 = np.array([0.0, 0.0])
    xf = np.array([1.0, 1.0])
    mesh = mesh_discretization.MeshDiscretizationFE(x0=x0, xf=xf, nbf=nbf, fe_type=fe_type)
    nel_dim = np.array([5,5])
    p = np.array([1,1])
    mesh.set_rectangular_mesh_fe(nel_dim=nel_dim, p=p, major_order='row')

    # for i in range(mesh.nel):
    #     gnodes = mesh.connectivity[i,:]
    #     print('-------------------')
    #     print(i, gnodes)
    #     for gnod in gnodes:
    #         print(gnod, mesh.x_mesh[gnod,:])

    # for i in mesh.boundary_nodes:
    #     print(i)
    # for i in mesh.multi_boundary_nodes:
    #     print(i)

    f_bc = [bc_x0, bc_xf, bc_y0, bc_yf]
    bc_type = [['dirichlet'], ['dirichlet'], ['neumann'], ['dirichlet']]
    boundary = boundary_conditions.BoundaryConditions(f_bc, bc_type)

    bvp_solver = bvp_setup.BVPSetupFE(neq, mesh, boundary, f_weak)

    ls_solver = direct_solver.LUSolver()
    bvp_solver.set_ls_solver(ls_solver)

    t0, tf, dt, dt_min, dt_max, atol, rtol = 0.0, 1.0, 1.0, 1.0, 1.0, 1e-3, 1e-2
    dt = tf

    u0 = initial_condition(t0, mesh.x_mesh, neq)

    bvp_solver.set_nr_solver(u0, k_max=100, tol=1e-8, r=1.0)

    bvp_solver.set_time_int(t0, tf, dt, dt_min, dt_max, atol, rtol, u0)

    if nbf == 3 and fe_type == 'triangle':
        ng = 4
    elif nbf == 6 and fe_type == 'triangle':
        ng = 6
    elif nbf == 4 and fe_type == 'quadrangle':
        ng = 2
    elif (nbf == 8 or nbf == 9) and fe_type == 'quadrangle':
        ng = 3
    bvp_solver.set_fe_gauss_int(ng=ng)

    bvp_solver.set_fe_stabilization(f_quantities=physical_quantities, f_strong=f_strong,
                                    f_operator=stabilization_operator)

    bvp_solver.set_fe_problem(theta=1.0)

    u = bvp_solver.solve(dtw=tf)

    x = mesh.x_mesh
    u_exact = u_analytical(tf, x, neq)
    err = np.max(np.abs(u_exact-u))
    print(f'err = {err:.4e}')

    plot_pde.plot_contour(mesh, u, u_exact)

if __name__ == '__main__':
    main()
