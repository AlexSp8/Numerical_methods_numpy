
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

def f_weak(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    dudt: NDArray[np.float64], grad_u: NDArray[np.float64], w: float,
    grad_w: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)

    beta = 0.5
    pi = np.pi

    xm, ym = x[0], x[1]

    dudx, dudy = grad_u[0,:], grad_u[1,:]
    dwdx, dwdy = grad_w[0], grad_w[1]

    sin_x, sin_y = np.sin(pi*xm), np.sin(pi*ym)
    cos_x, cos_y = np.cos(pi*xm), np.cos(pi*ym)

    s = (2.0*pi**2*sin_x*sin_y*np.cos(pi*t) - beta*pi**2*cos_x*cos_y*np.cos(pi*t)
         - pi*sin_x*sin_y*np.sin(pi*t) + 2.0*xm*ym*t - beta*t**2)

    res[0] =  - (dudx[0]*dwdx + beta*dudx[0]*dwdy + dudy[0]*dwdy) + (s - dudt[0])*w

    return res

def f_strong(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
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

    res[0] =  + (d2udx2[0] + beta*d2udxdy[0] + d2udy2[0]) + s - dudt[0]

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
    nnodes = x.shape[0]
    u = np.zeros(nnodes*neq)
    xm, ym = x[:,0], x[:,1]
    u = np.sin(np.pi*xm)*np.sin(np.pi*ym)*np.cos(np.pi*t) + xm*ym*(t**2)
    return u

def main():

    neq = 1
    nbf = 4
    fe_type = 'quadrangle'
    # nbf = 6
    # fe_type = 'triangle'
    x0 = np.array([0.0, 0.0])
    xf = np.array([1.0, 1.0])
    mesh = mesh_discretization.MeshDiscretizationFE(x0=x0, xf=xf, nbf=nbf, fe_type=fe_type)
    nel_dim = np.array([5,5])
    p = np.array([1,1])
    mesh.set_rectangular_mesh_fe(nel_dim=nel_dim, p=p, major_order='row')

    f_bc = [bc_x0, bc_xf, bc_y0, bc_yf]
    bc_type = [['dirichlet'], ['dirichlet'], ['dirichlet'], ['dirichlet']]
    boundary = boundary_conditions.BoundaryConditions(f_bc, bc_type)

    bvp_solver = bvp_setup.BVPSetupFE(neq, mesh, boundary, f_weak)

    ls_solver = direct_solver.LUSolver()
    bvp_solver.set_ls_solver(ls_solver)

    t0, tf, dt, dt_min, dt_max, atol, rtol = 0.0, 0.1, 1e-2, 1e-2, 0.1, 1e-2, 1e-1
    # dt = tf

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
