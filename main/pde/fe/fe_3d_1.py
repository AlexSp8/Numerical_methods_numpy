
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

def bc_x0(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    grad_u: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)
    xm, ym, zm = x[0], x[1], x[2]
    v = t*(xm**2 + ym**2 + zm**2)
    res[0] = u[0] - v
    return res

def bc_xf(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    grad_u: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)
    xm, ym, zm = x[0], x[1], x[2]
    v = t*(xm**2 + ym**2 + zm**2)
    res[0] = u[0] - v
    return res

def bc_y0(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    grad_u: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)
    xm, ym, zm = x[0], x[1], x[2]
    v = t*(xm**2 + ym**2 + zm**2)
    res[0] = u[0] - v
    return res

def bc_yf(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    grad_u: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)
    xm, ym, zm = x[0], x[1], x[2]
    v = t*(xm**2 + ym**2 + zm**2)
    res[0] = u[0] - v
    return res

def bc_z0(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    grad_u: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)
    xm, ym, zm = x[0], x[1], x[2]
    v = t*(xm**2 + ym**2 + zm**2)
    res[0] = u[0] - v
    return res

def bc_zf(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    grad_u: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)
    xm, ym, zm = x[0], x[1], x[2]
    v = t*(xm**2 + ym**2 + zm**2)
    res[0] = u[0] - v
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

    xm, ym, zm = x[0], x[1], x[2]

    dudx = grad_u[0,:]
    dudy = grad_u[1,:]
    dudz = grad_u[2,:]

    dwdx = grad_w[0]
    dwdy = grad_w[1]
    dwdz = grad_w[2]

    s = (xm**2+ym**2+zm**2) - 6*t

    res[0] = dudt[0]*w + (dudx[0]*dwdx + dudy[0]*dwdy + dudz[0]*dwdz) - s*w
    return res

def f_strong(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    dudt: NDArray[np.float64], grad_u: NDArray[np.float64],
    hess_u: NDArray[np.float64]) -> NDArray[np.float64]:

    neq = u.shape[0]
    res = np.zeros(neq)

    xm, ym, zm = x[0], x[1], x[2]

    d2udx2 = hess_u[0,0,:]
    d2udy2 = hess_u[1,1,:]
    d2udz2 = hess_u[2,2,:]

    s = (xm**2+ym**2+zm**2) - 6*t

    res[0] = dudt[0] - (d2udx2[0] + d2udy2[0] + d2udz2[0]) - s

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
    total_nodes = x.shape[0]
    u_exact = np.zeros(total_nodes*neq)
    xm, ym, zm = x[:,0], x[:,1], x[:,2]
    u_exact[0::neq] = t*(xm**2+ym**2+zm**2)
    return u_exact

def main():

    neq = 1
    nbf = 8
    fe_type = 'hexahedron'
    # nbf = 10
    # fe_type = 'tetrahedron'
    x0 = np.array([0.0, 0.0, 0.0])
    xf = np.array([1.0, 1.0, 1.0])
    mesh = mesh_discretization.MeshDiscretizationFE(x0=x0, xf=xf, nbf=nbf, fe_type=fe_type)
    nel_dim = np.array([1,1,1])
    p = np.array([1,1,1])
    mesh.set_rectangular_mesh_fe(nel_dim=nel_dim, p=p, major_order='row')

    f_bc = [bc_x0, bc_xf, bc_y0, bc_yf, bc_z0, bc_zf]
    bc_type = [['dirichlet'], ['dirichlet'], ['dirichlet'], ['dirichlet'],
               ['dirichlet'], ['dirichlet']]
    boundary = boundary_conditions.BoundaryConditions(f_bc, bc_type)

    bvp_solver = bvp_setup.BVPSetupFE(neq, mesh, boundary, f_weak)

    ls_solver = direct_solver.LUSolver()
    # ls_solver = direct_solver.SparseSolver()
    bvp_solver.set_ls_solver(ls_solver)

    t0, tf, dt, dt_min, dt_max, atol, rtol = 0.0, 0.1, 1e-2, 1e-2, 0.1, 1e-2, 1e-1
    # dt = tf

    u0 = initial_condition(t0, mesh.x_mesh, neq)

    bvp_solver.set_nr_solver(u0, k_max=100, tol=1e-8, r=1.0)

    bvp_solver.set_time_int(t0, tf, dt, dt_min, dt_max, atol, rtol, u0)

    if nbf == 4 and fe_type == 'tetrahedron':
        ng = 4
    elif nbf == 10 and fe_type == 'tetrahedron':
        ng = 10
    elif nbf == 8 and fe_type == 'hexahedron':
        ng = 2
    elif (nbf == 20 or nbf == 27) and fe_type == 'hexahedron':
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

    for i in range(neq):
        ue = u_exact[i::neq]
        un = u[i::neq]
        plot_pde.plot_contour(mesh, un, ue)
        plot_pde.plot_3d_volume(mesh, un, ue)

if __name__ == '__main__':
    main()
