
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

# Non-Linear
def bc_x0(t: float, x: float, u: NDArray[np.float64],
    grad_u: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - 1.0
    return res

def bc_xf(t: float, x: float, u: NDArray[np.float64],
    grad_u: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = - u[0]**2 + 21.0 # dudx = -u**2 + 21
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
    dudx = grad_u[0,:]
    dwdx = grad_w[0]
    xm = x[0]
    s = xm**6 - 5.0*(xm**4) + 2.0*(xm**3) - 8.0*(xm**2) + 10.0*xm - 3.0
    res[0] = - dudx[0]*dwdx + (- (dudx[0])**2 + u[0]**2 - s)*w
    return res

def f_strong(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    dudt: NDArray[np.float64], grad_u: NDArray[np.float64],
    hess_u: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u.shape[0]
    res = np.zeros(neq)
    d2udx2 = hess_u[0,0,:]
    dudx = grad_u[0,:]
    xm = x[0]
    s = xm**6 - 5.0*(xm**4) + 2.0*(xm**3) - 8.0*(xm**2) + 10.0*xm - 3.0
    res[0] = d2udx2[0] - (dudx[0])**2 + u[0]**2 - s
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
    xm = x[:,0]
    nnodes = xm.shape[0]
    u = np.zeros(nnodes*neq)
    u = xm**3 + 2.0*xm + 1.0
    return u

def main():

    neq = 1
    nbf = 3
    fe_type = 'line'
    x0 = np.array([0.0])
    xf = np.array([1.0])
    mesh = mesh_discretization.MeshDiscretizationFE(x0=x0, xf=xf, nbf=nbf, fe_type=fe_type)
    nel_dim = np.array([10])
    p = np.array([1])
    mesh.set_rectangular_mesh_fe(nel_dim=nel_dim, p=p, major_order='row')

    f_bc = [bc_x0, bc_xf]
    bc_type = [['dirichlet'], ['robin']]
    boundary = boundary_conditions.BoundaryConditions(f_bc, bc_type)

    bvp_solver = bvp_setup.BVPSetupFE(neq, mesh, boundary, f_weak)

    ls_solver = direct_solver.LUSolver()
    bvp_solver.set_ls_solver(ls_solver)

    t0, tf, dt, dt_min, dt_max, atol, rtol = 0.0, 1.0, 1e-1, 1e-1, 10.0, 1e-3, 1e-2
    dt = tf

    u0 = initial_condition(t0, mesh.x_mesh, neq)

    bvp_solver.set_nr_solver(u0, k_max=100, tol=1e-8, r=1.0)

    bvp_solver.set_time_int(t0, tf, dt, dt_min, dt_max, atol, rtol, u0)

    if nbf == 2:
        ng = 2
    else:
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

    for ieq in range(neq):
        ui_num, ui_ex = u[ieq::neq].copy(), u_exact[ieq::neq].copy()
        plot.plot_ode_system_evolution(x, [ui_num, ui_ex])

if __name__ == '__main__':
    main()
