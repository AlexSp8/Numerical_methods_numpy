
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

# Convection-diffusion reaction
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

    v, D, s = physical_quantities(t, x, u)

    dudx = grad_u[0,:]
    dwdx = grad_w[0]

    # b = 1.0/np.tanh(Pe) - 1.0/Pe # optimal b (0 <= b <= 1)
    # D += b*v*h/2.0 # upwind added diffusion

    # w_mod = w + b*h*dwdx/2.0 # or modified weighting function for convection only (SU)
    # res[0] = v*dudx[0]*w_mod - s*w + D*dudx[0]*dwdx

    res[0] = v[0]*dudx[0]*w - s[0]*w + D[0]*dudx[0]*dwdx

    return res

def f_strong(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    dudt: NDArray[np.float64], grad_u: NDArray[np.float64],
    hess_u: NDArray[np.float64]) -> NDArray[np.float64]:

    neq = u.shape[0]
    res = np.zeros(neq)

    v, D, s = physical_quantities(t, x, u)

    dudx = grad_u[0,:]
    d2udx2 = hess_u[0,0,:]

    res[0] = v[0]*dudx[0] - D[0]*(d2udx2[0] + s[0])
    return res

def physical_quantities(t: float, x: NDArray[np.float64], u: NDArray[np.float64]
    ) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:

    neq = u.shape[0]

    v = np.ones(neq)

    D = np.zeros(neq)
    # h = 0.1
    # Pe = 5.0
    # D[0] = v[0]*h/(2*Pe)
    D[0] = 0.01

    s = np.zeros(neq)
    s[0] = 1.0

    return v, D, s

def stabilization_operator(t: float, x: NDArray[np.float64], u: NDArray[np.float64],
    w: float, dwdt: float, grad_w: NDArray[np.float64],
    hess_w: NDArray[np.float64]) -> NDArray[np.float64]:

    neq = u.shape[0]
    p_mat = np.zeros((neq, neq))

    v, D, s = physical_quantities(t, x, u)

    dDdx = 0

    dwdx = grad_w[0]

    p_mat[0,0] = v[0]*dwdx # SUPG
    # p_mat[0,0] = v[0]*dwdx - (dDdx*dwdx + D[0]*d2wdx2) + s[0]*w # GLS
    # p_mat[0,0] = v[0]*dwdx + (dDdx*dwdx + D[0]*d2wdx2) - s[0]*w # SGS

    return p_mat

def u_analytical(t: float, x: NDArray[np.float64], neq: int) -> NDArray[np.float64]:

    nnodes = x.shape[0]
    u_exact = np.zeros(nnodes*neq)

    v, D, s = physical_quantities(t, x, np.zeros(neq))

    g = v[0]/D[0]

    xm = x[:,0]

    u_exact[0::neq] = (1.0/v[0])*(xm-(1.0-np.exp(g*xm))/(1.0-np.exp(g)))
    return u_exact

def main():

    neq = 1
    nbf = 2
    fe_type = 'line'
    x0 = np.array([0.0])
    xf = np.array([1.0])
    mesh = mesh_discretization.MeshDiscretizationFE(x0=x0, xf=xf, nbf=nbf, fe_type=fe_type)
    nel_dim = np.array([10])
    p = np.array([1])
    mesh.set_rectangular_mesh_fe(nel_dim=nel_dim, p=p, major_order='row')

    f_bc = [bc_x0, bc_xf]
    bc_type = [['dirichlet'], ['dirichlet']]
    boundary = boundary_conditions.BoundaryConditions(f_bc, bc_type)

    bvp_solver = bvp_setup.BVPSetupFE(neq, mesh, boundary, f_weak)

    ls_solver = direct_solver.LUSolver()
    bvp_solver.set_ls_solver(ls_solver)

    t0, tf, dt, dt_min, dt_max, atol, rtol = 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
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
        # ui_num = u[ieq::neq].copy()
        # plot.plot_ode_system_evolution(x, [ui_num])

if __name__ == '__main__':
    main()
