
import sys
from pathlib import Path

import numpy as np

current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent.parent

if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from linear_systems import direct_solver
from ode import plot
from pde import bvp_setup

# Heat Transfer
def bc_left(t: float, x_b: float, u_b: np.ndarray[tuple[int]],
    dudx_b: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    neq = u_b.shape[0]
    res = np.zeros(neq)
    res[0] = u_b[0] - 40.0
    return res

def bc_right(t: float, x_b: float, u_b: np.ndarray[tuple[int]],
    dudx_b: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    neq = u_b.shape[0]
    res = np.zeros(neq)
    res[0] = u_b[0]- 200.0
    # res[0] = dudx_b[0] - 0.0
    return res

def initial_condition(t0: float, x: np.ndarray[tuple[int]], neq: int
    ) -> np.ndarray[tuple[int]]:

    nnodes = x.shape[0]

    u0 = np.zeros(nnodes*neq)

    u_left = np.zeros(neq)
    u_right = np.zeros(neq)

    ieq = 0
    u0[ieq::neq] = 0.0 # x

    u_left[ieq] = 40.0
    u_right[ieq] = 200.0

    u0[:neq] = u_left
    u0[-neq:] = u_right

    return u0

def f_weak(t: float, x: float, u: np.ndarray[tuple[int]], dudt: np.ndarray[tuple[int]],
    dudx: np.ndarray[tuple[int]], w: float, dwdx: float) -> np.ndarray[tuple[int]]:
    T = u[0]
    h, Ta = 0.01, 20.0
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = dudt[0]*w + dudx[0]*dwdx - h*(Ta-T)*w
    return res

def f_strong(t: float, x: float, u: np.ndarray[tuple[int]], dudt: np.ndarray[tuple[int]],
    dudx: np.ndarray[tuple[int]], d2udx2: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    T = u[0]
    h, Ta = 0.01, 20.0
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = dudt[0] - (d2udx2[0] + h*(Ta-T))
    return res

def physical_quantities(t: float, x: float, u: np.ndarray[tuple[int]]):

    neq = u.shape[0]
    D = np.ones(neq)

    v = np.zeros(neq)

    s = np.zeros(neq)

    return v, D, s

def stabilization_operator(t: float, x: float, u: np.ndarray[tuple[int]],
    w: float, dwdt: np.ndarray[tuple[int]], dwdx: float,
    d2wdx2: float) -> np.ndarray[tuple[int]]:

    neq = u.shape[0]
    p_mat = np.zeros((neq, neq))

    # v, D, s = physical_quantities(t, x, u)

    # dDdx = 0

    # p_mat[0,0] = v*dwdx # SUPG
    # p_mat[0,0] = dwdt + v*dwdx - (dDdx*dwdx + D*d2wdx2) + s*w # GLS
    # p_mat[0,0] = dwdt + v*dwdx + (dDdx*dwdx + D*d2wdx2) - s*w # SGS

    return p_mat

def u_analytical(t: float, x: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:

    Ta, L = 20.0, 10.0
    c = (180-20*np.cosh(0.1*L))/np.sinh(0.1*L)
    T_ss = Ta*np.cosh(0.1*x) + c*np.sinh(0.1*x) + Ta
    n_tot = 5
    L2 = L**2
    T_tot = T_ss
    for n in range(1, n_tot):
        n_p = np.pi*n
        c1 = 2.0/(n_p*(0.01*L2+n_p**2))
        c2 = ((-1)**n)*(0.2*L2+200*(n_p)**2) - (0.2*L2+40*n_p)
        an = c1*c2
        T_tot += an*np.sin(n_p*x/L)*np.exp(-(n_p**2+1)*t/L2)
    return T_tot

def main():

    neq = 1

    bvp_solver = bvp_setup.BVPSetupFE(neq)

    xd = np.array([0.0, 10.0])
    nel, nbf = 5, 3
    nnodes = (nbf-1)*nel+1
    bvp_solver.create_mesh(nnodes=nnodes, xd=xd, p=1.0)
    x = bvp_solver.x

    D, v = 1.0, 0.0

    bvp_solver.check_dx(D, v)

    ls_solver = direct_solver.LUSolver()
    bvp_solver.set_ls_solver(ls_solver)

    t0, tf, dt, dt_min, dt_max, atol, rtol = 0.0, 250.0, 1e-1, 1e-1, 10.0, 1e-3, 1e-2
    # dt = bvp_solver.check_dt(D, v, dt)
    
    u0 = initial_condition(t0, x, neq)

    bvp_solver.set_nr_solver(u0, k_max=100, tol=1e-8, r=1.0)

    bvp_solver.set_time_int(t0, tf, dt, dt_min, dt_max, atol, rtol, u0)
    
    ng = 3
    bvp_solver.set_fe_gauss_int(ng=ng, nbf=nbf)
    
    bvp_solver.set_fe_stabilization(f_quantities=physical_quantities, f_strong=f_strong,
                                    f_operator=stabilization_operator, order=2)

    bc = {'left': bc_left, 'right': bc_right}

    bvp_solver.set_fe_problem(f_weak, bc, theta=1.0)

    u = bvp_solver.solve(dtw=100.0)

    u_exact = u_analytical(tf, x)
    err = np.max(np.abs(u_exact-u))
    print(f'err = {err:.4e}')

    for ieq in range(neq):
        ui_num, ui_ex = u[ieq::neq].copy(), u_exact[ieq::neq].copy()
        plot.plot_ode_system_evolution(x, [ui_num, ui_ex])

if __name__ == '__main__':
    main()
