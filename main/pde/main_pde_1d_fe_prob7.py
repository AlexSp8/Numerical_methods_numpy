
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

# Problem 7: convection-diffusion reaction
def bc_left(t: float, x: float, u: np.ndarray[tuple[int]], dudx: np.ndarray[tuple[int]]
    ) -> np.ndarray[tuple[int]]:

    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - 0.0

    return res

def bc_right(t: float, x: float, u: np.ndarray[tuple[int]], dudx: np.ndarray[tuple[int]]
    ) -> np.ndarray[tuple[int]]:

    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - 1.0

    return res

def initial_condition(t0: float, x: np.ndarray[tuple[int]], neq: int
    ) -> np.ndarray[tuple[int]]:

    nnodes = x.shape[0]
    u0 = np.zeros(nnodes*neq)

    return u0

def f_weak(t: float, x: float, u: np.ndarray[tuple[int]], dudt: np.ndarray[tuple[int]],
    dudx: np.ndarray[tuple[int]], w: float, dwdx: float) -> np.ndarray[tuple[int]]:

    neq = u.shape[0]
    res = np.zeros(neq)

    v, D, s = physical_quantities(t, x, u)

    # b = 1.0/np.tanh(Pe) - 1.0/Pe # optimal b (0 <= b <= 1)
    # D += b*v*h/2.0 # upwind added diffusion

    # w_mod = w + b*h*dwdx/2.0 # or modified weighting function for convection only (SU)
    # res[0] = v*dudx[0]*w_mod - s*w + D*dudx[0]*dwdx

    res[0] = v[0]*dudx[0]*w - s[0]*w + D[0]*dudx[0]*dwdx

    return res

def f_strong(t: float, x: float, u: np.ndarray[tuple[int]], dudt: np.ndarray[tuple[int]],
    dudx: np.ndarray[tuple[int]], d2udx2: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:

    neq = u.shape[0]
    res = np.zeros(neq)

    v, D, s = physical_quantities(t, x, u)

    res[0] = v[0]*dudx[0] - D[0]*(d2udx2[0] + s[0])
    return res

def physical_quantities(t: float, x: float, u: np.ndarray[tuple[int]]):

    neq = u.shape[0]

    v = np.ones(neq)

    D = np.ones(neq)
    h = 0.1
    Pe = 5.0
    D[0] = v[0]*h/(2*Pe)
    # D[0] = 0.01

    s = np.zeros(neq)
    s[0] = 10.0*np.exp(-5.0*x) - 4.0*np.exp(-x)

    return v, D, s

def stabilization_operator(t: float, x: float, u: np.ndarray[tuple[int]],
    w: float, dwdt: np.ndarray[tuple[int]], dwdx: float,
    d2wdx2: float) -> np.ndarray[tuple[int]]:

    neq = u.shape[0]
    p_mat = np.zeros((neq, neq))

    v, D, s = physical_quantities(t, x, u)

    dDdx = 0

    p_mat[0,0] = v[0]*dwdx # SUPG
    # p_mat[0,0] = v[0]*dwdx - (dDdx*dwdx + D[0]*d2wdx2) + s[0]*w # GLS
    # p_mat[0,0] = v[0]*dwdx + (dDdx*dwdx + D[0]*d2wdx2) - s[0]*w # SGS

    return p_mat

def u_analytical(t: float, x: np.ndarray) -> np.ndarray:
    
    v, D, _ = physical_quantities(t, x[0], np.array([0.0]))
    L = 1.0
    
    alpha = v/D
    A = (2.0*D)/(5.0*D + v)
    B = (4.0*D)/(D + v)
    
    if abs(5.0*D + v) < 1e-12 or abs(D + v) < 1e-12:
        raise ValueError("Physics parameters trigger a mathematical resonance condition.")
        
    numerator = 1.0 - A*(1.0 - np.exp(-5.0*L)) + B*(1.0 - np.exp(-L))
    C2 = numerator / (np.exp(alpha*L) - 1.0)
    C1 = A - B - C2
    
    u_exact = C1 + C2*np.exp(alpha*x) - A*np.exp(-5.0*x) + B*np.exp(-x)

    return u_exact
    
def main():

    neq = 1

    bvp_solver = bvp_setup.BVPSetupFE(neq)

    xd = np.array([0.0, 1.0])
    nel, nbf = 10, 2
    nnodes = (nbf-1)*nel+1
    bvp_solver.create_mesh(nnodes=nnodes, xd=xd, p=1.0)
    x = bvp_solver.x

    # D, v = 1.0, 0.0

    # bvp_solver.check_dx(D, v)

    ls_solver = direct_solver.LUSolver()
    bvp_solver.set_ls_solver(ls_solver)

    t0, tf, dt, dt_min, dt_max, atol, rtol = 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
    # dt = bvp_solver.check_dt(D, v, dt)
    
    u0 = initial_condition(t0, x, neq)

    bvp_solver.set_nr_solver(u0, k_max=100, tol=1e-8, r=1.0)

    bvp_solver.set_time_int(t0, tf, dt, dt_min, dt_max, atol, rtol, u0)
    
    ng = 2
    bvp_solver.set_fe_gauss_int(ng=ng, nbf=nbf)

    bvp_solver.set_fe_stabilization(f_quantities=physical_quantities, f_strong=f_strong,
                                    f_operator=stabilization_operator, order=1)

    bc = {'left': bc_left, 'right': bc_right}

    bvp_solver.set_fe_problem(f_weak, bc, theta=1.0)

    u = bvp_solver.solve(dtw=0.2)

    u_exact = u_analytical(tf, x)
    err = np.max(np.abs(u_exact-u))
    print(f'err = {err:.4e}')

    for ieq in range(neq):
        ui_num, ui_ex = u[ieq::neq].copy(), u_exact[ieq::neq].copy()
        plot.plot_ode_system_evolution(x, [ui_num, ui_ex])
        # ui_num = u[ieq::neq].copy()
        # plot.plot_ode_system_evolution(x, [ui_num])

if __name__ == '__main__':
    main()
